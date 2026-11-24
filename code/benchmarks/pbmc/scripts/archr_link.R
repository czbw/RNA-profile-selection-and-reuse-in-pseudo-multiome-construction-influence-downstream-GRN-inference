#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ArchR); library(Seurat); library(SummarizedExperiment)
  library(Matrix); library(rhdf5); library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(k, d = NULL) {
  i <- which(args == paste0("--", k)); if (length(i)) args[i + 1] else d
}
DS   <- getarg("dataset", "3k")
ROOT <- getarg("root", Sys.getenv("AMLGRN_PROJECT", getwd()))
THREADS <- as.integer(getarg("threads", "6"))
SEED <- 555

DATA   <- file.path(ROOT, "data/pbmc")
OUTDIR <- file.path(ROOT, "benchmarks/pbmc", "work", DS, "archr")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
setwd(OUTDIR)

FRAG <- file.path(DATA, sprintf("pbmc_granulocyte_sorted_%s_atac_fragments.tsv.gz", DS))
H5   <- file.path(DATA, sprintf("pbmc_granulocyte_sorted_%s_filtered_feature_bc_matrix.h5", DS))
stopifnot(file.exists(FRAG), file.exists(H5))

set.seed(SEED)
addArchRThreads(threads = THREADS)
addArchRGenome("hg38")
cat(sprintf(">> dataset=%s threads=%d\n>> frag=%s\n", DS, THREADS, FRAG))

# 1. ATAC side -- createArrowFiles defaults of the original AML study
existing <- file.path(OUTDIR, sprintf("pbmc_%s.arrow", DS))
if (file.exists(existing) && file.info(existing)$size > 1e6) {
  cat(">> reusing existing arrow file (delete it to force a rebuild)
")
  ArrowFiles <- existing
} else {
  ArrowFiles <- createArrowFiles(
    inputFiles = setNames(FRAG, sprintf("pbmc_%s", DS)),
    sampleNames = sprintf("pbmc_%s", DS),
    minTSS = 4, minFrags = 2000,
    addTileMat = TRUE, addGeneScoreMat = TRUE, force = TRUE)
}
cat(">> arrow files:", ArrowFiles, "\n")

proj <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = file.path(OUTDIR, "proj"),
                     copyArrows = FALSE)
cat(sprintf(">> ArchR project: %d cells\n", nCells(proj)))

proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI",
                        force = TRUE, seed = SEED)
proj <- addImputeWeights(proj)          # the MAGIC imputation step of the original study's ATAC preprocessing

# 2. RNA side -- read the GEX half of the same 10x h5, then discard the barcodes
# Read10X_h5 returns a list when the file holds several feature types (multiome does).
tenx <- Seurat::Read10X_h5(H5)
gex  <- if (is.list(tenx)) tenx[["Gene Expression"]] else tenx
cat(sprintf(">> RNA matrix: %d genes x %d cells\n", nrow(gex), ncol(gex)))

true_bc <- colnames(gex)                                   # held back for scoring only
# Rename the RNA cells so ArchR cannot match on barcode identity (otherwise the integration
# could trivially recover the native reference).
colnames(gex) <- paste0("RNACELL", seq_len(ncol(gex)))
rna_map <- data.frame(rna_id = colnames(gex), true_barcode = true_bc,
                      stringsAsFactors = FALSE)

seu <- CreateSeuratObject(counts = gex, project = sprintf("pbmc_%s", DS))
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, verbose = FALSE)
seu$orig.ident <- sprintf("pbmc_%s", DS)                   # groupRNA, as in the original study

if (utils::packageVersion("Seurat") >= "5.0.0") {
  if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
    seu[["RNA"]] <- SeuratObject::JoinLayers(seu[["RNA"]])
    cat(">> Seurat 5 detected: layers joined
")
  }
}
seRNA <- seu

proj <- addGeneIntegrationMatrix(
  ArchRProj = proj, useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix", reducedDims = "IterativeLSI",
  seRNA = seRNA, addToArrow = FALSE, groupRNA = "orig.ident", force = TRUE)

cd <- as.data.frame(getCellColData(proj))
stopifnot("predictedCell" %in% colnames(cd))

# ArchR prefixes ATAC cells "<sample>#<barcode>"; recover the raw barcode for the join
atac_bc <- sub("^.*#", "", rownames(cd))
out <- data.frame(
  cistopic_id = paste0(atac_bc, "___", sprintf("pbmc_%s", DS)),
  barcode = atac_bc,
  sample_id = sprintf("pbmc_%s", DS),
  Linked_scRNA_barcode = rna_map$true_barcode[match(cd$predictedCell, rna_map$rna_id)],
  Linked_scRNA_Sample = sprintf("pbmc_%s", DS),
  Linked_scRNA_confidence_score = cd$predictedScore,
  predictedCell_internal = cd$predictedCell,
  stringsAsFactors = FALSE)
# the native reference, held back until now and used only for scoring in p_03
out$true_barcode <- atac_bc
out$link_is_correct <- out$Linked_scRNA_barcode == out$true_barcode

write.csv(out, file.path(OUTDIR, sprintf("link_%s.csv", DS)), row.names = FALSE)

keep <- !is.na(out$Linked_scRNA_confidence_score) & out$Linked_scRNA_confidence_score >= 0.5
acc_all  <- mean(out$link_is_correct, na.rm = TRUE)
acc_conf <- mean(out$link_is_correct[keep], na.rm = TRUE)
reuse_all  <- length(out$Linked_scRNA_barcode) / length(unique(out$Linked_scRNA_barcode))
reuse_conf <- sum(keep) / length(unique(out$Linked_scRNA_barcode[keep]))

rep <- list(
  dataset = DS, seed = SEED,
  versions = list(R = paste(R.version$major, R.version$minor, sep = "."),
                  ArchR = as.character(packageVersion("ArchR")),
                  Seurat = as.character(packageVersion("Seurat")),
                  lambo_pinned = "R 4.0.2 / Seurat 4.0.1 / ArchR 1.0.2 -- version drift declared",
                  addToArrow = "FALSE (Lambo used TRUE); only the discarded imputed GeneIntegrationMatrix differs, predictedCell/predictedScore identical",
                  seurat_major = as.integer(substr(as.character(packageVersion("Seurat")),1,1))),
  n_atac_cells = nrow(out),
  n_rna_cells = nrow(rna_map),
  conf_median = median(out$Linked_scRNA_confidence_score, na.rm = TRUE),
  frac_conf_ge_0.5 = mean(keep),
  accuracy_all = acc_all,
  accuracy_conf_ge_0.5 = acc_conf,
  distinct_rna_used_all = length(unique(out$Linked_scRNA_barcode)),
  distinct_rna_used_conf = length(unique(out$Linked_scRNA_barcode[keep])),
  reuse_factor_all = reuse_all,
  reuse_factor_conf = reuse_conf,
  aml_comparator = list(reuse_dx = 15.97, reuse_rel = 14.59,
                        note = "AML link reuse. TRUE pairing is 1.0x by definition."))
write_json(rep, file.path(OUTDIR, sprintf("link_report_%s.json", DS)),
           auto_unbox = TRUE, pretty = TRUE, digits = 6)

cat("\n==== LINK ARM SUMMARY ====\n")
cat(sprintf("  ATAC cells                : %d\n", nrow(out)))
cat(sprintf("  confidence >= 0.5         : %.1f%% (median %.3f)\n", 100*mean(keep),
            median(out$Linked_scRNA_confidence_score, na.rm = TRUE)))
cat(sprintf("  ACCURACY vs true barcode  : %.4f (all) | %.4f (conf>=0.5)\n", acc_all, acc_conf))
cat(sprintf("  distinct RNA cells used   : %d of %d\n",
            length(unique(out$Linked_scRNA_barcode)), nrow(rna_map)))
cat(sprintf("  REUSE FACTOR              : %.2fx (all) | %.2fx (conf>=0.5)   [truth = 1.00x]\n",
            reuse_all, reuse_conf))
cat(sprintf("  AML comparator                       : 15.97x Dx / 14.59x Rel\n"))
cat("==== STEP 02 COMPLETE ====\n")
