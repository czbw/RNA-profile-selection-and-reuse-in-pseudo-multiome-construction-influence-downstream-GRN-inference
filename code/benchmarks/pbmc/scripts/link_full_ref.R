#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ArchR); library(Seurat); library(SummarizedExperiment)
  library(Matrix); library(rhdf5); library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(k, d = NULL) { i <- which(args == paste0("--", k)); if (length(i)) args[i + 1] else d }
DS      <- getarg("dataset", "10k")
ROOT    <- getarg("root", Sys.getenv("AMLGRN_PROJECT", getwd()))
THREADS <- as.integer(getarg("threads", "6"))
SEED    <- 555

DATA   <- file.path(ROOT, "data/pbmc")
SRCDIR <- file.path(ROOT, "benchmarks/pbmc", "work", DS, "archr")
OUTDIR <- file.path(ROOT, "benchmarks/pbmc", "work", DS, "archr_fullref")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
setwd(OUTDIR)

H5 <- file.path(DATA, sprintf("pbmc_granulocyte_sorted_%s_filtered_feature_bc_matrix.h5", DS))
ARROW <- file.path(SRCDIR, sprintf("pbmc_%s.arrow", DS))
stopifnot(file.exists(H5), file.exists(ARROW))

set.seed(SEED)
addArchRThreads(threads = THREADS)
addArchRGenome("hg38")
cat(sprintf(">> p_16 dataset=%s threads=%d\n>> reusing arrow: %s\n", DS, THREADS, ARROW))

proj <- ArchRProject(ArrowFiles = ARROW, outputDirectory = file.path(OUTDIR, "proj"),
                     copyArrows = FALSE)
cat(sprintf(">> ArchR project: %d cells\n", nCells(proj)))
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI",
                        force = TRUE, seed = SEED)
proj <- addImputeWeights(proj)

# --- RNA side: same as p_02 (true barcodes discarded) -----------------------
tenx <- Seurat::Read10X_h5(H5)
gex  <- if (is.list(tenx)) tenx[["Gene Expression"]] else tenx
true_bc <- colnames(gex)
colnames(gex) <- paste0("RNACELL", seq_len(ncol(gex)))
rna_map <- data.frame(rna_id = colnames(gex), true_barcode = true_bc, stringsAsFactors = FALSE)
cat(sprintf(">> RNA matrix: %d genes x %d cells\n", nrow(gex), ncol(gex)))

seu <- CreateSeuratObject(counts = gex, project = sprintf("pbmc_%s", DS))
seu <- NormalizeData(seu, verbose = FALSE); seu <- FindVariableFeatures(seu, verbose = FALSE)
seu$orig.ident <- sprintf("pbmc_%s", DS)
if (utils::packageVersion("Seurat") >= "5.0.0" &&
    "JoinLayers" %in% getNamespaceExports("SeuratObject")) {
  seu[["RNA"]] <- SeuratObject::JoinLayers(seu[["RNA"]])
}

N_RNA <- ncol(seu)
N_ATAC <- nCells(proj)
cat(sprintf(">> sampleCellsRNA raised to %d (was ArchR default 10000); ATAC=%d\n", N_RNA, N_ATAC))
if (N_RNA < N_ATAC) {
  cat(sprintf("!! WARNING: %d RNA < %d ATAC -- a bijection is STILL impossible; %d cells must reuse\n",
              N_RNA, N_ATAC, N_ATAC - N_RNA))
}

proj <- addGeneIntegrationMatrix(
  ArchRProj = proj, useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix", reducedDims = "IterativeLSI",
  seRNA = seu, addToArrow = FALSE, groupRNA = "orig.ident",
  sampleCellsATAC = N_ATAC, sampleCellsRNA = N_RNA,   # <- the only change from p_02
  plotUMAP = TRUE, force = TRUE)

cd <- as.data.frame(getCellColData(proj))
stopifnot("predictedCell" %in% colnames(cd))
atac_bc <- sub("^.*#", "", rownames(cd))

out <- data.frame(
  cistopic_id = paste0(atac_bc, "___", sprintf("pbmc_%s", DS)),
  barcode = atac_bc, sample_id = sprintf("pbmc_%s", DS),
  Linked_scRNA_barcode = rna_map$true_barcode[match(cd$predictedCell, rna_map$rna_id)],
  Linked_scRNA_Sample = sprintf("pbmc_%s", DS),
  Linked_scRNA_confidence_score = cd$predictedScore,
  predictedCell_internal = cd$predictedCell,
  stringsAsFactors = FALSE)
out$true_barcode <- atac_bc
out$link_is_correct <- out$Linked_scRNA_barcode == out$true_barcode
write.csv(out, file.path(OUTDIR, sprintf("link_fullref_%s.csv", DS)), row.names = FALSE)

# --- export the joint CCA so the assignment can be solved in Python ---------
cca_fp <- file.path(OUTDIR, "proj", "RNAIntegration", "GeneIntegrationMatrix",
                    "Save-Block1-JointCCA.rds")
if (!file.exists(cca_fp)) {
  cand <- list.files(file.path(OUTDIR, "proj"), pattern = "JointCCA\\.rds$",
                     recursive = TRUE, full.names = TRUE)
  if (length(cand)) cca_fp <- cand[1]
}
if (file.exists(cca_fp)) {
  j <- as.data.frame(readRDS(cca_fp))
  keep <- c(grep("^CC_", colnames(j), value = TRUE),
            intersect(c("Assay", "Group", "Score"), colnames(j)))
  j <- j[, keep, drop = FALSE]
  j$cell <- rownames(j)
  write.csv(j, file.path(OUTDIR, sprintf("jointcca_%s.csv", DS)), row.names = FALSE)
  cat(sprintf(">> joint CCA exported: %d cells x %d cols -> jointcca_%s.csv\n",
              nrow(j), ncol(j) - 1, DS))
  n_rna_cca  <- sum(grepl("^RNACELL", rownames(j)))
  n_atac_cca <- sum(grepl("#", rownames(j)))
} else {
  cat("!! joint CCA rds not found -- p_17 cannot run\n"); n_rna_cca <- NA; n_atac_cca <- NA
}

reuse <- nrow(out) / length(unique(out$Linked_scRNA_barcode))
rep <- list(
  script = "link_full_ref.R", dataset = DS, seed = SEED,
  change = "sampleCellsRNA raised from ArchR default 10000 to all RNA cells; everything else identical to p_02",
  n_atac_cells = nrow(out), n_rna_cells = nrow(rna_map),
  n_rna_in_jointcca = n_rna_cca, n_atac_in_jointcca = n_atac_cca,
  bijection_possible = isTRUE(n_rna_cca >= n_atac_cca),
  accuracy_all = mean(out$link_is_correct, na.rm = TRUE),
  distinct_rna_used = length(unique(out$Linked_scRNA_barcode)),
  reuse_factor = reuse,
  link_original_comparator = list(accuracy_all = 0.004227, reuse_factor_all = 6.095395,
                                  distinct_rna_used = 1824, sampleCellsRNA = 10000))
write_json(rep, file.path(OUTDIR, sprintf("link_fullref_report_%s.json", DS)),
           auto_unbox = TRUE, pretty = TRUE, digits = 6)

cat("\n==== p_16 SUMMARY (greedy argmax, FULL reference) ====\n")
cat(sprintf("  ATAC cells        : %d\n", nrow(out)))
cat(sprintf("  ACCURACY          : %.4f   [link original: 0.004227]\n",
            mean(out$link_is_correct, na.rm = TRUE)))
cat(sprintf("  distinct RNA used : %d of %d   [link original: 1824 of 10000 sampled]\n",
            length(unique(out$Linked_scRNA_barcode)), nrow(rna_map)))
cat(sprintf("  REUSE             : %.2fx   [link original: 6.10x]\n", reuse))
cat("==== p_16 COMPLETE ====\n")
