#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2); library(grid); library(patchwork); library(scales)
})
source("figures/shared/theme_paper.R")     # sources palette.R
source("figures/shared/gridkit.R")
source("figures/shared/manifest.R")
source("figures/shared/number_labels.R")

OUT <- "figures/main/fig01_construction"
STAGES <- c("Diagnosis", "Relapse")
message("\n=== Figure 1 (v3): the two constructions, compared in three currencies ===")

PANEL_SO_WHAT <- list(
  a = "The published per-cell construction assigns one transcriptome to each chromatin profile, so one transcriptome can serve several and others serve none.",
  b = "The field-standard alternative never forms a per-cell assignment at all: it averages k cells on each side within a cell-type label and pairs the two averages.",
  c = "The metacell arm's regulator list is largely a subset of link's, and the containment is 20 points weaker once the sign mode is part of the regulator's identity.",
  d = "Over the identical set of shared regulons, the two constructions agree about where a regulator is active and disagree almost completely about what it regulates."
)

## THE SET'S COLOURS AND TYPE
INK    <- "#111111"
RED    <- "#D11F1F"                      # published link      == Figure 4
PURPLE <- "#5C3696"                      # metacell k = 5      == Figures 2 and 4
GREY   <- "#8C8C8C"                      # not an arm          == the set's reference grey
## schematic tints: 20% of the arm colour over white, so a highlighted cell box reads as
## belonging to its arm without becoming a solid data mark
tint <- function(hex, a = 0.20) {
  m <- grDevices::col2rgb(hex) / 255
  grDevices::rgb(1 - a * (1 - m[1]), 1 - a * (1 - m[2]), 1 - a * (1 - m[3]))
}
ARM_NAME <- c(link = "published link", pseudomultiome = "metacell k = 5")

TXT <- 1.18                              # +18%, the top of the owner's 15-20% range
ts  <- function(pt) pt * TXT             # theme and grob sizes, in points
gs  <- function(mm) mm * TXT             # geom_text sizes, in ggplot2's mm

theme_tight <- function() theme_paper() +
  theme(axis.title.x = element_text(size = ts(6.0), margin = margin(t = 1.5)),
        axis.title.y = element_text(size = ts(6.0), margin = margin(r = 1.5)),
        axis.text = element_text(size = ts(5.2)),
        legend.box.spacing = unit(1.5, "pt"), legend.margin = margin(0, 0, 0, 0),
        legend.background = element_blank(), legend.key.size = unit(ts(5.5), "pt"),
        legend.text = element_text(size = ts(5.2)))

## SHARED SCHEMATIC GRAMMAR -- v1 panel (a)'s geometry, unchanged
FS    <- ts(5.8)  # title size: the set's strip-title size, plain face
CELL  <- "#8C8C8C"                       # an ordinary cell's outline (neutral, not an arm)
FILL_N <- "#F7F7F7"                      # an ordinary cell
FILL_R <- tint(RED)                      # the reused transcriptome
FILL_A <- tint(PURPLE)                   # an aggregate
LWD_N  <- 0.6                            # ordinary box
LWD_E  <- 0.9                            # emphasised box / every arrow
DASH   <- "22"                           # v1's faded outline
SCH   <- 1.56
ROW_H  <- 0.050 * SCH
GAP    <- 0.013 * SCH
BOX_W  <- 0.205
BOX_W_B <- 0.170                         # panel (b) only -- four columns
ARR    <- 1.7                            # arrow-head length, mm
CURV   <- 0.20                           # curvature magnitude
N_RNA  <- 6L
N_ATAC <- 5L
BLOCK_H <- N_RNA * ROW_H + (N_RNA - 1) * GAP
Y0      <- 0.455 + BLOCK_H / 2 - ROW_H / 2

cellbox <- function(x, y, w = BOX_W, fill = FILL_N, col = CELL, lwd = LWD_N,
                    lty = "solid", h = ROW_H)
  rectGrob(x, y, w, h, gp = gpar(fill = fill, col = col, lwd = lwd, lty = lty))

## the two panel names. Plain face, as every other title in the set.
namelab <- function(txt, y = 0.955)
  textGrob(txt, 0.5, y, gp = gpar(fontsize = FS, fontface = "plain", col = INK))

## COLUMN LABELS (added 2026-09-15). Both schematics drew two columns of cell boxes and named
## neither, so which side held transcriptomes and which held chromatin profiles was recoverable
## only from the caption. These are the schematic's axis titles: one size BELOW the panel title,
## so the hierarchy reads title > column label, and in INK because grey encodes "not an arm"
## everywhere else in this set and must not be spent on a structural label.
FS_COL <- ts(5.2)                        # == the set's axis.text / legend.text size
LAB_Y  <- 0.828                          # between the block top (0.740) and the title (0.955)
collab <- function(txt, x, y = LAB_Y)
  textGrob(txt, x, y, gp = gpar(fontsize = FS_COL, fontface = "plain", col = INK))
RNA_LAB  <- "RNA cells"
ATAC_LAB <- "ATAC cells"

## (a) per-cell matching (published link)
panel_a <- function() {
  .g <- list(); add <- function(...) { .g[[length(.g) + 1L]] <<- gList(...); invisible(NULL) }

  ys_r <- Y0 - (seq_len(N_RNA)  - 1) * (ROW_H + GAP)
  ys_a <- Y0 - (seq_len(N_ATAC) - 1) * (ROW_H + GAP)
  rna_x <- 0.205; atac_x <- 0.690
  reused <- 2L
  unused <- c(4L, 5L, 6L)
  assign_map <- c(reused, reused, reused, 1L, 3L)

  for (i in seq_len(N_RNA)) {
    is_un <- i %in% unused
    add(cellbox(rna_x, ys_r[i],
                fill = if (i == reused) FILL_R else if (is_un) "#FFFFFF" else FILL_N,
                col  = if (i == reused) RED else CELL,
                lwd  = if (i == reused) LWD_E else LWD_N,
                lty  = if (is_un) DASH else "solid"))
  }
  for (j in seq_len(N_ATAC)) add(cellbox(atac_x, ys_a[j]))

  for (j in seq_len(N_ATAC)) {
    i <- assign_map[j]
    add(gcurve(atac_x - BOX_W / 2 - 0.006, ys_a[j],
               rna_x  + BOX_W / 2 + 0.006, ys_r[i],
               col = RED, lwd = LWD_E,
               curvature = if (ys_r[i] > ys_a[j]) CURV else -CURV, len = ARR))
  }
  add(collab(RNA_LAB, rna_x), collab(ATAC_LAB, atac_x))
  add(namelab(ARM_NAME[["link"]]))
  gTree(children = do.call(gList, unlist(.g, recursive = FALSE)))
}

## (b) aggregation (metacell)
panel_b <- function(k = 3L) {
  .g <- list(); add <- function(...) { .g[[length(.g) + 1L]] <<- gList(...); invisible(NULL) }

  ys_r <- Y0 - (seq_len(N_RNA)  - 1) * (ROW_H + GAP)
  ys_a <- Y0 - (seq_len(N_ATAC) - 1) * (ROW_H + GAP)
  rna_x <- 0.105; agg_r_x <- 0.375; agg_a_x <- 0.625; atac_x <- 0.895
  y_mid <- mean(range(c(ys_r, ys_a)))

  sel_r <- c(1L, 3L, 5L)[seq_len(k)]
  sel_a <- c(1L, 3L, 4L)[seq_len(k)]

  for (i in seq_len(N_RNA))
    add(cellbox(rna_x, ys_r[i], w = BOX_W_B,
                fill = if (i %in% sel_r) FILL_N else "#FFFFFF",
                lty  = if (i %in% sel_r) "solid" else DASH))
  for (j in seq_len(N_ATAC))
    add(cellbox(atac_x, ys_a[j], w = BOX_W_B,
                fill = if (j %in% sel_a) FILL_N else "#FFFFFF",
                lty  = if (j %in% sel_a) "solid" else DASH))

  ## the two aggregates, taller than a cell and tinted in the metacell colour
  add(cellbox(agg_r_x, y_mid, w = BOX_W_B, h = ROW_H * 1.45,
              fill = FILL_A, col = PURPLE, lwd = LWD_E),
      cellbox(agg_a_x, y_mid, w = BOX_W_B, h = ROW_H * 1.45,
              fill = FILL_A, col = PURPLE, lwd = LWD_E))

  for (i in sel_r)
    add(gcurve(rna_x + BOX_W_B / 2 + 0.006, ys_r[i],
               agg_r_x - BOX_W_B / 2 - 0.006, y_mid,
               col = PURPLE, lwd = LWD_E,
               curvature = if (ys_r[i] > y_mid) CURV else -CURV, len = ARR))
  for (j in sel_a)
    add(gcurve(atac_x - BOX_W_B / 2 - 0.006, ys_a[j],
               agg_a_x + BOX_W_B / 2 + 0.006, y_mid,
               col = PURPLE, lwd = LWD_E,
               curvature = if (ys_a[j] > y_mid) -CURV else CURV, len = ARR))

  ## the join between the two averages; arrowheads both ends, it carries no direction
  add(linesGrob(c(agg_r_x + BOX_W_B / 2 + 0.008, agg_a_x - BOX_W_B / 2 - 0.008),
                c(y_mid, y_mid),
                arrow = arrow(length = unit(ARR, "mm"), type = "closed", ends = "both"),
                gp = gpar(col = PURPLE, lwd = LWD_E, fill = PURPLE)))

  add(collab(RNA_LAB, rna_x), collab(ATAC_LAB, atac_x))
  add(namelab("aggregation  (metacell)"))
  gTree(children = do.call(gList, unlist(.g, recursive = FALSE)))
}

## (c) set composition and containment -- ALL SHARED REGULATORS, no module
LEVELS_C <- c(TF = "transcription factors", regulon = "regulons")
ARMS_C   <- ARM_NAME
SEG_SHARED <- "shared by both constructions"

cc <- do.call(rbind, lapply(STAGES, function(s)
  do.call(rbind, lapply(names(LEVELS_C), function(lv)
    do.call(rbind, lapply(names(ARMS_C), function(a) {
      tot <- mf_num(sprintf("amlcmp.%s.n_%s.%s", lv, if (a == "link") "link" else "pm", s))
      shr <- mf_num(sprintf("amlcmp.%s.n_shared.%s", lv, s))
      data.frame(
        stage = factor(s, levels = STAGES),
        level = factor(LEVELS_C[[lv]], levels = unname(LEVELS_C)),
        arm   = factor(ARMS_C[[a]],    levels = unname(ARMS_C)),
        segment = factor(c("shared", "unique"), levels = c("shared", "unique")),
        n = c(shr, tot - shr), total = tot, shared = shr,
        stringsAsFactors = FALSE)
    }))))))

## SOLID SEGMENTS (2026-09-13). The shared block is the set's non-arm grey in both bars; the
## unique block takes its arm's colour. Fill key: shared, then one "unique to <arm>" per arm.
UNIQ_LAB <- setNames(sprintf("unique to %s", unname(ARMS_C)), unname(ARMS_C))
cc$fillkey <- ifelse(cc$segment == "shared", SEG_SHARED, UNIQ_LAB[as.character(cc$arm)])
FILL_KEYS <- c(SEG_SHARED, unname(UNIQ_LAB))
cc$fillkey <- factor(cc$fillkey, levels = FILL_KEYS)
FIL_C <- setNames(c(GREY, RED, PURPLE), FILL_KEYS)

ends <- unique(cc[, c("stage", "level", "arm", "total", "shared")])
ends$lab <- ifelse(is.na(ends$total), "[VALUE OWED]",
                   sprintf("%d / %d", as.integer(ends$shared), as.integer(ends$total)))
ends$pct <- ifelse(is.na(ends$total), "",
                   sprintf("%.0f%%", 100 * ends$shared / ends$total))

## shared block at the ORIGIN in every bar, so the two identical blocks compare by eye
p_c <- ggplot(cc, aes(n, arm, fill = fillkey, colour = fillkey)) +
  geom_col(width = 0.62, linewidth = 0.25, position = position_stack(reverse = TRUE)) +
  geom_text(data = ends, inherit.aes = FALSE,
            aes(x = shared / 2, y = arm, label = pct), colour = "#FFFFFF", size = gs(1.75)) +
  geom_text(data = ends, inherit.aes = FALSE,
            aes(x = total, y = arm, label = lab), hjust = -0.14, size = gs(1.75), colour = INK) +
  scale_fill_manual(values = FIL_C, breaks = FILL_KEYS, name = NULL) +
  scale_colour_manual(values = FIL_C, breaks = FILL_KEYS, guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.24)),
                     breaks = scales::pretty_breaks(4)) +
  scale_y_discrete(limits = rev(unname(ARMS_C))) +
  facet_grid(stage ~ level, scales = "free_x", switch = "y") +
  labs(x = NULL, y = NULL) +
  theme_tight() +
  theme(strip.placement = "outside",
        strip.text.x = element_text(size = ts(5.8), colour = INK, margin = margin(0, 0, 1.5, 0)),
        strip.text.y.left = element_text(angle = 0, size = ts(5.8), colour = INK,
                                         margin = margin(0, 2.5, 0, 0)),
        panel.spacing.x = unit(16, "pt"), panel.spacing.y = unit(3, "pt"),
        axis.line.x = element_line(colour = INK, linewidth = 0.3),
        axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.key.height = unit(6, "pt"), legend.key.width = unit(7, "pt"),
        legend.margin = margin(3, 0, 0, 0), legend.box.spacing = unit(2, "pt"),
        plot.margin = margin(9, 4, 2, 4))

## (d) the dissociation, over the IDENTICAL regulon set
ACT_CSV <- "analysis/out/f1_celltype_agreement_allshared.csv"
TGT_CSV <- "analysis/out/f1_shared_regulon_jaccard.csv"

act <- read.csv(ACT_CSV, stringsAsFactors = FALSE)
act <- act[act$level == "regulon", c("stage", "unit", "pearson_r")]
names(act) <- c("stage", "regulon", "value")

tgt <- read.csv(TGT_CSV, stringsAsFactors = FALSE)
tgt <- tgt[tgt$arm_A == "link" & tgt$arm_B == "pseudomultiome",
           c("stage", "regulon", "jaccard_genes")]
names(tgt) <- c("stage", "regulon", "value")

## set identity, asserted on the SETS, not their sizes
for (s in STAGES) {
  ra <- sort(act$regulon[act$stage == s])
  rt <- sort(tgt$regulon[tgt$stage == s])
  if (!setequal(ra, rt))
    stop(sprintf(paste0("panel (d): REGULON-SET IDENTITY ASSERTION FAILED at %s. ",
                        "activity n = %d, target n = %d, symmetric difference = %d."),
                 s, length(ra), length(rt),
                 length(union(setdiff(ra, rt), setdiff(rt, ra)))))
  message(sprintf("  [panel d] regulon-set identity %s: PASS (%d regulons, sets equal)",
                  s, length(ra)))
}

for (s in STAGES) {
  chk <- list(
    list(k = sprintf("amlcmp.act.regulon.median.%s", s),
         v = median(act$value[act$stage == s]), what = "activity median"),
    list(k = sprintf("amlcmp.act.regulon.n.%s", s),
         v = sum(act$stage == s), what = "activity n"),
    list(k = sprintf("amlcmp.targetJ.regulon.median.%s", s),
         v = median(tgt$value[tgt$stage == s]), what = "target median"),
    list(k = sprintf("amlcmp.targetJ.regulon.n.%s", s),
         v = sum(tgt$stage == s), what = "target n"))
  for (x in chk) {
    m <- mf_num(x$k)
    if (!isTRUE(all.equal(as.numeric(x$v), m, tolerance = 1e-9)))
      stop(sprintf("panel (d): %s at %s is %.10g in the CSV but %.10g in the manifest (%s)",
                   x$what, s, x$v, m, x$k))
  }
  message(sprintf("  [panel d] manifest gate %s: PASS (medians and n agree)", s))
}

CUR <- c(act = "cell-type activity agreement  (Pearson r)",
         tgt = "target-set agreement  (Jaccard)")
dd <- rbind(cbind(act, currency = CUR[["act"]]), cbind(tgt, currency = CUR[["tgt"]]))
dd$currency <- factor(dd$currency, levels = unname(CUR))
NLAB <- sapply(STAGES, function(s) sum(act$stage == s))
STAGE_LAB <- sprintf("%s\nn = %d", STAGES, NLAB[STAGES])
dd$stage <- factor(sprintf("%s\nn = %d", dd$stage, NLAB[dd$stage]), levels = STAGE_LAB)

## ONE NEUTRAL COLOUR FOR BOTH CURRENCIES (2026-09-13). Both histograms are agreements BETWEEN
## the two arms, so neither may wear an arm's colour, and the currency is already named by the
## strip title and separated by facet position. A thin white edge keeps the bins legible.
HIST_FILL <- GREY
BW_ACT <- 0.05
BW_TGT <- 0.02

## fixed ranges: activity spans at least -0.5 to 1.0 (a floor, not a clip); target runs to 1.0
XLIM <- do.call(rbind, lapply(STAGE_LAB, function(sl) data.frame(
  stage = factor(sl, levels = STAGE_LAB),
  currency = factor(c(CUR[["act"]], CUR[["act"]], CUR[["tgt"]], CUR[["tgt"]]),
                    levels = unname(CUR)),
  value = c(-0.5, 1.0, 0.0, 1.0))))

p_d <- ggplot(dd, aes(value)) +
  geom_blank(data = XLIM, aes(x = value)) +
  geom_histogram(data = dd[dd$currency == CUR[["act"]], ], binwidth = BW_ACT,
                 boundary = 0, fill = HIST_FILL, colour = "#FFFFFF", linewidth = 0.15) +
  geom_histogram(data = dd[dd$currency == CUR[["tgt"]], ], binwidth = BW_TGT,
                 boundary = 0, fill = HIST_FILL, colour = "#FFFFFF", linewidth = 0.15) +
  scale_x_continuous(breaks = scales::pretty_breaks(6), labels = lab_nolead,
                     expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.06))) +
  facet_grid(stage ~ currency, scales = "free_x", switch = "y") +
  labs(x = NULL, y = "regulons") +
  theme_tight() +
  theme(strip.placement = "outside",
        strip.text.x = element_text(size = ts(5.8), colour = INK, margin = margin(0, 0, 2, 0)),
        strip.text.y.left = element_text(angle = 0, size = ts(5.8), colour = INK,
                                         lineheight = 0.95, margin = margin(0, 3, 0, 0)),
        panel.spacing.x = unit(20, "pt"), panel.spacing.y = unit(13, "pt"),
        axis.line.x = element_line(colour = INK, linewidth = 0.3),
        axis.line.y = element_line(colour = INK, linewidth = 0.3),
        legend.position = "none",
        plot.margin = margin(9, 6, 2, 4))

## compose


fig <- (wrap_elements(full = panel_a()) | wrap_elements(full = panel_b())) /
       p_c / p_d +
  plot_layout(heights = c(0.66, 0.94, 0.86))
fig <- tag_paper(fig)

writable <- function(f) {
  if (!file.exists(f)) return(TRUE)
  con <- suppressWarnings(try(file(f, open = "ab"), silent = TRUE))
  if (inherits(con, "try-error")) return(FALSE)
  close(con); TRUE
}
target <- OUT
if (!writable(paste0(OUT, ".pdf"))) {
  target <- paste0(OUT, ".NEW")
  message("  !! ", OUT, ".pdf is LOCKED by another process (a PDF viewer?). Writing ", target, ".pdf")
}
save_paper(fig, target, w = COL2, h = 172 * MM)
mf_write_keys(OUT)
message("=== Figure 1 done ===")
