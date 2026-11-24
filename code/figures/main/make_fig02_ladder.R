#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2); library(grid); library(patchwork); library(scales)
})
source("figures/shared/theme_paper.R")
source("figures/shared/manifest.R")
source("figures/shared/pairing_styles.R")   # tag_v2(), save_paper(), COL2, MM
source("figures/shared/number_labels.R")   # lab_nolead(): the set's no-leading-zero convention

OUT <- "figures/main/fig02_ladder"
message("\n=== Figure 2: the ladder on two true multiomes (four panels) ===")
INK <- "#111111"

ARMS <- data.frame(
  slot  = c("rnd_1to1", "inf_matched", "inf_free", "inf_reuse", "met_k5", "met_k10",
            "inf_r201", "rnd_r201"),
  fam   = c("rnd", "inf", "inf", "inf", "met", "met", "inf", "rnd"),
  core  = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE),
  label = c("random 1to1 in stratum", "inference 1to1 in stratum",
            "inference 1to1", "inference with reuse",
            "metacell k = 5", "metacell k = 10",
            "inference with 2.01x reuse  (BMMC)", "random with 2.01x reuse  (BMMC)"),
  pbmc  = c("rand_subcluster_1to1", "link_subcluster_1to1", "link_1to1", "link",
            "pseudomultiome", "pseudomultiome_k10", NA, NA),
  bmmc  = c("rand_celltype_x_donor_1to1", "link_x_celltype_donor_1to1", "link_1to1",
            "link", "pseudomultiome", "pseudomultiome_k10",
            "link_x_celltype_donor", "rand_celltype_donor_matched"),
  stringsAsFactors = FALSE)

## THE ORDER. Defined once. Asserted in every panel; ggplot's factor inference is
## never trusted to reproduce it, because with the labels gone order is the only
## within-group cue a reader has.
SLOT_ORDER <- ARMS$slot                       # legend order, and bar order
CORE_ORDER <- ARMS$slot[ARMS$core]            # the six drawn in (a), (b), (d)
SLOT_LAB   <- setNames(ARMS$label, ARMS$slot)
SLOT_FAM   <- setNames(ARMS$fam,   ARMS$slot)

## x POSITION: adjacent inside a family, a gap between families. Computed from the order
## itself so it cannot fall out of step with it.
fam_pos <- function(slots) {
  f <- SLOT_FAM[slots]
  x <- numeric(length(slots)); cur <- 0; prev <- NA_character_
  for (i in seq_along(slots)) {
    cur <- cur + if (is.na(prev) || f[i] == prev) 1 else 1.9
    x[i] <- cur; prev <- f[i]
  }
  setNames(x, slots)
}
XPOS_CORE <- fam_pos(CORE_ORDER)

BAR_W <- 0.80

FAM_HUE <- c(inf = "red", rnd = "blue", met = "purple")
PAL_F2 <- c(
  rnd_1to1    = "#4292C6",   # L* 54
  inf_matched = "#66000C",   # L* 19
  inf_free    = "#D11F1F",   # L* 45
  inf_reuse   = "#FDC5AE",   # L* 84
  met_k5      = "#5C3696",   # L* 32
  met_k10     = "#ACAAD1",   # L* 71
  inf_r201    = "#FB6B4B",   # L* 63
  rnd_r201    = "#99C7E0")   # L* 78
FAM_DARK <- c(inf = "#66000C", rnd = "#08519C", met = "#3F007D")
## OUTLINE <- setNames(unname(FAM_DARK[SLOT_FAM[SLOT_ORDER]]), SLOT_ORDER)   # v2
## LWD_F2  <- c(rnd_1to1 = 0.80, inf_matched = 0.80, inf_free = 0.52, inf_reuse = 0.26,
##              met_k5 = 0.80, met_k10 = 0.32, inf_r201 = 0.38, rnd_r201 = 0.38)  # v2
OUTLINE <- PAL_F2[SLOT_ORDER]                                  # outline colour IS the fill
LWD_F2  <- setNames(rep(0.35, length(SLOT_ORDER)), SLOT_ORDER)
stopifnot(identical(unname(OUTLINE), unname(PAL_F2[SLOT_ORDER])),
          identical(names(OUTLINE), SLOT_ORDER))
SHP_F2 <- c(inf = 24, rnd = 21, met = 23)     # (c) is a scatter: family also takes a shape

## ACCESSIBILITY GATE -- recomputed on every build, printed, never asserted from memory.
## Within-family separation is what governs: between families, hue AND block position
## AND marker shape are all available, so colour is not carrying that distinction alone.
palette_check <- function(P, fam) {
  ok <- requireNamespace("colorspace", quietly = TRUE) &&
        requireNamespace("farver", quietly = TRUE)
  if (!ok) { message("  [a11y] colorspace/farver absent -- CHECK NOT RUN"); return(invisible(NULL)) }
  lab <- function(h) farver::convert_colour(t(grDevices::col2rgb(h)), "rgb", "lab")
  dE  <- function(a, b) farver::compare_colour(lab(a), lab(b), "lab", method = "cie2000")[1, 1]
  Ls  <- lab(P)[, 1]
  message(sprintf("  [a11y] L* by slot: %s",
                  paste(sprintf("%s=%.0f", names(P), Ls), collapse = "  ")))
  worst <- Inf
  for (mode in c("normal", "deutan", "protan", "tritan", "greyscale")) {
    Q <- switch(mode, normal = P, deutan = colorspace::deutan(P),
                protan = colorspace::protan(P), tritan = colorspace::tritan(P),
                greyscale = colorspace::desaturate(P))
    wi <- Inf; bi <- Inf; wp <- ""; bp <- ""
    for (i in seq_len(length(Q) - 1)) for (j in (i + 1):length(Q)) {
      d <- dE(Q[i], Q[j])
      if (fam[names(P)[i]] == fam[names(P)[j]]) {
        if (d < wi) { wi <- d; wp <- paste(names(P)[i], names(P)[j], sep = "/") }
      } else if (d < bi) { bi <- d; bp <- paste(names(P)[i], names(P)[j], sep = "/") }
    }
    worst <- min(worst, wi)
    message(sprintf("  [a11y] %-10s dE2000  within-family %5.1f (%-24s)  between %5.1f (%s)",
                    mode, wi, wp, bi, bp))
  }
  if (worst < 10)
    stop(sprintf("palette gate FAILED: within-family separation falls to %.1f dE2000; ",
                 worst), "lightness is the only channel there and 10 is the floor")
  message(sprintf("  [a11y] GATE PASS -- worst within-family separation %.1f dE2000 ",
                  worst), "over all five vision modes")
}
palette_check(PAL_F2, SLOT_FAM)

ORC <- "orc"
PAL_F2[ORC]   <- unname(PAL[["context"]])
OUTLINE[ORC]  <- unname(PAL[["context"]])
LWD_F2[ORC]   <- 0.35
SLOT_FAM[ORC] <- "orc"
SLOT_LAB[ORC] <- "RNA nearest-neighbor comparator"
FILL_ORDER <- c(SLOT_ORDER, ORC)                     # every level the fill scale must resolve
LEG_ORDER  <- c(CORE_ORDER, ORC)                     # what the legend shows, in order
B_ORDER    <- c("rnd_1to1", "inf_matched", ORC)      # panel (c), floor -> observed -> ceiling

## ORDER GATE -- run by every panel that draws bars.
assert_order <- function(d, panel) {
  stopifnot(all(as.character(d$slot) %in% SLOT_ORDER))
  for (ds in unique(as.character(d$dataset))) {
    g <- d[as.character(d$dataset) == ds, ]
    got <- as.character(g$slot)[order(g$xpos)]
    got <- got[!duplicated(got)]
    want <- CORE_ORDER[CORE_ORDER %in% got]
    if (!identical(got, want))
      stop(sprintf("BAR ORDER VIOLATED in panel (%s), dataset %s\n  drawn: %s\n  want : %s",
                   panel, ds, paste(got, collapse = " "), paste(want, collapse = " ")))
  }
  if (!identical(levels(d$slot), SLOT_ORDER))
    stop(sprintf("panel (%s): factor levels are not SLOT_ORDER", panel))
  message(sprintf("  [order] panel (%s): PASS -- %s, identical in %d facet column(s)",
                  panel, paste(CORE_ORDER, collapse = " > "),
                  length(unique(d$dataset))))
  invisible(TRUE)
}

DATASETS <- c("PBMC", "BMMC")
NS      <- c(PBMC = "pbmcx", BMMC = "bmmc")
NS_AVG  <- c(PBMC = "pbmc",  BMMC = "bmmc")

arm_of <- function(ds, slot) ARMS[[tolower(ds)]][match(slot, ARMS$slot)]

k_rung <- function(ds, stem, slot) sprintf("%s.%s.%s", NS[[ds]], stem, arm_of(ds, slot))
k_l6   <- function(ds, slot)       sprintf("%s.L6.within.%s", NS[[ds]], arm_of(ds, slot))
k_avg  <- function(ds, slot)       sprintf("%s.celltype_r.%s", NS_AVG[[ds]], arm_of(ds, slot))
k_trip <- function(ds, slot)       sprintf("%s.n_triplets.%s", NS[[ds]], arm_of(ds, slot))
k_ref  <- function(ds)             sprintf("%s.truth.n_triplets", NS[[ds]])
## accuracy: BMMC has one family; PBMC's full-set column lives in two namespaces
k_acc  <- function(ds, slot) {
  a <- arm_of(ds, slot)
  if (ds == "BMMC") return(sprintf("bmmc.selfhit.%s", a))
  k <- sprintf("pbmc.exact.%s.full", a)
  if (mf_has(k)) k else sprintf("pbmcx.exact.%s", a)
}
## a value, or NA when the manifest says the quantity does not exist for that arm
getv <- function(key) if (mf_has(key)) mf_num(key) else NA_real_

theme_tight <- function() theme_paper() +
  theme(axis.title.x = element_text(margin = margin(t = 1.5)),
        axis.title.y = element_text(size = 5.6, margin = margin(r = 1.5)),
        legend.box.spacing = unit(1.5, "pt"), legend.margin = margin(0, 0, 0, 0),
        legend.background = element_blank(),
        strip.text = element_text(size = 6.0, margin = margin(1.5, 1.5, 1.5, 1.5)))

## bars carry no x-axis text at all -- that is the point of the shared legend
no_x_text <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                   axis.title.x = element_blank())

SC_FILL <- scale_fill_manual(values = PAL_F2, breaks = LEG_ORDER, limits = FILL_ORDER,
                             labels = SLOT_LAB[LEG_ORDER], name = NULL, drop = FALSE,
                             guide = guide_legend(
                               order = 1, ncol = 2, byrow = TRUE,
                               override.aes = list(colour = unname(OUTLINE[LEG_ORDER]),
                                                   linewidth = unname(LWD_F2[LEG_ORDER]))))
SC_COL  <- scale_colour_manual(values = OUTLINE, breaks = FILL_ORDER, limits = FILL_ORDER,
                               guide = "none")
SC_LWD  <- scale_linewidth_manual(values = LWD_F2, breaks = FILL_ORDER,
                                  limits = FILL_ORDER, guide = "none")
## marks that are not arms. Both levels are declared in both panels that use the scale so
## the collected guide carries them once, in this order.
## Short by necessity: at 68 characters this key pushed the legend past the figure width and
## clipped the two rightmost keys off the page. The arm it belongs to and the exact value are
## in the caption; what the KEY has to carry is that it is structural and on another estimator.
MK_REF  <- "structural identity (group-mean estimator)"
MK_PAIR <- "capacity-matched pair, 2.01x reuse"
## DOTTED for the structural reference, solid for the panel-(c) connector. The dot is
## deliberately not the dash used anywhere else: this line is a different estimator on a
## different unit and must not read as one more member of any series.
SC_LTY  <- scale_linetype_manual(values = setNames(c("12", "solid"), c(MK_REF, MK_PAIR)),
                                 breaks = c(MK_REF, MK_PAIR), limits = c(MK_REF, MK_PAIR),
                                 drop = FALSE, name = NULL,
                                 guide = guide_legend(order = 2, ncol = 1))

PAIR_ORDER <- c("inf_r201", "rnd_r201")
stopifnot(all(PAIR_ORDER %in% SLOT_ORDER), !any(PAIR_ORDER %in% CORE_ORDER))
SC_SHP  <- scale_shape_manual(
  values = setNames(unname(SHP_F2[SLOT_FAM[PAIR_ORDER]]), PAIR_ORDER),
  breaks = PAIR_ORDER, limits = PAIR_ORDER, labels = SLOT_LAB[PAIR_ORDER],
  name = NULL, drop = FALSE,
  guide = guide_legend(order = 3, ncol = 1,
                       override.aes = list(fill = unname(PAL_F2[PAIR_ORDER]),
                                           colour = unname(OUTLINE[PAIR_ORDER]),
                                           size = 1.5, stroke = 0.42)))

## (a) the ladder -- five rungs, two datasets, six arms
RUNGS <- data.frame(
  stem = c("L1.topn_jaccard", "L2.topn_jaccard", "L3.genes", "L4.jaccard", "L5.jaccard"),
  lab  = c("L1  region-gene", "L2  TF-gene", "L3  target sets", "L4  eRegulons", "L5  TFs"),
  stringsAsFactors = FALSE)

A <- do.call(rbind, lapply(DATASETS, function(ds)
  do.call(rbind, lapply(seq_len(nrow(RUNGS)), function(j)
    do.call(rbind, lapply(CORE_ORDER, function(s) data.frame(
      dataset = ds, rung = RUNGS$lab[j], slot = s, xpos = XPOS_CORE[[s]],
      v = getv(k_rung(ds, RUNGS$stem[j], s)), stringsAsFactors = FALSE)))))))
A$slot    <- factor(A$slot, levels = SLOT_ORDER)
A$dataset <- factor(A$dataset, levels = DATASETS)
A$rung    <- factor(A$rung, levels = RUNGS$lab)
assert_order(A, "a")

p_a <- ggplot(A, aes(xpos, v, fill = slot, colour = slot, linewidth = slot)) +
  geom_col(width = BAR_W) +
  geom_text(aes(label = ifelse(is.na(v), "", sub("^0", "", sprintf("%.2f", v)))),
            vjust = -0.42, size = 1.6, colour = INK, show.legend = FALSE) +
  facet_grid(dataset ~ rung, switch = "y", axes = "all_x", axis.labels = "margins") +
  SC_FILL + SC_COL + SC_LWD +
  scale_x_continuous(expand = expansion(add = 0.7)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                     labels = c("0", ".5", "1"),
                     expand = expansion(mult = c(0, 0.16))) +
  labs(x = NULL, y = "agreement with the native reference  (Jaccard)") +
  theme_tight() + no_x_text +
  ## BOTTOM MARGIN 5 pt, NOT 2. A y-axis tick label is centred on its break, so half the "0"
  ## glyph sits BELOW y = 0 -- and y = 0 is the panel's bottom edge here (lower expansion 0).
  ## At 2 pt the plot boundary cut the bottom off the BMMC "0": measured 28 raster rows against
  ## 30 for every other tick label in the panel at 600 dpi. Invisible in the 300 dpi PNG and
  ## visible in the PDF, which is the file that ships. 5 pt clears a 5.0 pt label with margin.
  theme(panel.spacing.x = unit(3.5, "pt"), panel.spacing.y = unit(4, "pt"),
        strip.placement = "outside", axis.text.y = element_text(size = 5.0),
        plot.margin = margin(10, 3, 5, 2))

D <- do.call(rbind, lapply(DATASETS, function(ds) {
  ref <- mf_num(k_ref(ds))
  do.call(rbind, lapply(CORE_ORDER, function(s) data.frame(
    dataset = ds, slot = s, xpos = XPOS_CORE[[s]],
    v = getv(k_trip(ds, s)) / ref, stringsAsFactors = FALSE)))
}))
D$slot    <- factor(D$slot, levels = SLOT_ORDER)
D$dataset <- factor(D$dataset, levels = DATASETS)
assert_order(D, "b")

p_b <- ggplot(D, aes(xpos, v, fill = slot, colour = slot, linewidth = slot)) +
  geom_hline(yintercept = 1, linetype = "22", colour = INK, linewidth = 0.3) +
  geom_col(width = BAR_W) +
  geom_text(aes(label = sub("^0", "", sprintf("%.2f", v))),
            vjust = -0.42, size = 1.6, colour = INK, show.legend = FALSE) +
  facet_wrap(~ dataset, nrow = 1) +
  SC_FILL + SC_COL + SC_LWD +
  scale_x_continuous(expand = expansion(add = 0.7)) +
  scale_y_continuous(breaks = c(0, 0.5, 1, 1.5), labels = lab_nolead,
                     expand = expansion(mult = c(0, 0.10))) +
  labs(x = NULL, y = "eRegulon triplets built,\nrelative to native reference") +
  theme_tight() + no_x_text +
  theme(panel.spacing.x = unit(5, "pt"), axis.text.y = element_text(size = 5.0),
        plot.margin = margin(10, 3, 2, 2))

B_SLOTS <- c("rnd_1to1", "inf_matched", "orc")
## x by hand, not by fam_pos(): the two ARMS sit adjacent because they are the matched pair the
## panel contrasts, and the ORACLE is set apart because it is not an arm.
B_XPOS  <- c(rnd_1to1 = 1, inf_matched = 2, orc = 3.4)

B <- do.call(rbind, lapply(DATASETS, function(ds) data.frame(
  dataset = ds, slot = B_SLOTS, xpos = unname(B_XPOS[B_SLOTS]),
  v = c(getv(k_l6(ds, "rnd_1to1")), getv(k_l6(ds, "inf_matched")),
        getv(sprintf("%s.L6.ceiling.rna_nn_within_stratum", NS[[ds]]))),
  stringsAsFactors = FALSE)))
B$slot    <- factor(B$slot, levels = B_ORDER)
B$dataset <- factor(B$dataset, levels = DATASETS)
if (any(is.na(B$v)))
  stop("panel (c): a bar value is missing; run figures/data/build_manifest.py")

for (ds in DATASETS) for (w in c("floor", "observed")) {
  .k <- sprintf("%s.L6.ceiling.%s_reproduces_manifest", NS[[ds]], w)
  if (!identical(mf(.k, quiet = TRUE), TRUE))
    stop(sprintf("panel (c): %s is not TRUE -- the ceiling was not computed with the same L6 ",
                 .k), "estimator as the bars, so the three are not on one axis")
}
message(sprintf(paste("  [panel c] one-estimator gate: PASS -- floor and observed reproduce the",
                      "manifest under the ceiling script's own L6 on both datasets")))
for (ds in DATASETS)
  message(sprintf("  [panel c] %s floor %+.4f | observed %+.4f | ceiling %+.4f (oracle reuse %.4fx over %s cells)",
                  ds, B$v[B$dataset == ds & B$slot == "rnd_1to1"],
                  B$v[B$dataset == ds & B$slot == "inf_matched"],
                  B$v[B$dataset == ds & B$slot == "orc"],
                  mf_num(sprintf("%s.L6.ceiling.oracle_reuse", NS[[ds]])),
                  format(mf_num(sprintf("%s.L6.ceiling.n_cells_evaluated", NS[[ds]])),
                         big.mark = ",")))

if (!identical(mf("bmmc.celltype_degenerate.rand_celltype_x_donor_1to1", quiet = TRUE), TRUE))
  stop("panel (c): BMMC's degeneracy flag is not TRUE; the reference line's status must be re-derived")
if (!isTRUE(all.equal(mf_num("pbmc.reuse.rand_subcluster_1to1.reported"), 1)))
  stop("panel (c): PBMC's rand_subcluster_1to1 is not a bijection, so its 1.0000 is not an identity")
RF <- do.call(rbind, lapply(DATASETS, function(ds) data.frame(
  dataset = ds, y = mf_num(k_avg(ds, "rnd_1to1")), key = MK_REF, stringsAsFactors = FALSE)))
RF$dataset <- factor(RF$dataset, levels = DATASETS)
if (!all(abs(RF$y - 1) < 1e-9))
  stop(sprintf("panel (c): a reference value is %.12f, not 1 -- it is being drawn as an identity",
               RF$y[which.max(abs(RF$y - 1))]))
message(sprintf(paste("  [panel c] identity gate: PASS -- group-mean agreement of the",
                      "random-matched arm is exactly 1.000000 on both datasets",
                      "(BMMC flagged in-file; PBMC asserted from reuse = 1.0 within subcluster)")))

p_c <- ggplot(B, aes(xpos, v, fill = slot, colour = slot, linewidth = slot)) +
  geom_hline(yintercept = 0, colour = "#BFBFBF", linewidth = 0.3) +
  ## DOTTED, and never a bar. A different estimator on a different unit must not be able to be
  ## read as one more member of the bars' series.
  geom_hline(data = RF, aes(yintercept = y, linetype = key), colour = INK,
             linewidth = 0.4, inherit.aes = FALSE, key_glyph = "path") +
  geom_col(width = BAR_W) +
  geom_text(aes(label = sub("^0", "", sprintf("%.2f", v))),
            vjust = -0.42, size = 1.6, colour = INK, show.legend = FALSE) +
  facet_wrap(~ dataset, nrow = 1) +
  SC_FILL + SC_COL + SC_LWD + SC_LTY +
  scale_x_continuous(expand = expansion(add = 0.8)) +
  ## 0 TO 1, the full range of the quantity. The bars sit low and the reference sits at the top;
  ## that distance is the panel and a truncated axis would erase it.
  ## SHORT TICK LABELS, as panel (a) uses. "0.25"/"0.75" are wider, and in a three-plot row
  ## patchwork charges that width to the row: the first render of this panel narrowed (d)
  ## enough to close its tightest bar seam from 3 px to 2 px, at the QC gate's failure edge.
  scale_y_continuous(limits = c(0, 1.04), breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", ".25", ".5", ".75", "1"), expand = expansion(0)) +
  labs(x = NULL, y = "within-cell-type activity agreement\n(Spearman \u03c1)") +
  theme_tight() + no_x_text +
  theme(panel.spacing.x = unit(5, "pt"), axis.text.y = element_text(size = 5.0),
        plot.margin = margin(10, 3, 5, 2))

C_PAIR <- SLOT_ORDER[SLOT_ORDER %in% c("inf_r201", "rnd_r201")]
stopifnot(identical(C_PAIR, c("inf_r201", "rnd_r201")))   # order, from SLOT_ORDER, once

## --- the match itself, asserted from the manifest -------------------------------------
## Equality is required to full double precision for the three quantities that define the
## match; the two cross-stratum rates must be exactly zero on both arms.
pair_arm <- function(s) arm_of("BMMC", s)
pv <- function(stem, s) mf_num(sprintf("bmmc.%s.%s", stem, pair_arm(s)))
for (stem in c("reuse", "distinct_rna", "pct_pool")) {
  a <- pv(stem, "inf_r201"); b <- pv(stem, "rnd_r201")
  if (!identical(a, b))
    stop(sprintf("panel (d): the pair is NOT matched on %s -- %.17g vs %.17g", stem, a, b))
  message(sprintf("  [panel d] matched on %-13s %.17g (both arms)", stem, a))
}
for (stem in c("cross_celltype", "cross_donor")) for (s in C_PAIR) {
  v <- pv(stem, s)
  if (v != 0) stop(sprintf("panel (d): %s is %g for %s, expected 0", stem, v, s))
}
message("  [panel d] cross_celltype = cross_donor = 0 on both arms")

## --- the two quantities ---------------------------------------------------------------
Q_ACC <- "exact-cell accuracy\n(cells per 1,000)"
Q_REC <- "L5 TF-set\nJaccard vs\nnative reference"
C_N   <- mf_num("bmmc.n_cells")

CQ <- do.call(rbind, lapply(C_PAIR, function(s) data.frame(
  slot  = s,
  fam   = SLOT_FAM[[s]],
  xpos  = match(s, C_PAIR),
  quant = c(Q_ACC, Q_REC),
  v     = c(getv(k_acc("BMMC", s)) * 1000, getv(k_rung("BMMC", "L5.jaccard", s))),
  stringsAsFactors = FALSE)))
CQ$slot  <- factor(CQ$slot,  levels = SLOT_ORDER)
CQ$quant <- factor(CQ$quant, levels = c(Q_ACC, Q_REC))
CQ$key   <- MK_PAIR
if (any(is.na(CQ$v))) stop("panel (d): a pair value is missing from the manifest")

## the order gate, run here too: same slot order in both facets, from SLOT_ORDER
for (q in levels(CQ$quant)) {
  g <- CQ[CQ$quant == q, ]
  got <- as.character(g$slot)[order(g$xpos)]
  if (!identical(got, C_PAIR))
    stop(sprintf("panel (d): slot order in facet %s is %s", q, paste(got, collapse = " > ")))
}
message(sprintf("  [order] panel (d): PASS -- %s, identical in %d facet(s)",
                paste(C_PAIR, collapse = " > "), length(levels(CQ$quant))))

.acc <- CQ$v[CQ$quant == Q_ACC]; .rec <- CQ$v[CQ$quant == Q_REC]
message(sprintf(paste("  [panel d] reuse 2.01x pair over n = %s cells: accuracy %.2f vs %.2f",
                      "per 1,000 (%.2fx) with L5 %.4f vs %.4f (%+.4f)"),
                format(C_N, big.mark = ","), .acc[1], .acc[2], .acc[1] / .acc[2],
                .rec[1], .rec[2], .rec[2] - .rec[1]))

## BOTH facets are anchored at 0, and Jaccard is additionally forced up to 1 so it is read
## on the same scale as panels (a) and (b). Neither slope is steepened by a truncated axis.
## Forcing the 0 on the Jaccard facet is not cosmetic: without it the facet minimum was the
## DATA minimum (0.5407) and, with zero bottom expansion, the lower point was drawn clipped
## by the panel border. Caught by rendering and looking, not by the build.
CBL <- data.frame(quant = factor(c(Q_ACC, Q_REC, Q_REC), levels = c(Q_ACC, Q_REC)),
                  xpos = 1, v = c(0, 0, 1))

p_d <- ggplot(CQ, aes(xpos, v)) +
  geom_line(aes(linetype = key), colour = INK, linewidth = 0.35, key_glyph = "path") +
  ## shape = slot, NOT fam. The shape VALUES are still SHP_F2[fam] -- family still picks the
  ## glyph -- but mapping the slot is what lets SC_SHP build a key per arm, so the legend
  ## shows this panel's orange triangle and light-blue circle as a triangle and a circle.
  geom_point(aes(fill = slot, colour = slot, shape = slot), size = 1.5, stroke = 0.42) +
  geom_text(aes(label = ifelse(quant == Q_ACC, sprintf("%.1f", v),
                               sub("^0", "", sprintf("%.2f", v)))),
            vjust = -1.05, size = 1.6, colour = INK, show.legend = FALSE) +
  geom_blank(data = CBL, aes(xpos, v)) +
  ## strip.position "left" + strip.placement "outside" puts the facet strip exactly where
  ## an axis title goes, so the strip IS this facet's y-axis label -- the same device
  ## panel (a) already uses for its dataset names. There is no panel title here.
  facet_wrap(~ quant, ncol = 1, scales = "free_y", strip.position = "left") +
  SC_FILL + SC_COL + SC_LTY + SC_SHP +
  scale_x_continuous(expand = expansion(add = 0.55)) +
  scale_y_continuous(breaks = scales::pretty_breaks(3), labels = lab_nolead,
                     expand = expansion(mult = c(0, 0.20))) +
  labs(x = NULL, y = NULL) +
  theme_tight() + no_x_text +
  theme(panel.spacing.y = unit(6, "pt"), axis.text.y = element_text(size = 5.0),
        strip.placement = "outside", strip.text.y.left = element_text(size = 5.6,
                                                                      angle = 90),
        plot.margin = margin(10, 3, 2, 2))



LEG <- data.frame(slot = factor(LEG_ORDER, levels = FILL_ORDER),
                  x = seq_along(LEG_ORDER), v = 1)
p_leg <- ggplot(LEG, aes(x, v, fill = slot, colour = slot, linewidth = slot)) +
  geom_col(width = BAR_W) +
  geom_hline(data = data.frame(y = 1, key = MK_REF), aes(yintercept = y, linetype = key),
             colour = INK, linewidth = 0.4, inherit.aes = FALSE, key_glyph = "path") +
  ## MK_PAIR is declared in SC_LTY's limits but is drawn only in panel (d). Without a
  ## layer for it here its key rendered BLANK -- a legend entry with a label and no glyph.
  ## This segment exists solely to give that key something to draw.
  geom_segment(data = data.frame(x = 1, xend = 2, y = 1, yend = 1, key = MK_PAIR),
               aes(x = x, xend = xend, y = y, yend = yend, linetype = key),
               colour = INK, linewidth = 0.35, inherit.aes = FALSE, key_glyph = "path") +
  ## The two panel-(c) arms. This layer maps ONLY shape: if it also mapped fill or colour it
  ## would match the fill guide too and stamp a point into all six rectangle keys. Their fill
  ## and outline reach the key through SC_SHP's per-key override.aes instead.
  geom_point(data = data.frame(x = 1, v = 1,
                               slot = factor(PAIR_ORDER, levels = SLOT_ORDER)),
             aes(x, v, shape = slot), colour = INK, inherit.aes = FALSE) +
  SC_FILL + SC_COL + SC_LWD + SC_LTY + SC_SHP +
  theme_void() +
  ## "top", not "bottom": the key sits in an empty panel, and at the bottom of that panel
  ## it left a band of white between the figure and its own legend.
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.justification = "left", legend.box = "horizontal",
        legend.box.just = "top", legend.key.size = unit(5.5, "pt"),
        legend.text = element_text(size = 5.0), legend.spacing.x = unit(5, "pt"),
        legend.margin = margin(0, 0, 0, 0),
        ## 6 pt of TOP margin, set on this plot alone. The 2026-08-26 build had 12 pt
        ## applied to every subplot with `&`, which put a white band across the figure;
        ## this is the legend's own margin and reaches nothing else. It is needed because
        ## (d) no longer carries an x-axis title, so the b|c|d row grew and the key's
        ## first row came to sit directly on panel (c)'s baseline.
        plot.margin = margin(6, 0, 0, 0))
## The key-supplying panel must not draw its own artwork. Done by moving the VIEWPORT off
## the data rather than by setting alpha on the layers: a legend key inherits the layer's
## aes_params, so alpha = 0 would erase the key along with the bar.
p_leg <- p_leg + coord_cartesian(xlim = c(100, 101), ylim = c(10, 11), expand = FALSE)

## The four data panels draw no legend of their own.
for (.nm in c("p_a", "p_b", "p_c", "p_d"))
  assign(.nm, get(.nm) + theme(legend.position = "none"))

FIG <- (p_a / (p_b | p_c | p_d) / p_leg) +
  plot_layout(heights = c(1.42, 1.06, 0.22)) +
  plot_annotation(tag_levels = list(c("a", "b", "c", "d", "")))
FIG <- FIG & theme(plot.tag = element_text(size = 9, face = "bold", family = "sans",
                                           colour = PAL[["truth"]], hjust = 0, vjust = 1),
                   plot.tag.position = c(0, 1.02))

writable <- function(f) {
  if (!file.exists(f)) return(TRUE)
  con <- suppressWarnings(try(file(f, open = "ab"), silent = TRUE))
  if (inherits(con, "try-error")) return(FALSE)
  close(con); TRUE
}
target <- OUT
if (!writable(paste0(OUT, ".pdf"))) {
  target <- paste0(OUT, ".NEW")
  message("  !! ", OUT, ".pdf is LOCKED. Writing ", target, ".pdf")
}
save_paper(FIG, target, w = COL2, h = 132 * MM)
mf_write_keys(OUT)
message("=== Figure 2 done ===\n")
