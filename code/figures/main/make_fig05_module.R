#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2); library(grid); library(patchwork); library(scales); library(jsonlite)
})
source("figures/shared/theme_paper.R")
source("figures/shared/manifest.R")
source("figures/shared/number_labels.R")   # lab_nolead(), fmt_nolead(): the set's no-leading-zero convention

## cairo_pdf cannot overwrite a PDF that another process holds open -- it dies with
## "error while writing to output stream" and LEAVES THE STALE FILE IN PLACE, which is worse
## than failing loudly (figures/BUILD_LOG.md records the same trap). AMLGRN_FIG5_OUT lets a build
## be written and verified beside the locked one rather than silently not happening.
OUT <- Sys.getenv("AMLGRN_FIG5_OUT", "figures/main/fig05_module")
message("\n=== Figure 5: the corrected arm's module of record (the 55) ===")
INK <- "#111111"
CORR <- "#4292C6"          # corrected arm, identical to Figures 2-4
VOCABS <- c("Diagnosis", "Relapse")
## "Diagnosis" / "Relapse" NAME A VOCABULARY in (a), (b), (d), (e) and (f), and a TIMEPOINT only in
## (c). Before 2026-09-13 both read as bare "Diagnosis" / "Relapse", and (c) drew Relapse in the
## navy that (e) uses for translation elongation. Every vocabulary label now says "vocabulary".
VOCLAB <- setNames(paste(VOCABS, "vocabulary"), VOCABS)

theme_tight <- function() theme_paper() +
  theme(axis.title.x = element_text(size = 6.0, margin = margin(t = 1.5)),
        axis.title.y = element_text(size = 6.0, margin = margin(r = 1.5)),
        axis.text = element_text(size = 5.2),
        strip.text = element_text(size = 5.8, margin = margin(1.5, 1.5, 1.5, 1.5)),
        legend.box.spacing = unit(1.5, "pt"), legend.margin = margin(0, 0, 0, 0),
        legend.background = element_blank(), legend.key.size = unit(5.5, "pt"),
        legend.text = element_text(size = 5.2),
        plot.margin = margin(10, 4, 5, 3))

## THE MODULE, DEFINED ONCE AND ASSERTED
MODJ <- jsonlite::fromJSON("modules/out/corrected_module.json", simplifyVector = FALSE)
M55  <- unlist(MODJ$regulons)
if (length(M55) != 55L) stop(sprintf("module of record is %d regulons, expected 55", length(M55)))
if (length(unlist(MODJ$TFs)) != 31L) stop("module of record is not 31 TFs")
message(sprintf("  [module] %d regulons / %d TFs from modules/out/corrected_module.json",
                length(M55), length(unlist(MODJ$TFs))))

N_EXT  <- mf_num("set55.gate5.n_paired_patients")     # 59, external replication cohort
N_DISC <- mf_num("s3r.cohort.n_paired_patients")      # 19, paired discovery cohort
message(sprintf("  [cohort] (a)-(c): external replication n = %d | (d)-(e): discovery n = %d",
                N_EXT, N_DISC))

## per-regulon attributes -- carries ax4_would_remove, the retired gate's counterfactual
S55 <- read.csv("analysis/out/f6_set55_per_regulon.csv", stringsAsFactors = FALSE,
                check.names = FALSE)
if (!setequal(S55$regulon, M55)) stop("f6_set55_per_regulon.csv does not match the module")

FORBIDDEN <- c("ax4", "tier1", "tier 1", "T4", "34-regulon", "retired")
assert_no_retired_module <- function(strs, where) {
  hit <- unlist(lapply(FORBIDDEN, function(k)
    strs[grepl(k, strs, fixed = TRUE)]))
  if (length(hit))
    stop(sprintf("Figure 5 standing rule: %s names the retired module or its gate -- %s",
                 where, paste(unique(hit), collapse = " ; ")))
  invisible(TRUE)
}

A <- do.call(rbind, lapply(VOCABS, function(v) {
  d <- read.csv(sprintf("adult/corrected/out/per_regulon_external_T3_repair_%s.csv", v),
                stringsAsFactors = FALSE)
  d$vocabulary <- v
  d
}))
if (!setequal(unique(A$regulon), M55)) stop("panel (a): the external file is not the 55")

## ORIENTATION GATE. If this ever fails, the axis label is wrong and the panel must not be drawn.
.agree <- toupper(as.character(A$agrees_with_internal_sign)) %in% c("TRUE", "T")
if (!identical(.agree, A$oriented_mean_delta > 0))
  stop("panel (a): agrees_with_internal_sign is NOT equivalent to oriented_mean_delta > 0; ",
       "the deltas may be raw rather than oriented and the axis label would be wrong")
message(sprintf("  [panel a] ORIENTATION GATE: PASS -- deltas are oriented; agreement == (delta > 0) for all %d rows",
                nrow(A)))
for (v in VOCABS)
  message(sprintf("  [panel a] %-9s %d of 55 agree with the discovery direction (%.1f%%)",
                  v, sum(.agree[A$vocabulary == v]), 100 * mean(.agree[A$vocabulary == v])))

ORD <- A$regulon[A$vocabulary == "Diagnosis"][order(A$oriented_mean_delta[A$vocabulary == "Diagnosis"])]
A$regulon <- factor(A$regulon, levels = ORD)
A$vocabulary <- factor(A$vocabulary, levels = VOCABS)

## the connecting segment makes the two vocabularies one object per regulon
ASEG <- data.frame(
  regulon = factor(ORD, levels = ORD),
  lo = pmin(A$oriented_mean_delta[match(ORD, as.character(A$regulon[A$vocabulary == "Diagnosis"]))],
            A$oriented_mean_delta[A$vocabulary == "Relapse"][
              match(ORD, as.character(A$regulon[A$vocabulary == "Relapse"]))]),
  hi = pmax(A$oriented_mean_delta[match(ORD, as.character(A$regulon[A$vocabulary == "Diagnosis"]))],
            A$oriented_mean_delta[A$vocabulary == "Relapse"][
              match(ORD, as.character(A$regulon[A$vocabulary == "Relapse"]))]))

LAB3 <- c("ERG_direct_-/-", "TFDP2_direct_-/-", "GABPB1_direct_-/+")
ALAB <- A[A$vocabulary == "Diagnosis" & as.character(A$regulon) %in% LAB3, ]
ALAB$txt <- sprintf("%s  %s", sub("_direct", "", as.character(ALAB$regulon)),
                    fmt_nolead(ALAB$oriented_mean_delta, 4, signed = TRUE))
ALAB$vj <- c(`ERG_direct_-/-` = 0.35, `TFDP2_direct_-/-` = 1.55,
             `GABPB1_direct_-/+` = 2.75)[as.character(ALAB$regulon)]
message(sprintf("  [panel a] labelled: %s", paste(sprintf(
  "%s %+.4f (p %.2g)", ALAB$regulon, ALAB$oriented_mean_delta, ALAB$p), collapse = " | ")))

p_a <- ggplot(A, aes(oriented_mean_delta, regulon)) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.3) +
  geom_segment(data = ASEG, aes(x = lo, xend = hi, y = regulon, yend = regulon),
               colour = "#BFBFBF", linewidth = 0.25, inherit.aes = FALSE) +
  geom_point(aes(shape = vocabulary), fill = CORR, colour = CORR, size = 0.85, stroke = 0.3) +
  geom_text(data = ALAB, aes(x = oriented_mean_delta, y = regulon, label = txt, vjust = vj),
            hjust = -0.12, size = 1.5, colour = INK, inherit.aes = FALSE) +
  ## filled = the gate keeps it; hollow = the gate would remove it. Identical convention to
  ## Figure 4(c), where the imported ax4 gate is drawn unfilled.
  scale_shape_manual(values = c(Diagnosis = 21, Relapse = 24), labels = VOCLAB, name = NULL,
                     guide = guide_legend(override.aes = list(fill = CORR, colour = CORR,
                                                              size = 1.4, stroke = 0.4))) +
  scale_x_continuous(limits = c(-0.075, 0.235), breaks = seq(-0.05, 0.20, 0.05),
                     labels = c("-.05", "0", ".05", ".10", ".15", ".20"),
                     expand = expansion(0)) +
  ## HEADROOM ABOVE THE TOP ROW (2026-09-13). At ggplot's default 0.6-row expansion the ERG label's
  ## ascenders were cut by the panel edge and the top regulon's marks ran into the legend above.
  scale_y_discrete(expand = expansion(add = c(0.6, 2.2))) +
  labs(y = NULL,
       x = sprintf("agreement with the direction predicted by the discovery cohort\n(oriented change in regulon activity; %d paired external patients)",
                   N_EXT)) +
  theme_tight() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = "top", legend.direction = "vertical",
        ## bottom margin 2, not -3: the longer "... vocabulary" keys otherwise sat on the panel's
        ## top edge (2026-09-13)
        legend.justification = "left", legend.box.margin = margin(0, 0, 2, -4))

message("=== panel (a) built ===")

## (b) per-patient module change -- 59 paired external patients, drawn per patient
## PER PATIENT, NOT A COHORT-MEAN BAR. The claim is that most patients move the same way, and a
## bar cannot show that; the previous build rejected a bar for the same reason.
B <- read.csv("adult/corrected/out/per_patient_deltas_T3_repair.csv", stringsAsFactors = FALSE)
B$vocabulary <- factor(B$vocabulary, levels = VOCABS)
if (nrow(B) != 2 * N_EXT)
  stop(sprintf("panel (b): expected %d rows (2 vocabularies x %d patients), found %d",
               2 * N_EXT, N_EXT, nrow(B)))

G5 <- jsonlite::fromJSON("adult/corrected/out/gate5_result_T3_repair.json", simplifyVector = FALSE)
BSUM <- do.call(rbind, lapply(VOCABS, function(v) data.frame(
  vocabulary = factor(v, levels = VOCABS),
  mean_delta = as.numeric(G5$by_vocabulary[[v]]$mean_delta),
  frac_pos   = as.numeric(G5$by_vocabulary[[v]]$frac_patients_positive),
  p_wil      = as.numeric(G5$by_vocabulary[[v]]$p_wilcoxon), stringsAsFactors = FALSE)))

for (v in VOCABS) {
  .o <- mean(B$delta[B$vocabulary == v]); .m <- BSUM$mean_delta[BSUM$vocabulary == v]
  if (abs(.o - .m) > 1e-9)
    stop(sprintf("panel (b): %s per-patient mean %.10f != committed %.10f", v, .o, .m))
}
message("  [panel b] per-patient file reproduces the committed means on both vocabularies")

## ALL THREE NULLS ARE DRAWN, not annotated in prose. Each is its own object at 2,000 draws:
## A permutes the Diagnosis/Relapse labels, B draws size-matched random gene sets, C draws random
## 55-regulon modules from the corrected arm's 201-regulon universe. C sits highest because a
## random module of this arm also moves; the observation clears it all the same.
NULLKEY <- c(nullA_label_permutation      = "null A: label permutation",
             nullB_random_gene_sets       = "null B: random gene sets",
             nullC_random_regulon_modules = "null C: random 55-regulon modules")
BNULL <- do.call(rbind, lapply(VOCABS, function(v)
  do.call(rbind, lapply(names(NULLKEY), function(k) {
    z <- G5$by_vocabulary[[v]][[k]]
    data.frame(vocabulary = factor(v, levels = VOCABS), null = unname(NULLKEY[k]),
               m = as.numeric(z$mean), s = as.numeric(z$sd),
               p = as.numeric(z$p_empirical), stringsAsFactors = FALSE)
  }))))
BNULL$null <- factor(BNULL$null, levels = unname(NULLKEY))
BNULL$x0 <- as.numeric(BNULL$vocabulary) - 0.36
BNULL$x1 <- as.numeric(BNULL$vocabulary) + 0.36
message(sprintf("  [panel b] %s", paste(sprintf(
  "%s mean %+.4f, %.1f%% positive, Wilcoxon %.3g", BSUM$vocabulary, BSUM$mean_delta,
  100 * BSUM$frac_pos, BSUM$p_wil), collapse = " | ")))
message(sprintf("  [panel b] nulls (mean +- SD, empirical p): %s", paste(sprintf(
  "%s/%s %.4f+-%.4f p %.4g", BNULL$vocabulary, sub(":.*", "", BNULL$null), BNULL$m, BNULL$s,
  BNULL$p), collapse = " | ")))
if (!all(abs(BNULL$p - BNULL$p[1]) < 1e-12))
  message("  [panel b] NOTE: the three nulls no longer share one empirical p; the caption must say so")

## THE SCALE MUST CONTAIN EVERY DATUM. A hard-coded limit silently dropped one patient
## (AML130, Diagnosis, +0.16362) on the first render -- ggplot reports that only as a warning,
## and a clipped point is data loss, not a display choice. The limit is computed from the data
## and the containment is asserted.
BLIM <- range(c(B$delta, BNULL$m, 0))
BLIM <- c(floor(BLIM[1] * 40) / 40, ceiling(BLIM[2] * 40) / 40)
if (any(B$delta < BLIM[1] | B$delta > BLIM[2]))
  stop("panel (b): the y scale would clip a patient")
message(sprintf("  [panel b] y scale [%.3f, %.3f] contains all %d patient values and all null means",
                BLIM[1], BLIM[2], nrow(B)))

## NULLS A AND B SHARE ONE LABEL, AND THAT IS ASSERTED RATHER THAN ASSUMED. Their dotted lines
## coincide at this scale; one label is only honest if the two lines are genuinely
## indistinguishable -- less than 1% of the panel's y range apart, in both columns. If a rerun ever
## moves them apart, the build stops and they must be labelled separately.
.ab <- vapply(VOCABS, function(v) {
  z <- BNULL[BNULL$vocabulary == v, ]
  abs(z$m[as.integer(z$null) == 1] - z$m[as.integer(z$null) == 2])
}, numeric(1))
if (any(.ab >= 0.01 * diff(BLIM)))
  stop(sprintf("panel (b): nulls A and B are %s apart, not indistinguishable at this scale; label them separately",
               paste(sprintf("%.4f", .ab), collapse = " / ")))
message(sprintf("  [panel b] nulls A and B coincide at this scale (|mean A - mean B| %s against 1%% of the y range, %.4f); one label",
                paste(sprintf("%.5f", .ab), collapse = " / "), 0.01 * diff(BLIM)))
## the labels sit right of the Relapse column and name both columns' lines, which are drawn in the
## same order in each: C above, A and B together below. That order is asserted.
BNLAB <- data.frame(
  y = c(mean(BNULL$m[BNULL$vocabulary == "Relapse" & as.integer(BNULL$null) %in% 1:2]),
        BNULL$m[BNULL$vocabulary == "Relapse" & as.integer(BNULL$null) == 3]),
  lab = c("nulls A, B", "null C"))
for (v in VOCABS) {
  z <- BNULL[BNULL$vocabulary == v, ]
  if (!(z$m[as.integer(z$null) == 3] > max(z$m[as.integer(z$null) %in% 1:2])))
    stop(sprintf("panel (b): %s null C no longer sits above nulls A and B; relabel", v))
}

set.seed(555)
p_b <- ggplot(B, aes(vocabulary, delta)) +
  geom_hline(yintercept = 0, colour = "#BFBFBF", linewidth = 0.3) +
  geom_segment(data = BNULL, aes(x = x0, xend = x1, y = m, yend = m),
               colour = "#8C8C8C", linewidth = 0.45, linetype = "11", inherit.aes = FALSE) +
  ## A, B sits just BELOW its line (the grey zero rule runs through it otherwise); C just above
  geom_text(data = BNLAB, aes(x = 2.40, y = y, label = lab, vjust = c(1.35, -0.35)), hjust = 0,
            size = 1.5, colour = "#8C8C8C", inherit.aes = FALSE) +
  geom_point(position = position_jitter(width = 0.10, height = 0), colour = CORR,
             size = 0.55, alpha = 0.75) +
  geom_segment(data = BSUM, aes(x = as.numeric(vocabulary) - 0.02,
                                xend = as.numeric(vocabulary) + 0.30,
                                y = mean_delta, yend = mean_delta, linetype = "cohort mean"),
               colour = INK, linewidth = 0.45, inherit.aes = FALSE) +
  geom_text(data = BSUM, aes(as.numeric(vocabulary) + 0.30, mean_delta,
                             label = fmt_nolead(mean_delta, 4, signed = TRUE)),
            hjust = 1, vjust = -0.7, size = 1.7, colour = INK, inherit.aes = FALSE) +
  scale_linetype_manual(values = c(`cohort mean` = "solid"), name = NULL,
                        guide = guide_legend(order = 1)) +
  scale_x_discrete(labels = c(Diagnosis = "Diagnosis
vocabulary",
                              Relapse = "Relapse
vocabulary"), expand = expansion(add = c(0.5, 0.85))) +
  scale_y_continuous(limits = BLIM, breaks = c(-0.05, 0, 0.05, 0.10, 0.15),
                     labels = c("-.05", "0", ".05", ".10", ".15"), expand = expansion(0)) +
  labs(x = NULL, y = "module change, Diagnosis to Relapse") +
  theme_tight() +
  theme(legend.position = "top", legend.direction = "vertical", legend.box = "vertical",
        legend.justification = "left", legend.box.just = "left",
        legend.box.margin = margin(0, 0, -3, -4), legend.spacing.y = unit(0.5, "pt"))

message("=== panel (b) built ===")

SUP <- jsonlite::fromJSON("analysis/out/f6_fig06_support.json", simplifyVector = FALSE)
MIX <- read.csv("analysis/out/f6_fig06_mixture.csv", stringsAsFactors = FALSE)
MIX <- MIX[MIX$vocabulary == "Diagnosis", ]          # the mixture is a property of the cells
MIX$stage <- factor(MIX$stage, levels = VOCABS)
ADM <- unlist(SUP$mixture_and_activity$admitted_celltypes)
N_ADM <- length(ADM)
N_ALL_CT <- as.numeric(SUP$mixture_and_activity$n_celltypes_in_the_substrate)
SPR <- SUP$mixture_and_activity$activity_spread_Diagnosis
MIX$celltype <- factor(MIX$celltype,
                       levels = MIX$celltype[MIX$stage == "Diagnosis"][
                         order(MIX$module_mean_AUC[MIX$stage == "Diagnosis"])])
message(sprintf("  [panel c] admitted %d of %d cell types; module activity at Diagnosis %.4f-%.4f (%.0f%% of mean)",
                N_ADM, N_ALL_CT, as.numeric(SPR$min), as.numeric(SPR$max),
                as.numeric(SPR$range_as_pct_of_mean)))
message(sprintf("  [panel c] largest mixture change: %s %.4f -> %.4f; reweighting alone moves the module %+.5f against a global activity fall of %+.5f",
                SUP$mixture_and_activity$largest_share_change$celltype,
                as.numeric(SUP$mixture_and_activity$largest_share_change$from),
                as.numeric(SUP$mixture_and_activity$largest_share_change$to),
                as.numeric(SUP$mixture_and_activity$reweighting_only_effect),
                as.numeric(SUP$mixture_and_activity$global_activity_fall)))

p_c_top <- ggplot(MIX, aes(module_mean_AUC, celltype)) +
  geom_line(aes(group = celltype), colour = "#BFBFBF", linewidth = 0.3) +
  geom_point(aes(size = cell_share, fill = stage), shape = 21, colour = "#FFFFFF", stroke = 0.25) +
  ## TIMEPOINTS BY LIGHTNESS, not in navy (2026-09-13): #0B3B5C is panel (e)'s "translation
  ## elongation", so Relapse in navy gave one colour two meanings in one figure. #99C7E0 is the
  ## light member of the set's blue family (Figures 2 and 3).
  scale_fill_manual(values = c(Diagnosis = CORR, Relapse = "#99C7E0"), name = NULL,
                    labels = c(Diagnosis = "at diagnosis", Relapse = "at relapse"),
                    guide = guide_legend(order = 1, override.aes = list(size = 1.8))) +
  scale_size_area(max_size = 4.2, breaks = c(0.05, 0.15, 0.30),
                  labels = c("5%", "15%", "30%"),
                  name = "share of malignant cells\n(six cell-types)",
                  ## The size keys were INVISIBLE (2026-10-21): the points are shape 21 with fill mapped
                  ## to `stage`, which is a different legend, so the size legend's keys had no fill and
                  ## a white outline on a white key. Give them a neutral grey fill so the dot area reads.
                  guide = guide_legend(order = 2, override.aes = list(fill = "#8C8C8C", colour = "#FFFFFF"))) +
  ## ANCHORED AT ZERO. On an axis that starts at zero the six cell types are visibly one tight
  ## band: that is the flatness, and it is the reason composition cannot move this module.
  scale_x_continuous(limits = c(0, 0.28), breaks = c(0, 0.1, 0.2),
                     labels = lab_nolead, expand = expansion(0)) +
  labs(x = "module activity within the cell-type (AUCell)",
       y = sprintf("admitted cell-types (%d of %d)", N_ADM, N_ALL_CT)) +
  theme_tight() +
  theme(axis.text.y = element_text(size = 5.0), panel.grid = element_blank(),
        legend.position = "top", legend.direction = "horizontal", legend.box = "vertical",
        legend.justification = "left", legend.box.just = "left",
        legend.title = element_text(size = 5.2),
        legend.box.margin = margin(0, 0, -3, -4), legend.spacing.y = unit(0.5, "pt"),
        plot.margin = margin(10, 4, 2, 3))

## the decomposition itself, as the identity naive = within-cell-type activity change + composition
DEC <- jsonlite::fromJSON("analysis/out/f6_decomposition.json", simplifyVector = FALSE)
DC <- do.call(rbind, lapply(VOCABS, function(v) {
  z <- DEC$module_decomposition$oriented_offsetcentred[[v]]
  data.frame(vocabulary = factor(v, levels = VOCABS),
             component = c("rewiring", "composition"),
             value = c(as.numeric(z$rewiring_mean), as.numeric(z$composition_mean)),
             naive = as.numeric(z$naive_mean),
             share = as.numeric(z$composition_share), stringsAsFactors = FALSE)
}))
DC$component <- factor(DC$component, levels = c("composition", "rewiring"))
for (v in VOCABS) {
  z <- DC[DC$vocabulary == v, ]
  if (abs(sum(z$value) - z$naive[1]) > 1e-12)
    stop(sprintf("panel (c): %s rewiring + composition does not equal naive", v))
}
message(sprintf("  [panel c] identity holds; %s",
                paste(sprintf("%s naive %+.6f = rewiring %+.6f + composition %+.6f (share %.4f)",
                              unique(DC$vocabulary),
                              tapply(DC$naive, DC$vocabulary, `[`, 1),
                              DC$value[DC$component == "rewiring"],
                              DC$value[DC$component == "composition"],
                              tapply(DC$share, DC$vocabulary, `[`, 1)), collapse = " | ")))

p_c_bot <- ggplot(DC, aes(value, vocabulary, fill = component)) +
  geom_col(width = 0.5, colour = NA) +
  geom_text(data = unique(DC[c("vocabulary", "naive")]),
            aes(naive, vocabulary, label = fmt_nolead(naive, 4, signed = TRUE)),
            hjust = -0.15, size = 1.6, colour = INK, inherit.aes = FALSE) +
  scale_y_discrete(labels = VOCLAB) +
  scale_fill_manual(values = c(rewiring = CORR, composition = "#C8C8C8"), name = NULL,
                    breaks = c("rewiring", "composition"),
                    labels = c("within-cell-type\nactivity change", "composition"),
                    guide = guide_legend(ncol = 2, reverse = FALSE)) +
  scale_x_continuous(limits = c(0, 0.0175), breaks = c(0, 0.005, 0.010, 0.015),
                     labels = c("0", ".005", ".010", ".015"), expand = expansion(0)) +
  labs(y = NULL, x = "module change, orientation-signed") +
  theme_tight() +
  theme(axis.text.y = element_text(size = 5.0), panel.grid = element_blank(),
        legend.position = "top", legend.direction = "horizontal",
        legend.justification = "left", legend.box.margin = margin(0, 0, -3, -4),
        plot.margin = margin(6, 4, 5, 3))

p_c <- p_c_top / p_c_bot + plot_layout(heights = c(1, 0.62))
message("=== panel (c) built ===")

AX <- jsonlite::fromJSON("analysis/out/marker_exclusion/results.json",
                         simplifyVector = FALSE)
AXLAB <- c(A1_lambo_HSC_minus_Monocyte = "A1  HSC - Monocyte",
           A2_lambo_HSCProg_minus_MonoGMP = "A2  HSC+Prog - Mono+GMP")
N_AX <- as.numeric(AX$checks$n_patients)
if (N_AX != N_EXT || !isTRUE(AX$checks$original_axis_and_all_four_regressions_reproduced))
  stop("panel (d): cohort or baseline reproduction check failed")
CC <- do.call(rbind, lapply(names(AXLAB), function(ax)
  do.call(rbind, lapply(VOCABS, function(v) {
    pick <- function(kind) {
      z <- Filter(function(x) x$axis == ax && x$vocabulary == v && x$construction == kind,
                  AX$summary)
      if (length(z) != 1L) stop("panel (d): expected one fit per axis, vocabulary and construction")
      z[[1]]
    }
    original <- pick("original_module_intact_axis")
    excluded <- pick("marker_excluded_module_intact_axis")
    data.frame(axis = unname(AXLAB[ax]), vocabulary = v,
               r2_original = as.numeric(original$r_squared),
               r2_excluded = as.numeric(excluded$r_squared), stringsAsFactors = FALSE)
  }))))
CC$axis <- factor(CC$axis, levels = rev(unname(AXLAB)))
CC$vocabulary <- factor(CC$vocabulary, levels = VOCABS)
CCL <- rbind(
  data.frame(CC[c("axis", "vocabulary")], r2 = CC$r2_original,
             what = "Original module", stringsAsFactors = FALSE),
  data.frame(CC[c("axis", "vocabulary")], r2 = CC$r2_excluded,
             what = "Maturity-marker targets removed", stringsAsFactors = FALSE))
CCL$what <- factor(CCL$what, levels = c("Original module", "Maturity-marker targets removed"))
message(sprintf("  [panel d] original R2 %.3f-%.3f; marker-excluded module R2 %.3f-%.3f; intact axes throughout",
                min(CC$r2_original), max(CC$r2_original), min(CC$r2_excluded), max(CC$r2_excluded)))

p_d <- ggplot(CCL, aes(r2, axis)) +
  geom_segment(data = CC, aes(x = r2_excluded, xend = r2_original, y = axis, yend = axis),
               colour = "#BFBFBF", linewidth = 0.45, inherit.aes = FALSE) +
  geom_point(aes(shape = what), fill = CORR, colour = CORR, size = 1.25, stroke = 0.4) +
  facet_wrap(~ vocabulary, nrow = 1, labeller = as_labeller(VOCLAB)) +
  scale_shape_manual(values = c(`Original module` = 21,
                                `Maturity-marker targets removed` = 23), name = NULL,
                     guide = guide_legend(ncol = 1,
                                          override.aes = list(fill = CORR, colour = CORR, size = 1.4))) +
  scale_x_continuous(limits = c(-0.03, 1.0), breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", ".25", ".5", ".75", "1"), expand = expansion(0)) +
  labs(y = NULL,
       x = sprintf("variance in module change explained (R\u00b2)\n(%d paired external patients)", N_AX)) +
  theme_tight() +
  theme(axis.text.y = element_text(size = 5.0), panel.grid = element_blank(),
        legend.position = "top", legend.direction = "vertical", legend.box = "vertical",
        legend.justification = "left", legend.box.just = "left",
        legend.box.margin = margin(0, 0, -3, -4), legend.spacing.y = unit(0.5, "pt"))
message("=== panel (d) built ===")

DEG <- read.csv("analysis/out/f6_unit_correction_degree.csv", stringsAsFactors = FALSE)
DEG <- DEG[DEG$unit_kind == "per_TF_s1", ]
DEG$s1 <- sub("^.*_", "", DEG$unit)

## TFs carrying BOTH arms inside the module -- the set where a per-TF merge would have occurred
BOTH <- sort(Reduce(intersect, lapply(c("+", "-"), function(s)
  unique(DEG$TF[DEG$s1 == s]))))
message(sprintf("  [panel e] %d TFs carry both arms in the module: %s",
                length(BOTH), paste(BOTH, collapse = ", ")))

## TWO PRE-DECLARED MARKER TERMS, chosen before looking and applied identically to every arm, so
## the classification is reproducible rather than an impression. Reactome, the library with power
## at this unit.
MK <- c(inflammatory = "Neutrophil Degranulation",
        translation  = "Eukaryotic Translation Elongation")
.keep <- c("stage", "unit_kind", "unit", "library", "term", "n_overlap", "odds_ratio", "q")
.hdr <- names(read.csv("analysis/out/f6_unit_correction_enrichment.csv", nrows = 1))
E <- read.csv("analysis/out/f6_unit_correction_enrichment.csv", stringsAsFactors = FALSE,
              colClasses = ifelse(.hdr %in% .keep, NA, "NULL"))
E <- E[E$unit_kind == "per_TF_s1" & E$library == "REAC" & E$term %in% MK, ]
message(sprintf("  [panel e] enrichment rows for the two marker terms at this unit: %d", nrow(E)))

D <- do.call(rbind, lapply(VOCABS, function(v) do.call(rbind, lapply(BOTH, function(tf)
  do.call(rbind, lapply(c("+", "-"), function(s) {
    u <- paste0(tf, "_", s)
    g <- DEG[DEG$vocabulary == v & DEG$unit == u, ]
    qi <- E$q[E$stage == v & E$unit == u & E$term == MK[["inflammatory"]]]
    qt <- E$q[E$stage == v & E$unit == u & E$term == MK[["translation"]]]
    qi <- if (length(qi)) qi else NA_real_; qt <- if (length(qt)) qt else NA_real_
    cls <- if (!is.na(qi) && qi < 0.05 && (is.na(qt) || qi <= qt)) "inflammatory"
           else if (!is.na(qt) && qt < 0.05) "translation" else "neither"
    data.frame(vocabulary = v, TF = tf, s1 = s, n_regulons = g$n_regulons,
               degree = g$degree, q_infl = qi, q_tran = qt, cls = cls,
               stringsAsFactors = FALSE)
  })))))) 
D$vocabulary <- factor(D$vocabulary, levels = VOCABS)
D$x <- ifelse(D$s1 == "-", -D$degree, D$degree)
CLSLAB <- c(inflammatory = "neutrophil degranulation (q < 0.05)",
            translation  = "translation elongation (q < 0.05)",
            neither      = "neither marker term at q < 0.05")
D$cls <- factor(unname(CLSLAB[D$cls]), levels = unname(CLSLAB))
## TF order: by total module targets at Diagnosis, so the hubs sit at the top
DORD <- names(sort(tapply(D$degree[D$vocabulary == "Diagnosis"],
                          D$TF[D$vocabulary == "Diagnosis"], sum)))
D$TF <- factor(D$TF, levels = DORD)
for (v in VOCABS)
  message(sprintf("  [panel e] %-9s %s", v, paste(sprintf(
    "%s %s%d:%s", D$TF[D$vocabulary == v], D$s1[D$vocabulary == v], D$degree[D$vocabulary == v],
    substr(as.character(D$cls[D$vocabulary == v]), 1, 4)), collapse = " ")))

## SPI1 at Diagnosis is the worked case the correction turns on: 1,039 = 729 + 310, and the two
## halves enrich for different things. Asserted rather than asserted-in-prose.
.sp <- D[D$vocabulary == "Diagnosis" & D$TF == "SPI1", ]
if (sum(.sp$degree) != 1039L) stop("panel (e): SPI1 Diagnosis arms do not sum to 1,039")
message(sprintf("  [panel e] SPI1 Diagnosis 1039 = +:%d + -:%d, classes %s",
                .sp$degree[.sp$s1 == "+"], .sp$degree[.sp$s1 == "-"],
                paste(substr(as.character(.sp$cls), 1, 12), collapse = " / ")))

ARM_VJUST <- 9.7    # vjust grows downward here: 8.4 puts the labels under the tick labels, above the x title
p_e <- ggplot(D, aes(x, TF, fill = cls)) +
  geom_col(width = 0.68, colour = NA) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.4) +
  geom_text(aes(label = degree, hjust = ifelse(s1 == "-", 1.18, -0.18)),
            size = 1.5, colour = INK) +
  facet_wrap(~ vocabulary, nrow = 1, labeller = as_labeller(VOCLAB)) +
  scale_fill_manual(values = setNames(c(CORR, "#0B3B5C", "#C8C8C8"), unname(CLSLAB)),
                    name = NULL, guide = guide_legend(ncol = 1)) +
  ## FEWER BREAKS AND MORE ROOM. At six breaks in a half-width facet the tick labels ran into one
  ## another ("4002000200400600800"), and at -460 the left bar labels were clipped by the facet
  ## edge -- both are data made unreadable, not styling.
  scale_x_continuous(limits = c(-620, 980), breaks = c(-400, 0, 400, 800),
                     labels = c("400", "0", "400", "800"), expand = expansion(0)) +
  ## The arm direction labels sit under EACH facet's own zero line (a facet_wrap has one shared x title,
  ## which centres on the whole plot and so landed between the two facets, aligned with neither).
  ## `target genes per arm` is the x title proper, on its own line beneath them.
  geom_text(data = data.frame(vocabulary = factor(VOCABS, levels = VOCABS)),
            aes(x = -15, y = -Inf, label = "negative (-)"),
            inherit.aes = FALSE, hjust = 1, vjust = ARM_VJUST, size = 1.55, colour = INK) +
  geom_text(data = data.frame(vocabulary = factor(VOCABS, levels = VOCABS)),
            aes(x = 15, y = -Inf, label = "positive (+)"),
            inherit.aes = FALSE, hjust = 0, vjust = ARM_VJUST, size = 1.55, colour = INK) +
  coord_cartesian(clip = "off") +
  labs(y = NULL,
       x = sprintf("target genes per arm
(module selected on %d paired discovery patients)", N_DISC)) +
  theme_tight() +
  theme(axis.text.y = element_text(size = 5.0), panel.grid = element_blank(),
        axis.title.x = element_text(size = 6.0, margin = margin(t = 13)),
        legend.position = "top", legend.direction = "vertical",
        legend.justification = "left", legend.box.just = "left",
        legend.box.margin = margin(0, 0, -3, -4))

message("=== panel (e) built ===")

ARMENR <- read.csv("analysis/out/f6_fig06_arm_enrichment.csv", stringsAsFactors = FALSE)
REPRO <- SUP$arm_enrichment$reproducible_across_vocabularies
SURV <- SUP$survives_ribosome_drop
ARMLAB <- c("activating (+)", "repressing (-)")

## TERM SELECTION, FIXED BEFORE LOOKING AT WHICH TERMS RESULT: within each arm, keep only terms
## significant in BOTH vocabularies, rank by the WORSE of the two q-values (so a term carried by
## one vocabulary cannot rise), take the top 6, and additionally force in every reproducible term
## that survives the ribosome drop, so the surviving minority on the negative-correlation target set cannot be
## crowded out by the terms that do not survive.
FULL <- ARMENR[ARMENR$gene_set == "full", ]
sel <- list()
for (a in ARMLAB) {
  rep_terms <- unlist(REPRO[[a]]$terms)
  s <- FULL[FULL$arm == a & FULL$term %in% rep_terms, ]
  worst <- aggregate(q ~ term + library, data = s, FUN = max)
  worst <- worst[order(worst$q), ]
  keep <- unique(c(head(worst$term, 6), intersect(unlist(SURV[[a]]$surviving_terms), rep_terms)))
  if (length(keep) > 10) keep <- unique(c(head(worst$term, 6),
                                          head(intersect(worst$term,
                                                         unlist(SURV[[a]]$surviving_terms)), 4)))
  sel[[a]] <- keep
  message(sprintf("  [panel f] %s: %d reproducible in both vocabularies, %d survive the ribosome drop, %d drawn",
                  a, length(rep_terms), length(unlist(SURV[[a]]$surviving_terms)), length(keep)))
}
FF <- do.call(rbind, lapply(ARMLAB, function(a)
  FULL[FULL$arm == a & FULL$term %in% sel[[a]], ]))
FF$arm <- factor(FF$arm, levels = ARMLAB)
FF$vocabulary <- factor(FF$vocabulary, levels = VOCABS)
FF$survives <- mapply(function(a, tm) tm %in% unlist(SURV[[a]]$surviving_terms),
                      as.character(FF$arm), FF$term)
FF$robust <- factor(ifelse(FF$survives,
                           "survives removal of ribosomal protein genes",
                           "lost when ribosomal protein genes are removed"),
                    levels = c("survives removal of ribosomal protein genes",
                               "lost when ribosomal protein genes are removed"))
## order terms by arm then by worst q, so the two arms read as blocks
.ord <- do.call(rbind, lapply(ARMLAB, function(a) {
  s <- FF[FF$arm == a, ]
  w <- aggregate(q ~ term, data = s, FUN = max)
  w$arm <- a
  w[order(-w$q), ]
}))
FF$term <- factor(FF$term, levels = .ord$term)
FF$neglog10q <- -log10(pmax(FF$q, 1e-300))
message(sprintf("  [panel f] %d dots drawn; background = the 201-regulon pool union (%d genes at Diagnosis)",
                nrow(FF), FF$n_background[1]))

p_f <- ggplot(FF, aes(vocabulary, term)) +
  geom_point(aes(size = n_overlap, fill = neglog10q, shape = robust),
             colour = INK, stroke = 0.28) +
  facet_grid(arm ~ ., scales = "free_y", space = "free_y", switch = "y",
             labeller = labeller(arm = setNames(c("positive TF-gene correlation", "negative TF-gene correlation"), ARMLAB))) +
  scale_fill_gradientn(colours = c("#DCE7F1", CORR, "#0B3B5C"),
                       name = expression(-log[10]~italic(q)),
                       guide = guide_colourbar(barwidth = unit(52, "pt"),
                                               barheight = unit(3.5, "pt"),
                                               title.position = "top", order = 1)) +
  scale_size_area(max_size = 3.4, breaks = c(50, 150, 300), name = "genes in common",
                  guide = guide_legend(order = 2)) +
  ## shape carries the ROBUSTNESS ASYMMETRY: a filled circle is a term that is still significant
  ## in both vocabularies after the 167 ribosomal protein genes are removed; a hollow diamond is
  ## one that is not. Without this the panel would present two findings of equal standing when
  ## one of them is fragile.
  scale_shape_manual(values = c(`survives removal of ribosomal protein genes` = 21,
                                `lost when ribosomal protein genes are removed` = 23),
                     labels = c("significant after ribosomal depletion", "not retained after ribosomal depletion"),
                     name = NULL,
                     guide = guide_legend(order = 3, ncol = 1,
                                          override.aes = list(size = 1.6, fill = "#FFFFFF"))) +
  ## THE TERM NAMES ARE FOLDED ONTO TWO LINES. Reactome names run to 63 characters here
  ## ("L13a-mediated Translational Silencing of Ceruloplasmin Expression"), and on one line they
  ## take so much of the panel's width that the dots are squeezed into a strip. Wrapping trades
  ## width for height, which is what this panel has to spare.
  scale_x_discrete(labels = setNames(sub(" ", "\n", VOCLAB), VOCABS)) +
  scale_y_discrete(labels = function(x)
    vapply(x, function(s) paste(strwrap(s, width = 38), collapse = "
"), character(1))) +
  labs(x = NULL,
       y = "enriched terms, against the 201-regulon pool background") +
  theme_tight() +
  theme(axis.text.y = element_text(size = 4.4, lineheight = 0.85),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 90, size = 5.2),
        panel.grid.major.y = element_line(colour = "#EEEEEE", linewidth = 0.2),
        legend.position = "bottom", legend.direction = "horizontal", legend.box = "vertical",
        legend.justification = "left", legend.box.just = "left",
        legend.title = element_text(size = 5.2),
        legend.box.margin = margin(-2, 0, 0, -4), legend.spacing.y = unit(0.5, "pt"))

message("=== panel (f) built ===")

## SUPPLEMENTARY: per-patient x per-regulon activity -- 19 paired discovery patients, DESCRIPTIVE
## DESCRIPTIVE, AND THE CAPTION SAYS SO. No test is run on this panel and no summary is drawn
## from it; it is the substrate every other discovery-cohort number rests on, shown whole.
PANELS <- list(a = p_a, b = p_b, c = p_c, d = p_d, e = p_f, f = p_e)
if (!identical(names(PANELS), c("a", "b", "c", "d", "e", "f")))
  stop("panel order is not a-f")
message(sprintf("  [order] panels %s", paste(names(PANELS), collapse = " ")))

RENDERED <- unlist(lapply(PANELS, function(pp) {
  parts <- if (inherits(pp, "patchwork")) pp$patches$plots else list(pp)
  parts <- c(parts, list(pp))
  unlist(lapply(parts, function(q) {
    if (is.null(q$labels)) return(NULL)
    l <- q$labels
    c(unlist(l[vapply(l, is.character, logical(1))]),
      ## get_labels() throws when a scale's breaks fall outside the data's range, and a GATE
      ## MUST NOT BE ABLE TO CRASH THE BUILD IT GUARDS -- a scale whose labels cannot be read is
      ## reported as unreadable, not silently skipped and not fatal.
      unlist(lapply(q$scales$scales, function(s)
        c(if (is.character(s$name)) s$name else NULL,
          tryCatch(as.character(s$get_labels()),
                   error = function(e) {
                     message("  [rule] NOTE: a scale's labels could not be read and were not ",
                             "checked (", conditionMessage(e), ")")
                     character(0)
                   })))),
      as.character(unlist(lapply(q$layers, function(ly)
        if (!is.null(ly$data) && is.data.frame(ly$data))
          unlist(ly$data[vapply(ly$data, is.character, logical(1))]) else NULL))))
  }))
}))
RENDERED <- RENDERED[!is.na(RENDERED) & nzchar(RENDERED)]
assert_no_retired_module(unique(RENDERED), "a rendered string")
message(sprintf("  [rule] PASS -- %d rendered strings checked, none names the retired module or ax4",
                length(unique(RENDERED))))

TAGTHEME <- theme(plot.tag = element_text(size = 9, face = "bold", family = "sans",
                                          colour = PAL[["truth"]], hjust = 0, vjust = 1),
                  plot.tag.position = c(0, 1.03))
BOTTOM <- wrap_elements(full = ((p_f + labs(tag = "e") + TAGTHEME) |
                                (p_e + labs(tag = "f") + TAGTHEME)) +
                          plot_layout(widths = c(1.08, 0.92)))
FIG <- (p_a | p_b) / (p_c | p_d) / BOTTOM +
  plot_layout(heights = c(1.00, 0.86, 2.5))
FIG <- FIG + plot_annotation(tag_levels = list(c("a", "b", "c", "", "d"))) &
  theme(plot.tag = element_text(size = 9, face = "bold", family = "sans",
                                colour = PAL[["truth"]], hjust = 0, vjust = 1),
        plot.tag.position = c(0, 1.03))

save_paper(FIG, OUT, w = COL2, h = 232 * MM)
mf_write_keys(OUT)

message("Figure 5 complete.")
