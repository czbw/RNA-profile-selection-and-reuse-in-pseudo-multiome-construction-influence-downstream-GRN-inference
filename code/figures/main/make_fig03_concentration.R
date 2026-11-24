#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2); library(grid); library(patchwork); library(scales)
})
source("figures/shared/theme_paper.R")
source("figures/shared/manifest.R")
source("figures/shared/number_labels.R")   # lab_nolead(): the set's no-leading-zero convention
source("figures/shared/pairing_styles.R")

OUT <- "figures/main/fig03_concentration"
message("\n=== Figure 3: the concentration series (BMMC), v3 ===")

INK    <- "#111111"
LEVELS <- series_levels()                       # bijection, cap14, uncapped

## unname() is load-bearing -- series_arm() returns a NAMED character and the name leaks
## into every factor level built from it. This bit v2 and the note is kept.
link_arm <- unname(vapply(LEVELS, series_arm, character(1), side = "link"))
rand_arm <- unname(vapply(LEVELS, series_arm, character(1), side = "rand"))
REUSE    <- unname(vapply(LEVELS, series_reuse, numeric(1)))
names(link_arm) <- names(rand_arm) <- names(REUSE) <- LEVELS

## THE PALETTE. Hue = rule (Figure 2's own families), lightness = concentration.
## Hexes carried over from Figure 2 EXACTLY are marked; the rest are picked off the same
## ColorBrewer ramps at a target L*, which is how Figure 2's own palette was built.
PAL_F3 <- c(
  inference.bijection = "#66000C",   # L* 19  == Fig2 inf_matched      (SAME ARM)
  inference.cap14     = "#C0151A",   # L* 41  new -- Reds ramp
  inference.uncapped  = "#FB6B4B",   # L* 63  == Fig2 inf_r201         (SAME ARM)
  random.bijection    = "#0F5AA3",   # L* 38  new -- Blues ramp, see the note above
  random.cap14        = "#4292C6",   # L* 58  Fig2 hex, different level -- see README
  random.uncapped     = "#99C7E0")   # L* 78  == Fig2 rnd_r201         (SAME ARM)
RULE_HUE <- c(inference = "#D11F1F", random = "#2E7EBB")   # connecting lines only

## ACCESSIBILITY GATE -- recomputed on every build, printed, never asserted from memory.
## Same floor and same five vision modes as Figure 2. Within-family is what governs:
## between families, hue AND shape are both available.
palette_check <- function(P) {
  ok <- requireNamespace("colorspace", quietly = TRUE) &&
        requireNamespace("farver", quietly = TRUE)
  if (!ok) { message("  [a11y] colorspace/farver absent -- CHECK NOT RUN"); return(invisible()) }
  lab <- function(h) farver::convert_colour(t(grDevices::col2rgb(h)), "rgb", "lab")
  dE  <- function(a, b) farver::compare_colour(lab(a), lab(b), "lab", method = "cie2000")[1, 1]
  fam <- sub("\\..*$", "", names(P))
  message(sprintf("  [a11y] L* by arm: %s",
                  paste(sprintf("%s=%.0f", names(P), lab(P)[, 1]), collapse = "  ")))
  worst <- Inf
  for (mode in c("normal", "deutan", "protan", "tritan", "greyscale")) {
    Q <- switch(mode, normal = P, deutan = colorspace::deutan(P),
                protan = colorspace::protan(P), tritan = colorspace::tritan(P),
                greyscale = colorspace::desaturate(P))
    wi <- Inf; wp <- ""
    for (i in seq_len(length(Q) - 1)) for (j in (i + 1):length(Q))
      if (fam[i] == fam[j]) {
        d <- dE(Q[i], Q[j])
        if (d < wi) { wi <- d; wp <- paste(names(P)[i], names(P)[j], sep = "/") }
      }
    worst <- min(worst, wi)
    message(sprintf("  [a11y] %-10s within-family %5.1f  (%s)", mode, wi, wp))
  }
  if (worst < 10)
    stop(sprintf("palette gate FAILED: within-family separation falls to %.1f dE2000", worst))
  message(sprintf("  [a11y] GATE PASS -- worst within-family separation %.1f dE2000 over all five modes",
                  worst))
}
palette_check(PAL_F3)

key <- function(rule, level) paste(rule, level, sep = ".")

## AML'S OPERATING RANGE. Read, never typed. Drawn on all four panels, to scale.
AML_LO <- min(mf_num("aml.Diagnosis.achieved_reuse"), mf_num("aml.Relapse.achieved_reuse"))
AML_HI <- max(mf_num("aml.Diagnosis.achieved_reuse"), mf_num("aml.Relapse.achieved_reuse"))
message(sprintf("  [aml] operating range %.4fx - %.4fx; the cap14 pair sits at %.4fx",
                AML_LO, AML_HI, REUSE[["cap14"]]))

aml_band <- function() annotate("rect", xmin = AML_LO, xmax = AML_HI,
                                ymin = -Inf, ymax = Inf,
                                fill = PAL[["context"]], alpha = 0.16)

## the shared x scale -- ONE definition, so the AML band lands identically in all four
## WHERE THE 1.00x POINT SITS INSIDE THE PANEL, in npc. Both inside legends are anchored to
## it, so a legend's left edge lines up with the first datum rather than with the panel edge or
## the axis. Derived from the scale rather than typed: with a symmetric multiplicative expansion
## of m, the first datum sits at m/(1+2m) of the panel width whatever the data are.
X_PAD  <- 0.10
X_NPC  <- X_PAD / (1 + 2 * X_PAD)
message(sprintf("  [layout] the 1.00x datum sits at %.4f of panel width; both legends anchor there",
                X_NPC))

X_REUSE <- function(...) scale_x_continuous(
  breaks = REUSE, labels = sprintf("%.2fx", REUSE),
  expand = expansion(mult = c(X_PAD, X_PAD)), ...)

theme_tight <- function() theme_paper() +
  theme(axis.title.x = element_text(size = 6.0, margin = margin(t = 1.5)),
        axis.title.y = element_text(size = 6.0, margin = margin(r = 1.5)),
        axis.text = element_text(size = 5.4),
        legend.box.spacing = unit(1.5, "pt"), legend.margin = margin(0, 0, 0, 0),
        legend.background = element_blank(), legend.key.size = unit(5.5, "pt"),
        legend.text = element_text(size = 5.4),
        plot.margin = margin(10, 4, 5, 3))

## rows of (rule, level) with the arm and the shared plotting aesthetics attached
grid_arms <- function() do.call(rbind, lapply(LEVELS, function(l)
  do.call(rbind, lapply(c("inference", "random"), function(r) data.frame(
    level = l, rule = r, arm = if (r == "inference") link_arm[[l]] else rand_arm[[l]],
    reuse = REUSE[[l]], k = key(r, l), stringsAsFactors = FALSE)))))

## (a) the inference really does identify cells -- 5-6x chance at every level
CHANCE <- mf_num("bmmc.exact.chance")
A <- grid_arms()
A$v    <- vapply(A$arm, function(a) mf_num(sprintf("bmmc.exact.%s", a)), numeric(1))
A$lift <- A$v / CHANCE
message(sprintf("  [panel a] lift over chance: inference %s | random %s",
                paste(sprintf("%.1fx", A$lift[A$rule == "inference"]), collapse = " "),
                paste(sprintf("%.1fx", A$lift[A$rule == "random"]), collapse = " ")))

## the hand-built key for (a): first entry ON the 1.00x datum, both at one y just under the
## panel top. Positions are derived from the data range so they cannot fall out of step with it.
LEG_DX <- 0.030 * diff(range(REUSE))
LEG_A <- data.frame(
  rule = c("inference", "random"),
  lab  = c("inference-assigned", "random draw"),
  fill = unname(RULE_HUE[c("inference", "random")]),
  x    = REUSE[["bijection"]] + c(0, 0.60) * diff(range(REUSE)),
  y    = max(A$v) + 0.30 * diff(range(A$v)),
  stringsAsFactors = FALSE)

## MARKER OUTLINES ARE WHITE (2026-09-13), not the #111111 of the 2026-08-28 build. A dark ring made
## this the only figure in the set with outlined points. The outline cannot simply be dropped or
## matched to the fill, though: markers of the two rules OVERLAP here (both bijections at 1.000 in
## (d)), and without a separator qc_fig02.py finds runs holding two arm colours -- 23 of them on the
## fill-matched trial render. White keeps the marks flat and keeps them apart; Figure 5(c) uses it.
p_a <- ggplot(A, aes(reuse, v, group = rule)) +
  aml_band() +
  geom_hline(yintercept = CHANCE, linetype = "22", colour = PAL[["context"]],
             linewidth = 0.3) +
  ## a dashed rule with no name is a rule the reader has to guess at; this is a label
  ## attached to a mark, not a sentence
  annotate("text", x = min(REUSE), y = CHANCE, label = "chance", hjust = 0, vjust = 1.9,
           size = 1.75, colour = PAL[["context"]]) +
  geom_line(aes(colour = rule), linewidth = 0.4, show.legend = FALSE) +
  geom_point(aes(fill = k, shape = rule), colour = "#FFFFFF", size = 1.7, stroke = 0.30) +
  geom_text(aes(label = sprintf("%.1fx", lift)), vjust = -1.15, size = 1.75, colour = INK) +
  scale_colour_manual(values = RULE_HUE, guide = "none") +
  scale_fill_manual(values = PAL_F3, guide = "none") +
  scale_shape_manual(values = c(inference = 24, random = 21), guide = "none") +
  geom_point(data = LEG_A, aes(x, y, shape = rule), fill = LEG_A$fill, colour = "#FFFFFF",
             size = 1.7, stroke = 0.30, inherit.aes = FALSE) +
  geom_text(data = LEG_A, aes(x + LEG_DX, y, label = lab), hjust = 0, size = 1.85,
            colour = INK, inherit.aes = FALSE) +
  X_REUSE() +
  scale_y_continuous(breaks = seq(0.0025, 0.02, by = 0.0025),
                     labels = function(x) lab_nolead(round(x * 1000, 6)),
                     expand = expansion(mult = c(0.14, 0.10))) +
  ## the "5.7x / 1.0x" point labels are LIFT OVER CHANCE -- ratios -- and are unaffected by the
  ## axis unit; they are deliberately left exactly as they were.
  labs(x = "reuse", y = "cells given their own\ntranscriptome, per 1,000") +
  theme_tight() +
  theme(legend.position = "none")

## (b) the topology advantage INVERTS with concentration -- one line per rung
RUNGS <- data.frame(
  stem = c("bmmc.L1.rho", "bmmc.L1.topn_jaccard", "bmmc.L2.topn_jaccard",
           "bmmc.L3.genes", "bmmc.L3.regions", "bmmc.L4.jaccard", "bmmc.L5.jaccard"),
  ## terse on purpose: the rung vocabulary is defined once, in Figure 2, and repeated in
  ## this figure's caption.
  lab  = c("L1 rho", "L1 J", "L2 J", "L3 genes", "L3 regions", "L4", "L5"),
  stringsAsFactors = FALSE)

B <- do.call(rbind, lapply(LEVELS, function(l)
  do.call(rbind, lapply(seq_len(nrow(RUNGS)), function(j) data.frame(
    level = l, reuse = REUSE[[l]], rung = RUNGS$lab[j],
    delta = mf_num(sprintf("%s.%s", RUNGS$stem[j], link_arm[[l]])) -
            mf_num(sprintf("%s.%s", RUNGS$stem[j], rand_arm[[l]])),
    stringsAsFactors = FALSE)))))
B$sign <- ifelse(B$delta >= 0, "inference", "random")

## THE SIGN GATE. The panel's whole claim is that the sign flips, so it is asserted from
## the manifest rather than left to the eye: every rung positive at the bijection, every
## rung negative at 2.01x. If a rebuild ever breaks this the figure must not be drawn.
.b1 <- B$delta[B$level == "bijection"]; .b3 <- B$delta[B$level == "uncapped"]
if (!(all(.b1 > 0) && all(.b3 < 0)))
  stop(sprintf(paste("panel (b): the sign flip does not hold. %d/%d rungs positive at the",
                     "bijection, %d/%d negative at 2.01x -- the panel's claim has changed"),
               sum(.b1 > 0), length(.b1), sum(.b3 < 0), length(.b3)))
message(sprintf("  [panel b] sign gate: PASS -- all %d rungs positive at %.2fx, all %d negative at %.2fx",
                length(.b1), REUSE[["bijection"]], length(.b3), REUSE[["uncapped"]]))

LAB_GAP <- 0.0062
lab_pos <- function(d) {
  o <- order(-d$delta)
  y <- d$delta[o]
  for (i in seq_along(y)[-1]) y[i] <- min(y[i], y[i - 1] - LAB_GAP)
  d$laby[o] <- y
  d
}
BL <- lab_pos(B[B$level == "uncapped", ])
message(sprintf("  [panel b] label spread: max shift %.4f in y units (%d of %d labels moved)",
                max(abs(BL$laby - BL$delta)), sum(abs(BL$laby - BL$delta) > 1e-9), nrow(BL)))

LABX <- max(REUSE) + 0.055
p_b <- ggplot(B, aes(reuse, delta, group = rung)) +
  aml_band() +
  annotate("text", x = (AML_LO + AML_HI) / 2, y = Inf, label = "AML", angle = 90,
           hjust = 1.10, vjust = -0.40, size = 1.75, colour = PAL[["context"]]) +
  geom_line(colour = PAL[["context"]], linewidth = 0.30) +
  ## the zero line is heavier than the data lines: it is the thing being crossed
  geom_hline(yintercept = 0, colour = INK, linewidth = 0.45) +
  geom_point(aes(fill = sign), shape = 21, colour = "#FFFFFF", size = 1.25, stroke = 0.28) +
  geom_segment(data = BL, aes(x = reuse, xend = LABX - 0.008, y = delta, yend = laby),
               colour = PAL[["context"]], linewidth = 0.18) +
  geom_text(data = BL, aes(x = LABX, y = laby, label = rung),
            hjust = 0, size = 1.6, colour = INK) +
  scale_fill_manual(values = c(inference = RULE_HUE[["inference"]],
                               random = RULE_HUE[["random"]]), guide = "none") +
  scale_x_continuous(breaks = REUSE, labels = sprintf("%.2fx", REUSE),
                     expand = expansion(mult = c(0.10, 0.34))) +
  scale_y_continuous(labels = lab_nolead, expand = expansion(mult = c(0.10, 0.10))) +
  labs(x = "reuse", y = "inference - random\n(topology rungs)") +
  theme_tight()

C <- grid_arms()
C$v <- vapply(C$arm, function(a) mf_num(sprintf("bmmc.stratum_mean.rho.%s", a)), numeric(1))
C$lab <- sub("^0", "", sprintf("%.3f", C$v))
C$vj  <- ifelse(C$rule == "inference", -1.05, 1.85)
message(sprintf("  [panel c] medians: inference %s | random %s",
                paste(sprintf("%.4f", C$v[C$rule == "inference"]), collapse = " "),
                paste(sprintf("%.4f", C$v[C$rule == "random"]), collapse = " ")))

p_c <- ggplot(C, aes(reuse, v, group = rule)) +
  aml_band() +
  geom_hline(yintercept = 1, colour = INK, linewidth = 0.45) +
  geom_line(aes(colour = rule), linewidth = 0.4, show.legend = FALSE) +
  geom_point(aes(fill = k, shape = rule), colour = "#FFFFFF", size = 1.7, stroke = 0.30,
             show.legend = FALSE) +
  geom_text(aes(label = lab, vjust = vj), size = 1.75, colour = INK, show.legend = FALSE) +
  scale_colour_manual(values = RULE_HUE, guide = "none") +
  scale_fill_manual(values = PAL_F3, guide = "none") +
  scale_shape_manual(values = c(inference = 24, random = 21), guide = "none") +
  X_REUSE() +
  scale_y_continuous(labels = lab_nolead, expand = expansion(mult = c(0.12, 0.20))) +
  labs(x = "reuse", y = "stratum-mean rho across donors\n(median over cell-types, regulons)") +
  theme_tight()

## (d) BOTH variance components. Within falls; total rises. The relocation IS the panel.
D <- do.call(rbind, lapply(c("within", "total"), function(cp) {
  g <- grid_arms(); g$component <- cp
  g$v <- vapply(g$arm, function(a) mf_num(sprintf("bmmc.var.%s.%s", cp, a)), numeric(1))
  g
}))
D$component <- factor(D$component, levels = c("within", "total"))
message(sprintf("  [panel d] within %s | total %s",
                paste(sprintf("%.4f", D$v[D$component == "within"]), collapse = " "),
                paste(sprintf("%.4f", D$v[D$component == "total"]), collapse = " ")))

## THE BIJECTION IS AN IDENTITY, NOT A MEASUREMENT, and the panel is only readable if that
## is true: a permutation cannot change a multiset, so both bijective arms must return
## exactly 1. This is the estimator's own positive control and it is asserted, not hoped for.
.bij <- D$v[D$level == "bijection"]
if (any(abs(.bij - 1) > 1e-9))
  stop(sprintf("panel (d): a bijective arm returns %.12f, not 1 -- the estimator's positive control fails",
               .bij[which.max(abs(.bij - 1))]))
message("  [panel d] positive control: PASS -- both bijections return 1.0000 exactly")

p_d <- ggplot(D, aes(reuse, v, group = interaction(rule, component))) +
  aml_band() +
  geom_hline(yintercept = 1, colour = INK, linewidth = 0.45) +
  geom_line(aes(colour = rule, linetype = component), linewidth = 0.4) +
  geom_point(aes(fill = k, shape = rule), colour = "#FFFFFF", size = 1.7, stroke = 0.30,
             show.legend = FALSE) +
  scale_colour_manual(values = RULE_HUE, guide = "none") +
  scale_fill_manual(values = PAL_F3, guide = "none") +
  scale_shape_manual(values = c(inference = 24, random = 21), guide = "none") +
  scale_linetype_manual(values = c(within = "solid", total = "42"), name = NULL,
                        labels = c(within = "within-cell-type variance",
                                   total = "total variance"),
                        guide = guide_legend(order = 1, keywidth = unit(16, "pt"),
                                             override.aes = list(colour = INK,
                                                                 linewidth = 0.45))) +
  X_REUSE() +
  ## headroom for the legend, which sits INSIDE the panel
  scale_y_continuous(labels = lab_nolead, expand = expansion(mult = c(0.12, 0.36))) +
  labs(x = "reuse", y = "variance retained,\nfraction of native reference") +
  theme_tight() +
  ## INSIDE the panel and anchored on the first datum, exactly as (a). Drawn on TOP with
  ## justification "left" the legend aligned to the left edge of the PLOT, which put the solid
  ## key straight on top of this panel's own y-axis line: the two rendered as a single T and the
  ## key stopped reading as a key at all. Anchoring at X_NPC moves it clear of the axis.
  theme(legend.position = "inside",
        legend.position.inside = c(X_NPC, 1),
        legend.justification.inside = c(0, 1),
        legend.direction = "horizontal", legend.margin = margin(0, 0, 0, 0),
        legend.key.spacing.x = unit(7, "pt"))



FIG <- (p_a | p_b) / (p_c | p_d) +
  plot_layout(widths = c(1, 1.20), heights = c(1, 1))
FIG <- FIG + plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold", family = "sans",
                                colour = PAL[["truth"]], hjust = 0, vjust = 1),
        plot.tag.position = c(0, 1.03))

save_paper(FIG, OUT, w = COL2, h = 108 * MM)
mf_write_keys(OUT)
message("=== Figure 3 done ===\n")
