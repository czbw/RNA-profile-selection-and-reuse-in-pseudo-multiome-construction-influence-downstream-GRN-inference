#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ggplot2); library(grid); library(patchwork); library(scales); library(jsonlite)
})
source("figures/shared/theme_paper.R")
source("figures/shared/manifest.R")
source("figures/shared/pairing_styles.R")
source("figures/shared/number_labels.R")   # lab_nolead(): the set's no-leading-zero convention

OUT <- "figures/main/fig04_consequence"
message("\n=== Figure 4: the correction on the AML cohort (merged 4+5) ===")
INK <- "#111111"
STAGES <- c("Diagnosis", "Relapse")

ARMS4 <- data.frame(
  slot  = c("link", "repair", "pm"),
  label = c("published link", "corrected", "metacell k = 5"),
  fam   = c("inference", "corrected", "metacell"),
  col   = c("#D11F1F", "#4292C6", "#5C3696"),   # Fig2 inf_free / rnd_1to1 / met_k5
  key   = c("link", "repair", "pseudomultiome"),  # the token the manifest keys use
  stringsAsFactors = FALSE)
ARM_ORDER <- ARMS4$slot
PAL4  <- setNames(ARMS4$col,   ARMS4$slot)
LAB4  <- setNames(ARMS4$label, ARMS4$slot)
KEY4  <- setNames(ARMS4$key,   ARMS4$slot)

assert_arm_order <- function(d, panel, col = "slot") {
  got <- levels(d[[col]])
  want <- ARM_ORDER[ARM_ORDER %in% got]
  if (!identical(got, want))
    stop(sprintf("panel (%s): arm order is %s, want %s", panel,
                 paste(got, collapse = " > "), paste(want, collapse = " > ")))
  message(sprintf("  [order] panel (%s): PASS -- %s", panel, paste(got, collapse = " > ")))
}

theme_tight <- function() theme_paper() +
  theme(axis.title.x = element_text(size = 6.0, margin = margin(t = 1.5)),
        axis.title.y = element_text(size = 6.0, margin = margin(r = 1.5)),
        axis.text = element_text(size = 5.2),
        strip.text = element_text(size = 5.8, margin = margin(1.5, 1.5, 1.5, 1.5)),
        legend.box.spacing = unit(1.5, "pt"), legend.margin = margin(0, 0, 0, 0),
        legend.background = element_blank(), legend.key.size = unit(5.5, "pt"),
        legend.text = element_text(size = 5.2),
        plot.margin = margin(10, 4, 5, 3))

## n ON THE PANEL, not only in the caption, and read per panel because the populations in
## this figure are different.
N_DX <- mf_num("aml.cohort.Diagnosis.n_patients")
N_RL <- mf_num("aml.cohort.Relapse.n_patients")
N_S3 <- mf_num("s3r.cohort.n_paired_patients")
stage_lab_ab <- c(Diagnosis = "Diagnosis", Relapse = "Relapse")

A <- do.call(rbind, lapply(STAGES, function(s) data.frame(
  stage = s, slot = c("deposited", "repair"),
  v = c(mf_num(sprintf("aml.%s.old_reuse", s)), mf_num(sprintf("aml.%s.achieved_reuse", s))),
  stringsAsFactors = FALSE)))
A$stage <- factor(A$stage, levels = STAGES)
A$slot  <- factor(A$slot, levels = c("deposited", "repair"))
FL <- data.frame(stage = factor(STAGES, levels = STAGES),
                 y = vapply(STAGES, function(s) mf_num(sprintf("aml.floor.%s", s)), numeric(1)))

for (s in STAGES) {
  ach <- mf_num(sprintf("aml.%s.achieved_reuse", s)); flo <- mf_num(sprintf("aml.floor.%s", s))
  if (!(abs(ach - flo) < 1e-9))
    stop(sprintf("panel (a): %s achieved reuse %.4f does not equal its OWN floor %.4f", s, ach, flo))
  message(sprintf("  [panel a] %-9s deposited %.4fx -> corrected %.4fx, own floor %.4fx (imposed share of excess %.1f%%)",
                  s, mf_num(sprintf("aml.%s.old_reuse", s)), ach, flo,
                  100 * mf_num(sprintf("aml.excess.imposed_fraction.%s", s))))
}

p_a <- ggplot(A, aes(slot, v, fill = slot)) +
  geom_col(width = 0.62, colour = NA) +
  geom_hline(data = FL, aes(yintercept = y, linetype = "scarcity floor, this timepoint"),
             colour = INK, linewidth = 0.4, key_glyph = "path") +
  geom_text(aes(label = ifelse(slot == "repair", sprintf("%.4fx", v), sprintf("%.2fx", v))),
            vjust = -0.55, size = 1.75, colour = INK) +
  geom_text(data = FL, aes(x = 0.62, y = y, label = sprintf("floor %.4fx", y)),
            hjust = 0, vjust = -0.62, size = 1.6, colour = INK, inherit.aes = FALSE) +
  facet_wrap(~ stage, nrow = 1, labeller = as_labeller(stage_lab_ab)) +
  ## THE DEPOSITED PAIRING *IS* THE LINK ARM, so it takes the link colour. Giving it a
  ## neutral grey would break the one thing this figure's colour has to do: name the
  ## construction, identically, in every panel.
  scale_fill_manual(values = c(deposited = unname(PAL4[["link"]]),
                               repair = unname(PAL4[["repair"]])),
                    breaks = c("deposited", "repair"),
                    labels = c(deposited = unname(LAB4[["link"]]),
                               repair = unname(LAB4[["repair"]])), name = NULL) +
  scale_linetype_manual(values = c("scarcity floor, this timepoint" = "22"), name = NULL) +
  scale_y_log10(breaks = c(1, 2, 5, 10, 20), labels = c("1", "2", "5", "10", "20"),
                expand = expansion(mult = c(0.02, 0.16))) +
  labs(x = NULL, y = "reuse") +
  theme_tight() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "top", legend.direction = "vertical",
        legend.justification = "left", legend.box.margin = margin(0, 0, -3, -4))

DIRS <- data.frame(
  dir = c("repair_in_link", "link_in_repair"),
  lab = c("of corrected,\nalso called by published link", "of published link,\nsurviving correction"),
  slot = c("repair", "link"), stringsAsFactors = FALSE)
CC <- do.call(rbind, lapply(seq_len(nrow(DIRS)), function(i)
  do.call(rbind, lapply(STAGES, function(s) {
    num <- mf_num(sprintf("aml.contain.TFs.%s.%s.numerator", s, DIRS$dir[i]))
    den <- mf_num(sprintf("aml.contain.TFs.%s.%s.denominator", s, DIRS$dir[i]))
    data.frame(dir = DIRS$lab[i], slot = DIRS$slot[i], stage = s,
               num = num, den = den, v = num / den, stringsAsFactors = FALSE)
  }))))
CC$dir   <- factor(CC$dir, levels = DIRS$lab)
CC$stage <- factor(CC$stage, levels = STAGES)
CC$slot  <- factor(CC$slot, levels = ARM_ORDER)
message(sprintf("  [panel b] TF containment: %s",
                paste(sprintf("%s %s %d/%d = %.3f", CC$stage, sub("\n", " ", CC$dir),
                              CC$num, CC$den, CC$v), collapse = " | ")))

p_b <- ggplot(CC, aes(stage, v, fill = slot)) +
  geom_col(width = 0.62, colour = NA) +
  geom_text(aes(label = sprintf("%d/%d", num, den)), vjust = -0.5, size = 1.7, colour = INK) +
  facet_wrap(~ dir, nrow = 1) +
  scale_fill_manual(values = PAL4, breaks = ARM_ORDER, labels = LAB4[ARM_ORDER],
                    name = NULL, guide = "none") +
  scale_x_discrete(labels = c(Diagnosis = "Diagnosis",
                              Relapse = "Relapse")) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c("0", ".25", ".5", ".75", "1"), expand = expansion(0)) +
  labs(x = NULL, y = "fraction of that arm's TFs\nalso called by the other") +
  theme_tight()

FUN_AX <- data.frame(
  key = c("ax0_all", "ax1_directional", "ax2_selective", "ax3_celltype_robust"),
  lab = c("all shared regulons", "directional", "selective", "cell-type robust"),
  stringsAsFactors = FALSE)
FUN_AX$lab <- factor(FUN_AX$lab, levels = rev(FUN_AX$lab))

.pubjson <- "longitudinal/published/out/s3_nomination.json"
if (!file.exists(.pubjson)) stop("panel (c): missing published funnel source: ", .pubjson)
.pf <- jsonlite::fromJSON(.pubjson, simplifyVector = FALSE)$funnel
pub_n <- vapply(.pf, function(x) as.integer(x$n_regulons), integer(1))
if (length(pub_n) < nrow(FUN_AX))
  stop(sprintf("panel (c): published funnel has %d steps, need at least %d",
               length(pub_n), nrow(FUN_AX)))
## AXIS 5 IS DROPPED HERE, EXPLICITLY. The published funnel carries a step 5 and the corrected
## one does not, so the fifth entry is discarded rather than silently unused, and the build stops
## if the corrected arm ever gains one -- at which point this decision must be revisited.
if (mf_has("long.funnel.ax5_chip.n_regulons"))
  stop("panel (c): the corrected arm now has an axis-5 count; the drop-axis-5 decision must be re-derived")
message(sprintf("  [panel c] axis 5 dropped: published has %d steps, corrected has %d; the ChIP gate was never computed on the corrected arm",
                length(pub_n), nrow(FUN_AX)))
pub_n <- pub_n[seq_len(nrow(FUN_AX))]

FN <- do.call(rbind, lapply(seq_len(nrow(FUN_AX)), function(i) data.frame(
  step = FUN_AX$lab[i], key = FUN_AX$key[i],
  slot = c("link", "repair"),
  n = c(pub_n[i], mf_num(sprintf("long.funnel.%s.n_regulons", FUN_AX$key[i]))),
  stringsAsFactors = FALSE)))
FN$slot <- factor(FN$slot, levels = c("link", "repair"))
assert_arm_order(FN, "c")

M55  <- mf_num("fig4.set55.n_regulons")
AX3  <- mf_num("long.funnel.ax3_celltype_robust.n_regulons")
AX4  <- mf_num("long.funnel.ax4_tier1.n_regulons")
NREM <- mf_num("fig4.ax4_gate.n_would_be_removed")
if (M55 != AX3)
  stop(sprintf("panel (c): the corrected module (%d) is not the ax3 count (%d)", M55, AX3))
if (AX3 - AX4 != NREM)
  stop(sprintf("panel (c): ax3 - ax4 = %d but the recorded ax4 removal count is %d",
               AX3 - AX4, NREM))
GATE_CUT <- NREM / M55
message(sprintf("  [panel c] corrected module = %d (ax3); ax4 is an IMPORTED GATE that would remove %d of them (%.1f%%), leaving the retired %d",
                M55, NREM, 100 * GATE_CUT, AX4))
message(sprintf("  [panel c] published %s | corrected %s",
                paste(FN$n[FN$slot == "link"], collapse = " -> "),
                paste(FN$n[FN$slot == "repair"], collapse = " -> ")))

## EVERY BAR IS NOW AN APPLIED STEP on both arms -- tier 1, the only gate, is gone. Asserted,
## so a step reintroduced without revisiting this panel fails loudly instead of drawing as a
## solid bar that silently claims to have been applied.
if (any(FN$key == "ax4_tier1"))
  stop("panel (c): tier 1 is back in FUN_AX; the gate layer and its dashed key were removed")
FN$applied <- TRUE
## share of that arm's OWN pool
FN$share <- ave(FN$n, FN$slot, FUN = function(x) x / x[1])
FN$xdir  <- ifelse(FN$slot == "link", -FN$share, FN$share)
FN$tip <- as.character(FN$n)
FN_STEP <- FN

p_c <- ggplot(mapping = aes(xdir, step)) +
  geom_col(data = FN_STEP, aes(fill = slot), width = 0.62, colour = NA) +
  geom_text(data = FN, aes(label = tip, hjust = ifelse(slot == "link", 1.2, -0.18)),
            size = 1.7, colour = INK) +
  geom_vline(xintercept = 0, colour = INK, linewidth = 0.4) +
  scale_fill_manual(values = PAL4, breaks = c("link", "repair"),
                    labels = LAB4[c("link", "repair")], name = NULL,
                    guide = guide_legend(order = 1)) +
  scale_colour_manual(values = PAL4, guide = "none") +
  ## PADDING IS FOR THE TIP LABELS AND NOTHING ELSE. At +-1.46 the two arms were pushed into
  ## the middle third of a full-width row and the panel read as mostly empty. The widest tip
  ## label is the gate's, which starts at 0.169 and needs about 0.25 of the axis; 1.28 covers
  ## it with the bars still reaching the full width they are entitled to.
  scale_x_continuous(limits = c(-1.28, 1.28),
                     breaks = c(-1, -0.5, 0, 0.5, 1), labels = c("1", ".5", "0", ".5", "1"),
                     expand = expansion(0)) +
  labs(x = "share of that arm's own starting pool", y = NULL) +
  theme_tight() +
  theme(axis.text.y = element_text(size = 5.2), panel.grid = element_blank(),
        legend.position = "top", legend.direction = "horizontal",
        legend.justification = "left", legend.box.margin = margin(0, 0, -3, -4))

D_SETS <- data.frame(
  set = c("corrected module", "LMM-significant"),
  n_reg = c(mf_num("fig4.set55.n_regulons"), mf_num("fig5.lmm26.n_regulons")),
  n_tf  = c(mf_num("fig4.set55.n_TFs"),      mf_num("fig5.lmm26.n_TFs")),
  ieg   = c(mf_num("fig4.set55.ieg.n_in_set55"), mf_num("fig4.set55.ieg.n_in_lmm26")),
  stringsAsFactors = FALSE)
D_SETS$set <- factor(D_SETS$set, levels = c("corrected module", "LMM-significant"))
D_SETS$frac_ieg <- D_SETS$ieg / D_SETS$n_reg
## COLOUR (2026-09-13). Both sets are the CORRECTED arm's objects -- the 26 are `repair_LMM26` in
## figures/data/build_manifest.py, the corrected arm's E1 LMM output; the published arm's E1 was never
## fitted (panel e's absence gate). Red means the published arm in every other panel of the set, so
## the 26 were drawn in the wrong arm's colour. They take the corrected arm's dark second tone, the
## same #0B3B5C Figure 5 uses for a second category inside the corrected arm.
LMM_COL <- "#0B3B5C"
D_INT <- mf_num("fig4.set55.intersection.n_regulons")
D_JAC <- mf_num("fig4.set55.intersection.jaccard_regulon")
## the two IEG counts are drawn from ONE frozen list, and the panel says so by asserting it
message(sprintf("  [panel d] module %d/%d vs LMM %d/%d | intersection %d regulons, Jaccard %.4f | AP-1/IEG %d/%d vs %d/%d (frozen list of %d TFs)",
                D_SETS$n_reg[1], D_SETS$n_tf[1], D_SETS$n_reg[2], D_SETS$n_tf[2],
                D_INT, D_JAC, D_SETS$ieg[1], D_SETS$n_reg[1], D_SETS$ieg[2], D_SETS$n_reg[2],
                mf_num("fig4.set55.ieg.list_size")))

## the set relationship as three disjoint counts, so the intersection is a drawn quantity
DV <- data.frame(
  part = factor(c("module only", "in both", "LMM only"),
                levels = c("module only", "in both", "LMM only")),
  n = c(D_SETS$n_reg[1] - D_INT, D_INT, D_SETS$n_reg[2] - D_INT))
DV$fam <- factor(c("corrected module", "in both", "LMM-significant"),
                 levels = c("corrected module", "in both", "LMM-significant"))

p_d <- ggplot(DV, aes(part, n, fill = fam)) +
  geom_col(width = 0.62, colour = NA) +
  geom_text(aes(label = n), vjust = -0.45, size = 1.75, colour = INK) +
  scale_fill_manual(values = c(`corrected module` = unname(PAL4[["repair"]]),
                               `in both` = "#8C8C8C",
                               `LMM-significant` = LMM_COL),
                    breaks = c("corrected module", "in both", "LMM-significant"), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  labs(x = NULL, y = "eRegulons") +
  theme_tight() +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.justification = "left", legend.box.margin = margin(0, 0, -3, -4),
        plot.margin = margin(10, 4, 2, 3))

CAL <- do.call(rbind, lapply(STAGES, function(s) data.frame(
  stage = s, slot = c("link", "repair"),
  v = c(mf_num(sprintf("s3.calib.E1.%s.frac_p_below_05", s)),
        mf_num(sprintf("s3r.calib.E1.%s.frac_p_below_05", s))),
  ks = c(mf_num(sprintf("s3.calib.E1.%s.ks_uniform_p", s)),
         mf_num(sprintf("s3r.calib.E1.%s.ks_uniform_p", s))),
  stringsAsFactors = FALSE)))
CAL$stage <- factor(CAL$stage, levels = STAGES)
CAL$slot  <- factor(CAL$slot, levels = ARM_ORDER)
ALPHA <- mf_num("s3.calib.nominal_alpha")
CAL$fails_ks <- CAL$ks < 0.05
message(sprintf("  [panel e] calibration: %s",
                paste(sprintf("%s %s %.4f (KS %.3g%s)", CAL$stage, CAL$slot, CAL$v, CAL$ks,
                              ifelse(CAL$fails_ks, ", FAILS uniformity", "")),
                      collapse = " | ")))

## THE ABSENCE IS ASSERTED FROM THE MANIFEST, not inferred from a missing key.
if (!identical(mf("s3.e1.results_block_is_read_false", quiet = TRUE), TRUE))
  stop("panel (e): s3.e1.results_block_is_read_false is not TRUE; the link arm's absence must be re-derived")
if (!identical(mf("s3.e1.fits_file_present", quiet = TRUE), FALSE))
  stop("panel (e): s3.e1.fits_file_present is not FALSE; the link arm may in fact have fits")
if (!identical(mf_entry("s3.e1.confirmed_regulons.n")$status, "NA"))
  stop("panel (e): s3.e1.confirmed_regulons.n is not marked NOT APPLICABLE")
message("  [panel e] absence gate: PASS -- link E1 block is {read:false}, no fits file, count marked NOT APPLICABLE")

p_e <- ggplot(CAL, aes(stage, v, fill = slot)) +
  geom_hline(aes(yintercept = ALPHA, linetype = "nominal 0.05"), colour = INK,
             linewidth = 0.4, key_glyph = "path") +
  geom_col(width = 0.62, position = position_dodge(width = 0.70), colour = NA) +
  geom_text(aes(label = sub("^0", "", sprintf("%.4f", v))),
            position = position_dodge(width = 0.70), vjust = -1.15, size = 1.6, colour = INK) +
  scale_fill_manual(values = PAL4, breaks = c("link", "repair"),
                    labels = LAB4[c("link", "repair")], name = NULL,
                    guide = guide_legend(order = 1)) +
  scale_linetype_manual(values = c("nominal 0.05" = "22"), name = NULL,
                        guide = guide_legend(order = 2)) +
  scale_x_discrete(labels = c(Diagnosis = "Diagnosis",
                              Relapse = "Relapse")) +
  scale_y_continuous(labels = lab_nolead, expand = expansion(mult = c(0, 0.30))) +
  labs(x = NULL, y = "false positives at p < 0.05,\nscrambled labels") +
  theme_tight() +
  theme(legend.position = "top", legend.direction = "horizontal", legend.box = "vertical",
        legend.justification = "left", legend.box.just = "left",
        legend.box.margin = margin(0, 0, -3, -4),
        legend.spacing.y = unit(0.5, "pt"), legend.spacing.x = unit(4, "pt"))

## THE DROPPED PANEL'S EVIDENCE, PRINTED SO THE MERGE CANNOT SILENTLY LOSE IT.
## Former 5d is gone from the artwork; the caption and the text must carry it.
message(sprintf("  [merge] EXTERNAL GATE 5 IS NO LONGER DRAWN. It must appear in the caption: %s",
                paste(vapply(STAGES, function(v) sprintf(
                  "%s %+.4f, %.1f%% of %d patients, Wilcoxon p %.3g", v,
                  mf_num(sprintf("long.gate5.%s.mean_delta", v)),
                  100 * mf_num(sprintf("long.gate5.%s.frac_positive", v)),
                  mf_num("long.gate5.n_paired_patients"),
                  mf_num(sprintf("long.gate5.%s.p_wilcoxon", v))), character(1)),
                  collapse = " | ")))



FIG <- (p_a | p_b) / p_c / (p_d | p_e) +
  plot_layout(heights = c(1, 0.95, 0.80))
FIG <- FIG + plot_annotation(tag_levels = list(c("a", "b", "c", "d", "e"))) &
  theme(plot.tag = element_text(size = 9, face = "bold", family = "sans",
                                colour = PAL[["truth"]], hjust = 0, vjust = 1),
        plot.tag.position = c(0, 1.03))

save_paper(FIG, OUT, w = COL2, h = 158 * MM)
mf_write_keys(OUT)
message("=== Figure 4 done ===\n")
