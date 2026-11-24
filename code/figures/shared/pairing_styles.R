
source("figures/shared/palette.R")

RULE_COL <- c(inference = unname(PAL[["inferred"]]),
              random    = unname(PAL[["baseline"]]),
              metacell  = unname(PAL[["context"]]),
              truth     = unname(PAL[["truth"]]))

RULE_FILL <- c(inference = unname(PAL_FILL[["inferred"]]),
               random    = unname(PAL_FILL[["baseline"]]),
               metacell  = unname(PAL_FILL[["context"]]),
               truth     = unname(PAL_FILL[["truth"]]))

RULE_SHAPE <- c(inference = 24, random = 21, metacell = 23, truth = 22)

RULE_NAME <- c(inference = "inference-assigned",
               random    = "random within stratum",
               metacell  = "metacell aggregation",
               truth     = "truth")

## Which rule each BMMC arm follows. Keyed on the arm token the rung report uses.
BM_RULE <- c(
  truth                             = "truth",
  link                              = "inference",
  link_1to1                         = "inference",
  link_x_celltype_donor             = "inference",
  link_x_celltype_donor_1to1        = "inference",
  link_x_celltype_donor_cap14       = "inference",
  rand_celltype_1to1                = "random",
  rand_celltype_x_donor_1to1        = "random",
  rand_celltype_donor_matched       = "random",
  rand_celltype_donor_matched_cap14 = "random",
  rand_global                       = "random",
  pseudomultiome                    = "metacell",
  pseudomultiome_k10                = "metacell"
)

## Short arm labels for the BMMC set, parallel in form so they read as a set.
BM_LABEL <- c(
  truth                             = "truth",
  link                              = "link (ArchR)",
  link_1to1                         = "link, 1-to-1",
  link_x_celltype_donor             = "inference, reuse 2.01x",
  link_x_celltype_donor_1to1        = "inference, 1-to-1",
  link_x_celltype_donor_cap14       = "inference, reuse 1.41x",
  rand_celltype_1to1                = "random, donor-blind",
  rand_celltype_x_donor_1to1        = "random, 1-to-1",
  rand_celltype_donor_matched       = "random, reuse 2.01x",
  rand_celltype_donor_matched_cap14 = "random, reuse 1.41x",
  rand_global                       = "random, no stratum",
  pseudomultiome                    = "metacell k = 5",
  pseudomultiome_k10                = "metacell k = 10"
)

rule_of  <- function(a) unname(BM_RULE[a])
rule_col <- function(a) unname(RULE_COL[BM_RULE[a]])
bm_label <- function(a) unname(BM_LABEL[a])

## The three concentration levels of the BMMC series, in ladder order. The arm tokens
## are read from the manifest, never typed, so a mis-paired level is impossible.
series_levels <- function() c("bijection", "cap14", "uncapped")

series_arm <- function(level, side)
  mf(sprintf("bmmc.series.%s.%s_arm", level, side), quiet = TRUE)

## Axis label for a concentration level: the measured reuse, read from the manifest.
series_reuse <- function(level)
  mf_num(sprintf("bmmc.reuse.%s", series_arm(level, "link")))

## Shared legend scale, so all figures name the assignment rules identically.
guide_rule <- function(keys = c("inference", "random")) {
  ggplot2::scale_colour_manual(
    values = RULE_COL[keys], breaks = keys, labels = RULE_NAME[keys], name = NULL)
}

tag_v2 <- function(p, levels = "a", top = 12) {
  p + patchwork::plot_annotation(tag_levels = levels) &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(size = 9, face = "bold", family = "sans",
                                       colour = PAL[["truth"]], hjust = 0, vjust = 1),
      plot.tag.position = c(0, 1.02),
      plot.margin = ggplot2::margin(top, 4, 2, 4))
}
