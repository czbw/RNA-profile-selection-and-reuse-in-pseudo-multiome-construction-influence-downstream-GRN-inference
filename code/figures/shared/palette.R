
PAL <- c(
  inferred = "#B2182B",
  baseline = "#2166AC",
  context  = "#8C8C8C",
  truth    = "#111111"
)

## pale fills for boxes and bars that must carry text on top
PAL_FILL <- c(
  inferred = "#F3DEDE",
  baseline = "#DCE7F1",
  context  = "#EDEDED",
  truth    = "#FFFFFF"
)

## the second channel. Filled shapes (21-25) so fill = role, outline = ink.
PAL_SHAPE <- c(
  inferred = 24,   # triangle up
  baseline = 21,   # circle
  context  = 23,   # diamond
  truth    = 22    # square
)

## ---------------------------------------------------------------- arm roles
## Which colour each arm takes. `link` and `pseudomultiome` are the inferred pairings in
## current use; the rand_* arms are oracle-label baselines that are given the true cell
## type and involve no inference.
ARM_ROLE <- c(
  truth                    = "truth",
  link                     = "inferred",
  link_1to1                = "inferred",
  pseudomultiome           = "inferred",
  pseudomultiome_k10       = "inferred",
  rand_subcluster_1to1     = "baseline",
  rand_subcluster_matched  = "baseline",
  rand_lineage_matched     = "baseline",
  rand_compartment_matched = "baseline",
  rand_global              = "context"
)

arm_colour <- function(arm) unname(PAL[ARM_ROLE[arm]])
arm_fill   <- function(arm) unname(PAL_FILL[ARM_ROLE[arm]])
arm_shape  <- function(arm) unname(PAL_SHAPE[ARM_ROLE[arm]])

## ---------------------------------------------------------------- arm names
## Two naming styles: precise names for the paper, plain-language names for talks.
ARM_NAME_PAPER <- c(
  truth                    = "truth",
  rand_subcluster_1to1     = "oracle 1-to-1",
  link_1to1                = "link, 1-to-1",
  rand_subcluster_matched  = "oracle, reuse-matched",
  link                     = "link (ArchR)",
  rand_lineage_matched     = "oracle lineage",
  rand_compartment_matched = "oracle compartment",
  rand_global              = "random, global",
  pseudomultiome           = "metacell k = 5",
  pseudomultiome_k10       = "metacell k = 10"
)

## Kept short and parallel in form, so y-axis labels stay narrow and read as a set.
ARM_NAME_TALK <- c(
  truth                    = "the right answer",
  rand_subcluster_1to1     = "random  (each cell once)",
  link_1to1                = "matched  (each cell once)",
  rand_subcluster_matched  = "random  (cells reused)",
  link                     = "standard method",
  rand_lineage_matched     = "random  (coarse groups)",
  rand_compartment_matched = "random  (very coarse)",
  rand_global              = "random  (no groups)",
  pseudomultiome           = "averaging  (5 cells)",
  pseudomultiome_k10       = "averaging  (10 cells)"
)

ARM_ROLE_CARD <- c(
  truth                    = "truth",
  link_1to1                = "one_to_one",
  rand_subcluster_1to1     = "one_to_one",
  link                     = "reuse",
  rand_subcluster_matched  = "reuse",
  rand_lineage_matched     = "reuse",
  rand_compartment_matched = "reuse",
  rand_global              = "reuse",
  pseudomultiome           = "metacell",
  pseudomultiome_k10       = "metacell"
)

PAL_CARD <- c(
  one_to_one = "#2166AC",   # blue
  reuse      = "#B2182B",   # red
  metacell   = "#8C8C8C",   # grey
  truth      = "#111111",   # black
  unavailable= "#8C8C8C"
)

PAL_CARD_FILL <- c(
  one_to_one = "#DCE7F1",
  reuse      = "#F3DEDE",
  metacell   = "#EDEDED",
  truth      = "#FFFFFF",
  unavailable= "#FFFFFF"     # open, so "no value" cannot read as a small value
)

PAL_CARD_SHAPE <- c(
  one_to_one = 21,   # circle
  reuse      = 24,   # triangle up
  metacell   = 23,   # diamond
  truth      = 22,   # square
  unavailable= 4     # cross -- deliberately not a plotting symbol with an area
)

card_role   <- function(arm) unname(ARM_ROLE_CARD[arm])
card_colour <- function(arm) unname(PAL_CARD[ARM_ROLE_CARD[arm]])
card_fill   <- function(arm) unname(PAL_CARD_FILL[ARM_ROLE_CARD[arm]])
card_shape  <- function(arm) unname(PAL_CARD_SHAPE[ARM_ROLE_CARD[arm]])
