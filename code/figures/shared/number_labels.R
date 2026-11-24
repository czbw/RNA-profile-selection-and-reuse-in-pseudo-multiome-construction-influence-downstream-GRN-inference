
## axis/legend labels: minimal decimals, no leading zero, bare "0" for zero
lab_nolead <- function(x) {
  s <- format(x, trim = TRUE, drop0trailing = TRUE, scientific = FALSE)
  s <- sub("^(-?)0\\.", "\\1.", s)
  s[!is.na(x) & abs(x) < 1e-12] <- "0"
  s[is.na(x)] <- NA_character_
  s
}

## the same, as a percent -- ".5%", "1%", "1.5%"
lab_nolead_pct <- function(x) {
  s <- lab_nolead(round(x * 100, 10))
  ifelse(is.na(s), NA_character_, paste0(s, "%"))
}

## a value label at fixed precision: fmt_nolead(0.0344, 4) -> ".0344", signed = TRUE -> "+.0344"
fmt_nolead <- function(x, digits, signed = FALSE) {
  s <- sprintf(paste0(if (signed) "%+." else "%.", digits, "f"), x)
  sub("^([+-]?)0\\.", "\\1.", s)
}

## self-test: the convention is asserted, not assumed
stopifnot(
  identical(lab_nolead(c(-0.05, 0, 0.05, 0.1, 1, 1.5)), c("-.05", "0", ".05", ".1", "1", "1.5")),
  identical(lab_nolead(c(0.95, 0.975, 1, 1.025)), c(".95", ".975", "1", "1.025")),
  identical(lab_nolead_pct(c(0.005, 0.01, 0.015)), c(".5%", "1%", "1.5%")),
  identical(fmt_nolead(c(0.0344, -0.00139), 4, signed = TRUE), c("+.0344", "-.0014")),
  identical(fmt_nolead(1.3788, 4), "1.3788"))
