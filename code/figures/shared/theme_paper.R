
suppressPackageStartupMessages({
  library(ggplot2); library(grid); library(patchwork)
})
if (!file.exists("figures/shared/palette.R"))
  stop("run from the repo root: Rscript figures/<folder>/make_*.R")
source("figures/shared/palette.R")

MM    <- 1 / 25.4
COL1  <- 85  * MM
COL15 <- 114 * MM
COL2  <- 174 * MM

theme_paper <- function(base_size = 7) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      text            = element_text(colour = PAL[["truth"]]),
      axis.title      = element_text(size = base_size),
      axis.text       = element_text(size = base_size - 0.6, colour = PAL[["truth"]]),
      axis.line       = element_line(colour = PAL[["truth"]], linewidth = 0.3),
      axis.ticks      = element_line(colour = PAL[["truth"]], linewidth = 0.3),
      axis.ticks.length = unit(1.6, "pt"),
      panel.grid      = element_blank(),
      legend.position = "none",
      ## No plot.title is used; the top margin leaves room for the panel tag.
      plot.margin     = margin(9, 3, 2, 2)
    )
}

declare_so_what <- function(x) {
  for (k in names(x)) message(sprintf("  [so what] (%s) %s", k, x[[k]]))
  invisible(x)
}

## patchwork tag theme: lowercase, bold, outside the panel
tag_paper <- function(p, levels = "a") {
  p + plot_annotation(tag_levels = levels) &
    theme(plot.tag = element_text(size = 9, face = "bold", family = "sans",
                                  colour = PAL[["truth"]], hjust = 0, vjust = 1),
          plot.tag.position = c(0, 1))
}

save_paper <- function(plot, path_noext, w, h, dpi = 300) {
  grDevices::cairo_pdf(paste0(path_noext, ".pdf"), width = w, height = h, onefile = FALSE)
  print(plot); grDevices::dev.off()
  ragg::agg_png(paste0(path_noext, ".png"), width = w, height = h, units = "in",
                res = dpi, background = "white")
  print(plot); grDevices::dev.off()
  message(sprintf("  wrote %s.{pdf,png}  (%.0f x %.0f mm, %d dpi)",
                  basename(path_noext), w / MM, h / MM, dpi))
}

save_paper_grob <- function(grob, path_noext, w, h, dpi = 300) {
  grDevices::cairo_pdf(paste0(path_noext, ".pdf"), width = w, height = h, onefile = FALSE)
  grid.newpage(); grid.draw(grob); grDevices::dev.off()
  ragg::agg_png(paste0(path_noext, ".png"), width = w, height = h, units = "in",
                res = dpi, background = "white")
  grid.newpage(); grid.draw(grob); grDevices::dev.off()
  message(sprintf("  wrote %s.{pdf,png}  (%.0f x %.0f mm, %d dpi)",
                  basename(path_noext), w / MM, h / MM, dpi))
}
