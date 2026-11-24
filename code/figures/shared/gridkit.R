
suppressPackageStartupMessages(library(grid))

## width/height a string will occupy, in npc, at a given size
str_w <- function(txt, fontsize, fontface = "plain") {
  max(vapply(strsplit(txt, "\n")[[1]], function(s)
    convertWidth(grobWidth(textGrob(s, gp = gpar(fontsize = fontsize, fontface = fontface))),
                 "npc", valueOnly = TRUE), numeric(1)))
}
str_h <- function(txt, fontsize, lineheight = 1.2) {
  n <- length(strsplit(txt, "\n")[[1]])
  n * convertHeight(unit(fontsize * lineheight, "points"), "npc", valueOnly = TRUE)
}

## A labelled box. `w`/`h` are minimum sizes: the box grows to fit its contents.
## Returns list(grob, x, y, w, h, ...) so callers can anchor arrows to the real edges.
gbox <- function(x, y, label, sub = NULL, w = 0, h = 0,
                 fill = "white", border = "#111111", lwd = 1.2,
                 fontsize = 13, subsize = 10.5, fontface = "bold",
                 col = "#111111", subcol = "#8C8C8C", pad_x = 0.018, pad_y = 0.030,
                 r = 0.012) {
  need_w <- max(str_w(label, fontsize, fontface),
                if (is.null(sub)) 0 else str_w(sub, subsize)) + 2 * pad_x
  need_h <- str_h(label, fontsize) +
            (if (is.null(sub)) 0 else str_h(sub, subsize) + 0.010) + 2 * pad_y
  W <- max(w, need_w); H <- max(h, need_h)

  ty <- if (is.null(sub)) y else y + H * 0.13
  g <- gList(
    roundrectGrob(x = x, y = y, width = W, height = H, r = unit(r, "npc"),
                  gp = gpar(fill = fill, col = border, lwd = lwd)),
    textGrob(label, x = x, y = ty,
             gp = gpar(fontsize = fontsize, fontface = fontface, col = col))
  )
  if (!is.null(sub))
    g <- gList(g, textGrob(sub, x = x, y = y - H * 0.22,
                           gp = gpar(fontsize = subsize, col = subcol, lineheight = 1.15)))
  list(grob = g, x = x, y = y, w = W, h = H,
       left = x - W / 2, right = x + W / 2, top = y + H / 2, bottom = y - H / 2)
}

## straight arrow between two points
garrow <- function(x1, y1, x2, y2, col = "#111111", lwd = 1.2, len = 2.6) {
  linesGrob(x = c(x1, x2), y = c(y1, y2),
            arrow = arrow(length = unit(len, "mm"), type = "closed"),
            gp = gpar(col = col, lwd = lwd, fill = col))
}

## curved arrow, for many-to-one fans where straight lines would overlap
gcurve <- function(x1, y1, x2, y2, col = "#111111", lwd = 1.2, curvature = 0.18,
                   len = 2.6) {
  curveGrob(x1, y1, x2, y2, default.units = "npc", curvature = curvature,
            square = FALSE, ncp = 12,
            arrow = arrow(length = unit(len, "mm"), type = "closed"),
            gp = gpar(col = col, lwd = lwd, fill = col))
}

## Stop if two boxes overlap. Called by the schematic scripts on every pair that
## shares a band, so a layout error halts the build.
assert_clear <- function(a, b, gap = 0.004, what = "boxes") {
  ox <- min(a$right, b$right) - max(a$left, b$left)
  oy <- min(a$top, b$top)   - max(a$bottom, b$bottom)
  if (ox > -gap && oy > -gap)
    stop(sprintf("%s overlap: x by %.4f, y by %.4f npc", what, ox, oy))
  invisible(TRUE)
}
