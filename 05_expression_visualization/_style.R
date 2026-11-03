# =============================================================================
# Shared style and layout helpers for the plotting scripts.
#
# Sourced by every plotting script so the plots cannot drift apart in font,
# palette, tile geometry or panel-tag placement.
# =============================================================================

# The container's default locale is "C", which mangles the Greek delta and the
# typographic minus at PARSE time - long before the graphics device sees them.
Sys.setlocale("LC_ALL", "C.utf8")

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
  library(patchwork); library(viridisLite)
  library(grid); library(gtable); library(systemfonts)
})

# Arial is not installed here. Liberation Sans is the metric-compatible clone
# (identical advance widths), so registering it under the name "Arial" keeps
# every measurement identical to a machine with real Arial. On a Mac the
# file.exists() guard fails, the registration is skipped, and FONT picks up the
# real system Arial.
LIB <- "/usr/share/fonts/truetype/liberation"
if (file.exists(file.path(LIB, "LiberationSans-Regular.ttf"))) {
  register_font("Arial",
                plain      = file.path(LIB, "LiberationSans-Regular.ttf"),
                bold       = file.path(LIB, "LiberationSans-Bold.ttf"),
                italic     = file.path(LIB, "LiberationSans-Italic.ttf"),
                bolditalic = file.path(LIB, "LiberationSans-BoldItalic.ttf"))
}

# Text is measured on the CURRENT device. The default postscript device does
# not know fonts registered with systemfonts, so every string would be
# mis-measured and every layout computed from those measurements would be wrong.
ragg::agg_capture(width = 900, height = 1200, units = "px", res = 300)

GREEN  <- "#5AAE61"   # control  (house palette)
PURPLE <- "#6B3FA0"   # glucose
BLUE   <- "#1F5FA9"   # Bacteria subtitle
BASE   <- 7           # hard floor: nothing on either figure is smaller
FONT   <- "Arial"
VIR    <- viridis(256)

GROUPS <- c("C24", "G24", "C48", "G48", "C72", "G72")
CONTRAST_LEVELS <- c("G24 vs C24", "G48 vs C48", "G72 vs C72",
                     "C48 vs C24", "C72 vs C24", "C72 vs C48",
                     "G48 vs G24", "G72 vs G24", "G72 vs G48")
CONTRAST_KEYS <- c("glucose_G24_vs_C24", "glucose_G48_vs_C48", "glucose_G72_vs_C72",
                   "time_C48_vs_C24", "time_C72_vs_C24", "time_C72_vs_C48",
                   "time_G48_vs_G24", "time_G72_vs_G24", "time_G72_vs_G48")

theme_fig <- function(base = BASE) {
  theme_classic(base_size = base, base_family = FONT) +
    theme(
      plot.title    = element_text(size = base, face = "bold", hjust = 0,
                                   margin = margin(b = 2)),
      plot.subtitle = element_text(size = base, hjust = 0.5, face = "bold"),
      axis.title    = element_text(size = base),
      axis.text     = element_text(size = base, colour = "black"),
      axis.line     = element_line(linewidth = 0.3),
      axis.ticks    = element_line(linewidth = 0.3),
      legend.title  = element_text(size = base, face = "bold"),
      legend.text   = element_text(size = base),
      legend.key.size = unit(0.28, "cm"),
      legend.margin = margin(1, 1, 1, 1),
      plot.margin   = margin(2, 3, 2, 3)
    )
}
theme_heat <- function(...) {
  theme_fig() + theme(axis.line = element_blank(), axis.ticks = element_blank(),
                      panel.grid = element_blank(), ...)
}

# --- layout machinery --------------------------------------------------------
# Tile size is fixed in millimetres by rewriting the panel cell of the gtable to
# an absolute unit (the trick egg::set_panel_size uses). aspect.ratio would only
# fix the SHAPE, not the size: a tile still stretches when the layout hands the
# panel more room.
set_panel_mm <- function(p, w_mm, h_mm) {
  g <- ggplotGrob(p)
  is_panel <- grepl("^panel", g$layout$name)
  g$widths[ unique(g$layout$l[is_panel])] <- unit(w_mm, "mm")
  g$heights[unique(g$layout$t[is_panel])] <- unit(h_mm, "mm")
  g
}
gh <- function(g) sum(convertHeight(g$heights, "mm", valueOnly = TRUE))
gw <- function(g) sum(convertWidth( g$widths,  "mm", valueOnly = TRUE))
# A guide box wider than the plot OVERFLOWS its cell instead of widening the
# gtable, so gw() alone under-reports and the layout crops the legend. Measure
# the guide grob itself and take whichever is wider.
gw_full <- function(g) {
  w  <- gw(g)
  gb <- which(g$layout$name == "guide-box")
  if (length(gb)) {
    lw <- convertWidth(grobWidth(g$grobs[[gb[1]]]), "mm", valueOnly = TRUE)
    w  <- max(w, lw + 3)
  }
  w
}

# Panel tags are drawn into the grob. patchwork's tag_levels puts the tag at the
# very edge of the layout cell, where a top-row tag falls off the canvas.
add_tag <- function(g, lab) {
  g2 <- gtable_add_rows(g, unit(4.2, "mm"), pos = 0)
  gtable_add_grob(g2,
    textGrob(lab, x = unit(0.4, "mm"), y = unit(0.2, "mm"), hjust = 0, vjust = 0,
             gp = gpar(fontsize = BASE + 3, fontface = "bold", fontfamily = FONT)),
    t = 1, l = 1, r = ncol(g2), name = "tag", clip = "off")
}
# wrap_elements() centres its grob; pad on the right instead so a narrow panel
# stays flush with the panels above and below it
pad_right <- function(g, target_mm) {
  d <- target_mm - gw(g)
  if (d <= 0) return(g)
  gtable_add_cols(g, unit(d, "mm"), pos = -1)
}
# Text colour that stays readable on a viridis tile.
on_viridis <- function(v, lo, hi) {
  t <- pmin(pmax((v - lo) / (hi - lo), 0), 1)
  rgbv <- col2rgb(VIR[pmax(1, round(t * 255) + 1)]) / 255
  lum  <- 0.2126 * rgbv[1, ] + 0.7152 * rgbv[2, ] + 0.0722 * rgbv[3, ]
  ifelse(lum > 0.55, "black", "white")
}

save_fig <- function(fig, stem, w_mm, h_mm) {
  if (w_mm > 210) message("!! ", basename(stem), " is ", round(w_mm, 1),
                          " mm wide - over the 210 mm limit")
  ggsave(paste0(stem, ".pdf"), fig, width = w_mm, height = h_mm,
         units = "mm", device = cairo_pdf)
  ggsave(paste0(stem, ".svg"), fig, width = w_mm, height = h_mm,
         units = "mm", device = svglite::svglite)
  ggsave(paste0(stem, ".png"), fig, width = w_mm, height = h_mm,
         units = "mm", dpi = 600, device = ragg::agg_png)
  cat(sprintf("%-34s %6.1f x %6.1f mm\n", basename(stem), w_mm, h_mm))
}
