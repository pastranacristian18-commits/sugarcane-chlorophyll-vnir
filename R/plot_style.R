# ==============================================================================
# R/plot_style.R -- journal-style figures, English labels, exported at PRINT size
#
# Why the manuscript figures looked tiny: they were exported at 7-14 in wide with
# ~10-pt text and then shrunk to ~3-5.5 in on the page, so text printed at ~3-5 pt.
# Here every figure is exported at its final width (fig_width_full / fig_width_half
# in config.R) so the point sizes below are the sizes that are actually printed.
#   standard figures : axis text 9 pt, axis titles 10 pt, legends 9 pt, tags 12 pt
#   Figure 6         : axis text 12 pt, titles 13 pt (>= +60 % vs. the other figures)
# ==============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
PAL <- c(orange = "#E69F00", sky = "#56B4E9", green = "#009E73", yellow = "#F0E442",
         blue = "#0072B2", vermillion = "#D55E00", purple = "#CC79A7", grey = "grey40")
PIG_COL <- c("Chl a" = "#009E73", "Chl b" = "#D55E00")

theme_paper <- function(base = 9) {
  theme_bw(base_size = base, base_family = "sans") +
    theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
          axis.text = element_text(size = base, colour = "black"), axis.title = element_text(size = base + 1),
          legend.text = element_text(size = base), legend.title = element_text(size = base + 1),
          legend.position = "bottom", legend.key.size = unit(0.9, "lines"),
          strip.text = element_text(size = base, face = "bold"), strip.background = element_rect(fill = "grey92", colour = "grey60"),
          plot.tag = element_text(size = base + 3, face = "bold"), plot.margin = margin(4, 6, 4, 4))
}
mm_pt <- function(pt) pt / 2.845              # ggplot geom_text size (mm) for a given point size

save_fig <- function(p, name, width = CFG$fig_width_full, height = 4, cfg = CFG) {
  dir.create(cfg$out_figures, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(cfg$out_figures, paste0(name, ".png")), p, width = width, height = height, dpi = cfg$dpi, bg = "white")
  ggsave(file.path(cfg$out_figures, paste0(name, ".pdf")), p, width = width, height = height, device = cairo_pdf)
  message("[figure] ", name, " (", width, " x ", height, " in)")
}
