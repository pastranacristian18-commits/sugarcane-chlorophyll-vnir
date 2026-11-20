# ==============================================================================
# 14_toc_graphic.R -- table-of-contents (TOC) graphic for the journal submission
# Size 3.25 x 1.75 in, 600 dpi. Values are read from the tables written by script 03:
#   RPD (Chl a) under random 5-fold, leave-one-plot-out and leave-one-campaign-out validation,
#   and the RMSE ratio (campaign held out / random) for Chl a and Chl b.
# Output: results/figures/toc_graphic.png and .pdf
# ==============================================================================
source("config.R")
suppressPackageStartupMessages(library(ggplot2))
tab  <- read.csv(file.path(CFG$out_tables, "table03_validation_schemes.csv"), check.names = FALSE)
rat  <- read.csv(file.path(CFG$out_tables, "table03b_scheme_ratios_vs_random.csv"), check.names = FALSE)
rpd  <- function(s) tab$RPD[tab$pigment == "Chl a" & tab$scheme == s & tab$predictor == "PLSR"]
rr   <- function(p) rat$RMSE_ratio_vs_random[rat$pigment == p & rat$scheme == "Leave-one-campaign-out"]
bars <- data.frame(
  scheme = factor(c("Random\n5-fold", "Plots\nheld out", "Campaigns\nheld out"),
                  levels = c("Random\n5-fold", "Plots\nheld out", "Campaigns\nheld out")),
  rpd = c(rpd("Random 5-fold (leaves)"), rpd("Leave-one-plot-out"), rpd("Leave-one-campaign-out")),
  col = c("#2a78d6", "#2a78d6", "#eb6834"))
up_a <- round(100 * (rr("Chl a") - 1)); up_b <- round(100 * (rr("Chl b") - 1))
INK <- "#0b0b0b"; INK2 <- "#52514e"; ORANGE <- "#eb6834"
pt <- function(x) x / 2.845   # point size -> ggplot text size (mm)

bar <- ggplot(bars, aes(scheme, rpd, fill = col)) +
  geom_col(width = 0.58) + scale_fill_identity() +
  geom_text(aes(label = sprintf("%.2f", rpd)), vjust = -0.45, size = pt(6.4), fontface = "bold", colour = INK) +
  scale_y_continuous(limits = c(0, 1.75), breaks = c(0, 0.5, 1, 1.5), expand = c(0, 0)) +
  labs(x = NULL, y = "RPD (Chl a)") +
  theme_minimal(base_size = 5.6, base_family = "sans") +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = "#d9d8d4", linewidth = 0.2),
        axis.text.x = element_text(colour = INK, size = 5.6, lineheight = 0.95),
        axis.text.y = element_text(colour = INK2, size = 5.6),
        axis.title.y = element_text(colour = INK2, size = 5.8),
        axis.line.x = element_line(colour = INK2, linewidth = 0.3),
        plot.margin = margin(0, 0, 0, 0))

txt <- ggplot() + xlim(0, 1) + ylim(0, 1) + theme_void() +
  annotate("text", 0.5, 0.80, label = paste0("+", up_a, "%"), size = pt(19), fontface = "bold", colour = ORANGE) +
  annotate("text", 0.5, 0.52, label = sprintf("RMSE for Chl a\nwith campaigns held out\n(+%d%% for Chl b)", up_b),
           size = pt(5.6), colour = INK, lineheight = 1.05) +
  annotate("text", 0.5, 0.25, label = "vs random validation", size = pt(5.4), colour = INK2)

p <- ggplot() + xlim(0, 1) + ylim(0, 1) + theme_void() +
  annotate("text", 0.03, 0.93, label = "Random validation overstates transferability", hjust = 0, size = pt(7.6), fontface = "bold", colour = INK) +
  annotate("text", 0.03, 0.835, label = "Leaf chlorophyll from VIS-NIR reflectance, sugarcane", hjust = 0, size = pt(5.6), colour = INK2) +
  annotate("text", 0.03, 0.045, label = "477 leaves, 20 plots, 6 campaigns, PLSR", hjust = 0, size = pt(5), colour = INK2) +
  annotation_custom(ggplotGrob(bar), xmin = 0.03, xmax = 0.62, ymin = 0.09, ymax = 0.76) +
  annotation_custom(ggplotGrob(txt), xmin = 0.64, xmax = 0.97, ymin = 0.10, ymax = 0.78)
dir.create(CFG$out_figures, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(CFG$out_figures, "toc_graphic.png"), p, width = 3.25, height = 1.75, dpi = 600, bg = "white",
       device = png, type = "cairo")
ggsave(file.path(CFG$out_figures, "toc_graphic.pdf"), p, width = 3.25, height = 1.75, device = cairo_pdf)
message("TOC graphic written to ", CFG$out_figures)
