# ==============================================================================
# 09_figures.R -- all manuscript figures (English, print size, colour-blind palette)
# Figure numbers follow the manuscript (Figs. 2-9); new diagnostic figures are 10-12.
# Run after scripts 01-08. Output: results/figures/*.png (600 dpi) and *.pdf (vector)
# ==============================================================================
source("config.R"); source("R/functions.R"); source("R/plot_style.R")
suppressPackageStartupMessages(library(readxl))
D <- load_soca(); meta <- D$meta; X <- D$X; wl <- D$wl
W <- CFG$fig_width_full; Wh <- CFG$fig_width_half
dac_lab <- function(d) factor(paste(d, "DAC"), levels = paste(sort(unique(d)), "DAC"))
RE <- annotate("rect", xmin = 700, xmax = 770, ymin = -Inf, ymax = Inf, fill = "steelblue", alpha = 0.10)

# ---- Fig. 2  Pearson correlation profiles ------------------------------------
prof <- readRDS(file.path(CFG$out_rds, "pearson_profiles.rds")); prof$panel <- dac_lab(prof$DAC)
p2 <- ggplot(prof, aes(wavelength, r_smooth, colour = pigment)) + RE +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.3) +
  geom_hline(yintercept = c(-0.44, 0.44), linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 0.6) + facet_wrap(~panel, ncol = 3) +
  scale_colour_manual(values = PIG_COL, name = NULL) + coord_cartesian(ylim = c(-0.6, 0.6)) +
  labs(x = "Wavelength (nm)", y = "Pearson r (LOESS-smoothed)") + theme_paper(9)
save_fig(p2, "fig02_pearson_profiles", W, 4.6)

# ---- Fig. 3  raw vs processed spectra (mean +/- SD, individual spectra in grey) -
raw <- as.data.frame(read_excel(CFG$data_raw)); rk <- paste(raw$Dia, raw$Unidad)
Xr <- as.matrix(raw[, -(1:2)])[match(meta$id, rk), ]
spec_df <- function(M, lab) { mu <- colMeans(M); sd <- apply(M, 2, sd)
  list(band = data.frame(w = wl, lo = mu - sd, hi = mu + sd, mu = mu, panel = lab),
       lines = data.frame(w = rep(wl, each = nrow(M)), id = rep(seq_len(nrow(M)), ncol(M)), v = as.vector(M), panel = lab)) }
a <- spec_df(Xr, "(a) Raw"); b <- spec_df(X, "(b) SG + MSC")
bands <- rbind(a$band, b$band); lines <- rbind(a$lines[a$lines$id %% 2 == 0, ], b$lines[b$lines$id %% 2 == 0, ])   # every 2nd spectrum (file size)
p3 <- ggplot() + geom_line(data = lines, aes(w, v, group = interaction(id, panel)), colour = "grey70", linewidth = 0.15, alpha = 0.5) +
  geom_ribbon(data = bands, aes(w, ymin = lo, ymax = hi), fill = PAL[["blue"]], alpha = 0.25) +
  geom_line(data = bands, aes(w, mu), colour = PAL[["blue"]], linewidth = 0.7) +
  facet_wrap(~panel, scales = "free_y") + labs(x = "Wavelength (nm)", y = "Reflectance (%)") + theme_paper(9)
save_fig(p3, "fig03_spectra_raw_vs_processed", W, 3.2)

# ---- Figs. 4-5  observed vs predicted, random 80/20 (numbering as in the manuscript) ---
fits <- readRDS(file.path(CFG$out_rds, "random_8020_fits.rds"))
obs_pred <- function(v, lab) { f <- fits[[v]]; mv <- metrics(f$obs_val, f$pred_val)
  d <- rbind(data.frame(obs = f$obs_cal, pred = f$pred_cal, set = "Calibration"), data.frame(obs = f$obs_val, pred = f$pred_val, set = "Validation"))
  txt <- sprintf("r\u00b2 = %.2f\nR\u00b2 = %.2f\nRMSE = %.3f\nRPD = %.2f\nBias = %.3f", mv$r2, mv$R2, mv$RMSE, mv$RPD, mv$Bias)
  lim <- range(c(d$obs, d$pred))
  ggplot(d, aes(obs, pred)) + geom_abline(slope = 1, intercept = 0, colour = "grey30", linetype = "dashed", linewidth = 0.4) +
    geom_smooth(data = d[d$set == "Validation", ], method = "lm", se = FALSE, colour = PAL[["vermillion"]], linewidth = 0.7, formula = y ~ x) +
    geom_point(aes(fill = set, shape = set), size = 1.9, alpha = 0.85, colour = "black", stroke = 0.25) +
    scale_shape_manual(values = c(Calibration = 21, Validation = 24), name = NULL) + scale_fill_manual(values = c(Calibration = PAL[["sky"]], Validation = PAL[["orange"]]), name = NULL) +
    annotate("text", x = lim[1], y = lim[2], label = txt, hjust = 0, vjust = 1, size = mm_pt(8.5), lineheight = 0.95) +
    coord_equal(xlim = lim, ylim = lim) + labs(title = NULL, x = bquote("Observed " * .(lab) ~ "(mg" ~ g^-1 ~ "FW)"), y = bquote("Predicted " * .(lab) ~ "(mg" ~ g^-1 ~ "FW)")) + theme_paper(9) }
save_fig(obs_pred("Chla", "Chl a"), "fig04_obs_vs_pred_chla", Wh + 0.4, 4.0)
save_fig(obs_pred("Chlb", "Chl b"), "fig05_obs_vs_pred_chlb", Wh + 0.4, 4.0)

# ---- Fig. 6  spectral importance (fonts >= +60 %) ----------------------------
imp <- readRDS(file.path(CFG$out_rds, "importance_profiles.rds"))
il <- rbind(data.frame(imp[, c("pigment", "wavelength")], value = imp$RF_importance, measure = "PCA-RF\nimportance"),
            data.frame(imp[, c("pigment", "wavelength")], value = imp$VIP, measure = "PLSR\nVIP"))
il$panel <- factor(ifelse(il$pigment == "Chl a", "(a) Chl a", "(b) Chl b")); il$measure <- factor(il$measure, levels = unique(il$measure))
p6 <- ggplot(il, aes(wavelength, value)) + RE +
  geom_hline(data = data.frame(measure = factor("PLSR\nVIP", levels = levels(il$measure)), y = 1), aes(yintercept = y), colour = PAL[["vermillion"]], linetype = "dashed", linewidth = 0.6) +
  geom_line(aes(colour = measure), linewidth = 0.8) + scale_colour_manual(values = c(PAL[["blue"]], PAL[["green"]]), guide = "none") +
  facet_grid(measure ~ panel, scales = "free_y") + labs(x = "Wavelength (nm)", y = "Importance score (PCA-RF normalised to 1)") + theme_paper(12) +
  theme(strip.text = element_text(size = 11, face = "bold"), axis.text = element_text(size = 12), axis.title = element_text(size = 13))
save_fig(p6, "fig06_spectral_importance", W, 6.0)

# ---- Fig. 7  like-for-like comparison of validation schemes -------------------
T3 <- read.csv(file.path(CFG$out_tables, "table03_validation_schemes.csv"), check.names = FALSE)
lv <- c("Random 5-fold (leaves)", "Leave-one-plot-out", "Leave-one-block-out", "Leave-one-campaign-out", "New campaign + new plots")
short <- c("Random\n5-fold", "Leave-one-\nplot-out", "Leave-one-\nblock-out", "Leave-one-\ncampaign-out", "Campaign +\nplots held out")
T3$scheme <- factor(T3$scheme, levels = lv, labels = short)
pl <- T3[T3$predictor == "PLSR", ]; pl$kind <- "PLSR on spectra (95% CI, plot bootstrap)"
bl <- T3[T3$predictor != "PLSR" & ((T3$predictor == "Campaign mean" & T3$scheme %in% short[1:3]) | (T3$predictor == "Training mean" & T3$scheme %in% short[4:5])), ]
bl$kind <- ifelse(bl$predictor == "Campaign mean", "Campaign-mean predictor (no spectra)", "Training-mean predictor (no spectra)")
p7 <- ggplot(pl, aes(scheme, RPD)) + geom_hline(yintercept = 1, linetype = "dashed", colour = PAL[["vermillion"]], linewidth = 0.5) +
  geom_pointrange(aes(ymin = RPD_lo, ymax = RPD_hi), colour = PAL[["blue"]], size = 0.5, linewidth = 0.6) +
  geom_point(data = bl, aes(shape = kind), colour = "grey20", size = 2.8, fill = "white", stroke = 0.9) +
  scale_shape_manual(values = c("Campaign-mean predictor (no spectra)" = 24, "Training-mean predictor (no spectra)" = 22), name = NULL) +
  facet_wrap(~pigment) + labs(x = NULL, y = "RPD (pooled out-of-fold predictions)", caption = "Blue: PLSR (95% CI from resampling the 20 plots). Open symbols: predictors that use no spectra.\nDashed line: RPD = 1 (no better than predicting the mean).") +
  theme_paper(9) + theme(axis.text.x = element_text(size = 7, lineheight = 0.9), plot.caption = element_text(size = 7.5, hjust = 0))
save_fig(p7, "fig07_validation_schemes", W, 4.2)

# ---- Fig. 7-S  Monte-Carlo random 80/20 vs LOCO (original style, supplementary) --------
mcf <- file.path(CFG$out_rds, "monte_carlo_8020.rds")
if (file.exists(mcf)) { mc <- readRDS(mcf)
  d <- do.call(rbind, lapply(PIGMENTS, function(v) data.frame(pigment = unname(PIG_LABEL[v]), RPD = mc[[v]][, "RPD"])))
  T4 <- read.csv(file.path(CFG$out_tables, "table04_loco_per_campaign.csv"), check.names = FALSE); T4$RPD_fold <- T4$RPD
  pS <- ggplot(d, aes(pigment, RPD)) + geom_violin(fill = PAL[["sky"]], colour = NA, alpha = 0.6, width = 0.8) + geom_boxplot(width = 0.12, outlier.shape = NA, fill = PAL[["sky"]]) +
    geom_jitter(data = T4, aes(pigment, RPD_fold, colour = "LOCO, per campaign"), width = 0.1, shape = 21, size = 2, stroke = 0.8, fill = "white") +
    geom_hline(yintercept = 1, linetype = "dashed", colour = PAL[["vermillion"]]) + scale_colour_manual(values = c("LOCO, per campaign" = PAL[["vermillion"]]), name = NULL) +
    labs(x = NULL, y = "RPD", caption = "Violin: repeated random 80/20 splits. Per-campaign RPD uses the SD of a single campaign, so it is not comparable with the pooled RPD (Fig. 7).") +
    theme_paper(9) + theme(plot.caption = element_text(size = 7.5, hjust = 0))
  save_fig(pS, "fig07s_monte_carlo_vs_loco_per_fold", Wh + 1.2, 3.6) }

# ---- Fig. 8  stage-specific models under three validation designs -------------
T9 <- read.csv(file.path(CFG$out_tables, "table09_stage_model_diagnosis.csv"), check.names = FALSE)
d8 <- rbind(data.frame(T9[, c("pigment", "DDC")], r2 = T9$manuscript_KS_r2_median, design = "Manuscript procedure (median of 20 CV seeds)"),
            data.frame(T9[, c("pigment", "DDC")], r2 = T9$naive_leaf_r2, design = "Random leaf-level 5-fold"),
            data.frame(T9[, c("pigment", "DDC")], r2 = T9$plotCV_groupedLV_r2, design = "Whole plots held out (plot-grouped CV)"))
d8$design <- factor(d8$design, levels = unique(d8$design))
p8 <- ggplot(d8, aes(factor(DDC), r2, colour = design, shape = design, group = design)) + geom_line(linewidth = 0.5, alpha = 0.7) + geom_point(size = 2.4) +
  facet_wrap(~pigment) + scale_colour_manual(values = c(PAL[["vermillion"]], PAL[["orange"]], PAL[["blue"]]), name = NULL) + scale_shape_manual(values = c(17, 15, 16), name = NULL) +
  guides(colour = guide_legend(nrow = 3, byrow = TRUE), shape = guide_legend(nrow = 3, byrow = TRUE)) + labs(x = "Days after cutting (DAC)", y = expression("Validation " * r^2 * " (within campaign)")) + theme_paper(9)
save_fig(p8, "fig08_stage_models_by_validation_design", W, 4.0)

# ---- Fig. 9  chlorophyll trajectory by N dose --------------------------------
tl <- c(T1 = "0 kg N/ha", T2 = "80 kg N/ha", T3 = "160 kg N/ha", T4 = "240 kg N/ha")
d9 <- do.call(rbind, lapply(PIGMENTS, function(v) { a <- aggregate(meta[[v]], list(DAC = meta$DDC, trt = meta$trt), function(z) c(m = mean(z), s = sd(z)))
  data.frame(pigment = unname(PIG_LABEL[v]), DAC = a$DAC, trt = factor(tl[a$trt], levels = tl), m = a$x[, "m"], s = a$x[, "s"]) }))
p9 <- ggplot(d9, aes(DAC, m, colour = trt, group = trt)) +
  annotate("rect", xmin = 45, xmax = 135, ymin = -Inf, ymax = Inf, fill = PAL[["sky"]], alpha = 0.12) + annotate("rect", xmin = 135, xmax = 225, ymin = -Inf, ymax = Inf, fill = PAL[["green"]], alpha = 0.12) +
  geom_vline(xintercept = 135, linetype = "dotted", colour = "grey40") + geom_errorbar(aes(ymin = m - s, ymax = m + s), width = 3.5, alpha = 0.5, linewidth = 0.35) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.8) + facet_wrap(~pigment, scales = "free_y") +
  scale_colour_manual(values = c(PAL[["orange"]], PAL[["sky"]], PAL[["green"]], PAL[["blue"]]), name = "N dose") + scale_x_continuous(breaks = seq(60, 210, 30)) +
  labs(x = "Days after cutting (DAC)", y = expression("Chlorophyll (mg " * g^-1 * " FW), mean " %+-% " SD"), caption = "Blue: tillering (60-120 DAC); green: elongation (150-210 DAC). Treatments are not balanced across blocks 4-5 (see audit).") +
  theme_paper(9) + theme(plot.caption = element_text(size = 7.5, hjust = 0))
save_fig(p9, "fig09_chlorophyll_trajectory_by_N_dose", W, 3.6)

# ---- Fig. 10  EPO: correlation vs. determination vs. error ---------------------
T5 <- read.csv(file.path(CFG$out_tables, "table05_epo_loco.csv"), check.names = FALSE)
keep <- c("No EPO", "EPO leaky (all campaigns), k=5", "EPO training-only, k=4", "EPO nested-selected k")
lab  <- c("No EPO", "EPO as in original\n(all campaigns, k=5)", "EPO training-only\n(k=4)", "EPO nested\n(k chosen in training)")
e <- T5[T5$variant %in% keep, ]; e$variant <- factor(lab[match(e$variant, keep)], levels = lab)
e <- rbind(data.frame(e[, c("pigment", "variant")], value = e$r2, metric = "r\u00b2 (squared correlation)"),
           data.frame(e[, c("pigment", "variant")], value = e$R2, metric = "R\u00b2 (determination)"),
           data.frame(e[, c("pigment", "variant")], value = e$RMSE, metric = "RMSE (mg/g)"))
e$metric <- factor(e$metric, levels = unique(e$metric))
p10 <- ggplot(e, aes(variant, value, fill = pigment)) + geom_col(position = position_dodge(0.75), width = 0.7, colour = "black", linewidth = 0.25) +
  geom_hline(data = data.frame(metric = factor("R\u00b2 (determination)", levels = levels(e$metric)), y = 0), aes(yintercept = y), linewidth = 0.3) +
  facet_wrap(~metric, scales = "free_y", ncol = 3) + scale_fill_manual(values = PIG_COL, name = NULL) + labs(x = NULL, y = NULL) +
  theme_paper(9) + theme(axis.text.x = element_text(size = 7, angle = 35, hjust = 1))
save_fig(p10, "fig10_epo_correlation_vs_error", W, 3.6)

# ---- Fig. 11  why r2, R2 and RPD disagree under LOCO ---------------------------
vs <- readRDS(file.path(CFG$out_rds, "validation_schemes.rds"))$res
d11 <- do.call(rbind, lapply(PIGMENTS, function(v) data.frame(pigment = unname(PIG_LABEL[v]), obs = meta[[v]], pred = vs[[v]]$preds[["PLSR|loco"]], DAC = dac_lab(meta$DDC))))
ann <- do.call(rbind, lapply(PIGMENTS, function(v) { m <- metrics(meta[[v]], vs[[v]]$preds[["PLSR|loco"]]); data.frame(pigment = unname(PIG_LABEL[v]), obs = -Inf, pred = Inf,
  txt = sprintf("r\u00b2 = %.2f\nR\u00b2 = %.2f\nRPD = %.2f", m$r2, m$R2, m$RPD)) }))
p11 <- ggplot(d11, aes(obs, pred)) + geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey30", linewidth = 0.4) +
  geom_point(aes(colour = DAC), size = 1.4, alpha = 0.8) + geom_text(data = ann, aes(label = txt), hjust = -0.1, vjust = 1.2, size = mm_pt(8.5), lineheight = 0.95) +
  facet_wrap(~pigment, scales = "free") + scale_colour_manual(values = unname(c(PAL[c("vermillion", "orange", "yellow", "green", "sky", "blue")])), name = NULL) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 2.5))) + labs(x = expression("Observed (mg " * g^-1 * " FW)"), y = expression("Predicted, LOCO (mg " * g^-1 * " FW)")) + theme_paper(9)
save_fig(p11, "fig11_loco_predictions_by_campaign", W, 3.7)

# ---- Fig. 12  within-campaign skill: campaign-centred models -------------------
T6 <- read.csv(file.path(CFG$out_tables, "table06_within_campaign_skill.csv"), check.names = FALSE)
T6 <- T6[!grepl("^  ", T6$scheme), ]; T6$scheme <- ifelse(grepl("plot", T6$scheme), "Held-out plots\n(campaign known)", "Held-out campaign\n(offset removed)")
p12 <- ggplot(T6, aes(scheme, r2_dev, fill = pigment)) + geom_col(position = position_dodge(0.7), width = 0.6, colour = "black", linewidth = 0.25) +
  scale_fill_manual(values = PIG_COL, name = NULL) + scale_y_continuous(limits = c(0, 0.2)) +
  labs(x = NULL, y = expression("Within-campaign " * r^2)) + theme_paper(9)
save_fig(p12, "fig12_within_campaign_skill", Wh + 1.2, 3.3)
message("[done] figures written to ", CFG$out_figures)
