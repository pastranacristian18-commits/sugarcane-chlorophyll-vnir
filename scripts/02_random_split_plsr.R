# ==============================================================================
# 02_random_split_plsr.R -- manuscript Table 2 (and the model behind Fig. 2, Fig. 4, Table S1):
# PLSR on ONE random 80/20 split stratified by campaign (literature-style comparison).
#   * split: set.seed(42), 80 % of the leaves of every campaign for calibration (381 / 96 leaves)
#   * preprocessing: Savitzky-Golay spectra (MATLAB export, no MSC) + MSC whose reference
#     is the mean of the CALIBRATION spectra only; the same reference is applied to the validation set
#   * latent variables (1-15): minimum RMSE of an inner 5-fold CV grouped by plot (seed = CFG$seed)
# ==============================================================================
source("config.R"); source("R/functions.R"); suppressPackageStartupMessages(library(readxl))
D <- load_soca(); meta <- D$meta; X <- D$X
sg <- as.data.frame(read_excel(CFG$data_sgonly)); key <- paste(sg$Dia, sg$Unidad)
Xs <- as.matrix(sg[, -(1:2)]); rownames(Xs) <- key; Xs <- Xs[meta$id, ]            # align by key, never by position
stopifnot(all(rownames(Xs) == meta$id), ncol(Xs) == ncol(X))

set.seed(CFG$seed)
idx_cal <- c()
for (cm in unique(as.character(meta$campaign))) {
  ic <- which(meta$campaign == cm); idx_cal <- c(idx_cal, sample(ic, floor(0.8 * length(ic)))) }
idx_val <- setdiff(seq_len(nrow(meta)), idx_cal)

pp <- pp_msc(Xs[idx_cal, ], Xs[idx_val, ], NULL, NULL)          # MSC reference from calibration spectra only
tab <- list(); store <- list()
for (v in PIGMENTS) {
  ycal <- meta[[v]][idx_cal]; yval <- meta[[v]][idx_val]
  f  <- fit_pls(pp$train, ycal, groups = meta$plot[idx_cal], seed = CFG$seed)
  lv <- f$lv; pc <- predict_pls(f, pp$train); pv <- predict_pls(f, pp$test)
  tab[[v]] <- rbind(cbind(pigment = PIG_LABEL[v], set = "Calibration", LV = lv, metrics(ycal, pc)),
                    cbind(pigment = PIG_LABEL[v], set = "Validation 80/20", LV = lv, metrics(yval, pv)))
  store[[v]] <- list(model = f$model, lv = lv, idx_cal = idx_cal, idx_val = idx_val,
                     obs_cal = ycal, obs_val = yval, pred_cal = pc, pred_val = pv)
}
T2 <- do.call(rbind, tab); rownames(T2) <- NULL
T2[, c("r2", "R2", "RMSE", "RPD", "RPIQ", "Bias")] <- round(T2[, c("r2", "R2", "RMSE", "RPD", "RPIQ", "Bias")], 3)
write_table(T2, "table02_plsr_random_8020.csv"); saveRDS(store, file.path(CFG$out_rds, "random_8020_fits.rds"))
cat("\nr2 = squared correlation; R2 = 1 - SSE/SST (see REPRODUCIBILITY.md). Expected: 11 LV for both pigments.\n")
print(T2, row.names = FALSE)
