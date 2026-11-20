# ==============================================================================
# 08_msc_sensitivity.R -- does estimating the MSC reference on ALL spectra leak information?
#   global : SG+MSC matrix exported by MATLAB (MSC reference = mean of all 478 spectra)
#   fold   : SG-only matrix + MSC re-estimated inside every training fold
# Requires the 'readxl' package to read the MATLAB export.
# ==============================================================================
source("config.R"); source("R/functions.R"); suppressPackageStartupMessages(library(readxl))
D <- load_soca(); meta <- D$meta; Xg <- D$X
sg <- as.data.frame(read_excel(CFG$data_sgonly)); key <- paste(sg$Dia, sg$Unidad)
Xs <- as.matrix(sg[, -(1:2)]); rownames(Xs) <- key; Xs <- Xs[meta$id, ]           # align by key, never by position
stopifnot(all(rownames(Xs) == meta$id), ncol(Xs) == ncol(Xg))
rows <- list()
for (s in c("loco", "random")) for (v in PIGMENTS) { y <- meta[[v]]; sp <- make_splits(meta, s)
  a <- cv_pls(Xg, y, meta, sp, seed = CFG$seed)$pred; b <- cv_pls(Xs, y, meta, sp, preprocess = pp_msc, seed = CFG$seed)$pred
  ma <- metrics(y, a); mb <- metrics(y, b)
  rows[[paste(s, v)]] <- data.frame(pigment = PIG_LABEL[v], scheme = s, r2_global = ma$r2, r2_fold = mb$r2, R2_global = ma$R2, R2_fold = mb$R2,
    RMSE_global = ma$RMSE, RMSE_fold = mb$RMSE, delta_RMSE_pct = 100 * (mb$RMSE - ma$RMSE) / ma$RMSE) }
T <- do.call(rbind, rows); rownames(T) <- NULL; T[, -(1:2)] <- round(T[, -(1:2)], 3)
write_table(T, "table11_msc_leakage_sensitivity.csv"); print(T, row.names = FALSE)
