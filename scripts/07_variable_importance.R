# ==============================================================================
# 07_variable_importance.R -- PLSR VIP and PCA-RF permutation importance (Fig. 6)
# Model = the manuscript's random 80/20 calibration set (seed 42), LV by LOO.
# The PCA-RF hyper-parameter grid is the one implemented in the original code.
# ==============================================================================
source("config.R"); source("R/functions.R"); suppressPackageStartupMessages(library(randomForest))
D <- load_soca(); meta <- D$meta; X <- D$X; wl <- D$wl
fits <- readRDS(file.path(CFG$out_rds, "random_8020_fits.rds"))

pcarf_grid <- function(Xc, yc, pcs = c(5, 8, 10, 12, 15), trees = c(100, 200, 300, 500), k = 5, seed = 42) {
  set.seed(seed); f <- sample(rep(1:k, length.out = nrow(Xc))); pca <- prcomp(Xc, center = TRUE, scale. = TRUE)
  best <- list(score = -Inf)
  for (np in pcs) { S <- pca$x[, 1:np, drop = FALSE]
    for (nt in trees) for (mt in unique(c(floor(np / 3), floor(np / 2), floor(sqrt(np)), np))) {
      r <- sapply(1:k, function(j) { m <- randomForest(S[f != j, , drop = FALSE], yc[f != j], ntree = nt, mtry = mt)
                                      sqrt(mean((yc[f == j] - predict(m, S[f == j, , drop = FALSE]))^2)) })
      if (-mean(r) > best$score) best <- list(score = -mean(r), np = np, nt = nt, mt = mt) } }
  best }

out <- list(); summ <- list()
for (v in PIGMENTS) {
  ic <- fits[[v]]$idx_cal; iv <- fits[[v]]$idx_val; yc <- meta[[v]][ic]; lv <- fits[[v]]$lv
  dc <- data.frame(y = yc); dc$X <- I(X[ic, ]); m <- plsr(y ~ X, ncomp = lv, data = dc, method = "oscorespls")
  vip <- as.vector(VIP(m)[, lv])
  bp <- pcarf_grid(X[ic, ], yc); pca <- prcomp(X[ic, ], center = TRUE, scale. = TRUE)
  rf <- randomForest(pca$x[, 1:bp$np, drop = FALSE], yc, ntree = bp$nt, mtry = bp$mt, importance = TRUE)
  imp <- as.vector(abs(pca$rotation[, 1:bp$np, drop = FALSE]) %*% importance(rf)[, 1]); imp <- imp / max(imp)
  pv <- predict(rf, predict(pca, X[iv, ])[, 1:bp$np, drop = FALSE]); pc <- predict(rf, pca$x[, 1:bp$np, drop = FALSE])
  cat(v, "PCA-RF: PCs =", bp$np, "ntree =", bp$nt, "mtry =", bp$mt, "| cal r2 =", round(cor(yc, pc)^2, 3), "val r2 =", round(cor(meta[[v]][iv], pv)^2, 3), "\n")
  out[[v]] <- data.frame(pigment = unname(PIG_LABEL[v]), wavelength = wl, VIP = vip, RF_importance = imp)
  summ[[v]] <- list(rf = c(cal_r2 = cor(yc, pc)^2, val_r2 = cor(meta[[v]][iv], pv)^2, np = bp$np, nt = bp$nt, mt = bp$mt))
}
imp_all <- do.call(rbind, out); saveRDS(imp_all, file.path(CFG$out_rds, "importance_profiles.rds"))
reg <- data.frame(region = c("450-520 (blue-green)", "520-600 (green)", "600-700 (red absorption)", "700-770 (red edge)", "770-820 (NIR plateau)"),
                  lo = c(450, 520, 600, 700, 770), hi = c(520, 600, 700, 770, 820))
T6 <- do.call(rbind, lapply(PIGMENTS, function(v) { d <- out[[v]]
  do.call(rbind, lapply(seq_len(nrow(reg)), function(i) { s <- d$wavelength >= reg$lo[i] & d$wavelength < reg$hi[i] + (i == 5)
    data.frame(pigment = PIG_LABEL[v], region = reg$region[i], mean_VIP = mean(d$VIP[s]), pct_bands_VIP_gt1 = 100 * mean(d$VIP[s] > 1), mean_RF_importance = mean(d$RF_importance[s])) })) }))
T6[, -(1:2)] <- round(T6[, -(1:2)], 3); write_table(T6, "table10_spectral_importance_by_region.csv")
cmp <- data.frame(quantity = c("cor(VIP Chl a, VIP Chl b)", "cor(RF Chl a, RF Chl b)", "PCA-RF r2 calibration Chl a", "PCA-RF r2 validation Chl a"),
  value = round(c(cor(out$Chla$VIP, out$Chlb$VIP), cor(out$Chla$RF_importance, out$Chlb$RF_importance), summ$Chla$rf["cal_r2"], summ$Chla$rf["val_r2"]), 3))
write_table(cmp, "table10b_importance_checks.csv"); print(T6, row.names = FALSE); print(cmp, row.names = FALSE)
