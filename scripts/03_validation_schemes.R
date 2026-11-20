# ==============================================================================
# 03_validation_schemes.R -- like-for-like comparison of validation schemes
#
# Every scheme yields out-of-fold predictions for the SAME 477 leaves, so R2,
# RMSE and RPD are computed identically (pooled, same SD denominator).
#   random  : 5-fold random, stratified by campaign, ignores plots  (naive 80/20)
#   lopo    : leave-one-plot-out   (new experimental units)
#   lobo    : leave-one-block-out  (new blocks)
#   loco    : leave-one-campaign-out (new time)
#   strict  : new campaign AND new plots
# Two reference predictors are evaluated on the same splits: the training mean,
# and the training mean OF THE SAME CAMPAIGN ("just knowing the campaign").
# Uncertainty: plot-level (cluster) bootstrap -- resamples the 20 plots.
# ==============================================================================
source("config.R"); source("R/functions.R")
D <- load_soca(); meta <- D$meta; X <- D$X
schemes <- c("random", "lopo", "lobo", "loco", "strict")
SCHEME_LABEL <- c(random = "Random 5-fold (leaves)", lopo = "Leave-one-plot-out", lobo = "Leave-one-block-out",
                  loco = "Leave-one-campaign-out", strict = "New campaign + new plots")

res <- list(); rows <- list(); ratios <- list(); lvtab <- list()
for (v in PIGMENTS) {
  y <- meta[[v]]; preds <- list()
  for (s in schemes) {
    sp <- make_splits(meta, s)
    r  <- cv_pls(X, y, meta, sp, seed = CFG$seed)
    preds[[paste0("PLSR|", s)]]          <- r$pred
    preds[[paste0("Training mean|", s)]] <- cv_baseline(y, meta, sp, "mean")
    preds[[paste0("Campaign mean|", s)]] <- cv_baseline(y, meta, sp, "campaign_mean")
    lvtab[[paste(v, s)]] <- data.frame(pigment = PIG_LABEL[v], scheme = SCHEME_LABEL[s], mean_LV = mean(r$lv), min_LV = min(r$lv), max_LV = max(r$lv))
  }
  bt <- boot_cluster(y, preds, meta$plot)
  for (nm in names(preds)) {
    p <- strsplit(nm, "|", fixed = TRUE)[[1]]
    m <- metrics(y, preds[[nm]], meta$campaign)
    if (p[2] %in% c("loco", "strict")) m$R2_within <- NA        # campaign never seen in training
    rows[[paste(v, nm)]] <- cbind(pigment = PIG_LABEL[v], scheme = SCHEME_LABEL[p[2]], predictor = p[1], m,
      R2_lo = ci(bt[, nm, "R2"])[1], R2_hi = ci(bt[, nm, "R2"])[2],
      RMSE_lo = ci(bt[, nm, "RMSE"])[1], RMSE_hi = ci(bt[, nm, "RMSE"])[2],
      RPD_lo = ci(bt[, nm, "RPD"])[1], RPD_hi = ci(bt[, nm, "RPD"])[2])
  }
  # paired comparison against the naive random scheme (same resampled plots)
  for (s in setdiff(schemes, "random")) {
    a <- bt[, paste0("PLSR|", s), "RMSE"] / bt[, "PLSR|random", "RMSE"]
    b <- bt[, "PLSR|random", "RPD"] / bt[, paste0("PLSR|", s), "RPD"]
    ratios[[paste(v, s)]] <- data.frame(pigment = PIG_LABEL[v], scheme = SCHEME_LABEL[s],
      RMSE_ratio_vs_random = mean(a), lo = ci(a)[1], hi = ci(a)[2],
      RPD_inflation_of_random = mean(b), inf_lo = ci(b)[1], inf_hi = ci(b)[2])
  }
  res[[v]] <- list(preds = preds, boot = bt)
}
T3 <- do.call(rbind, rows); rownames(T3) <- NULL
T3[, sapply(T3, is.numeric)] <- round(T3[, sapply(T3, is.numeric)], 3)
write_table(T3, "table03_validation_schemes.csv")
R3 <- do.call(rbind, ratios); rownames(R3) <- NULL; R3[, -(1:2)] <- round(R3[, -(1:2)], 3)
write_table(R3, "table03b_scheme_ratios_vs_random.csv")
write_table(do.call(rbind, lvtab), "table03c_selected_latent_variables.csv")

# ---- manuscript Table 4 analogue: LOCO per fold (campaign) -------------------
fold <- list()
for (v in PIGMENTS) { y <- meta[[v]]; p <- res[[v]]$preds[["PLSR|loco"]]
  for (cm in levels(meta$campaign)) { i <- meta$campaign == cm
    m <- metrics(y[i], p[i]); fold[[paste(v, cm)]] <- cbind(pigment = PIG_LABEL[v], DDC = DDC_MAP[cm], m[, c("n", "r2", "R2", "RMSE", "RPD", "Bias")]) } }
T4 <- do.call(rbind, fold); rownames(T4) <- NULL; T4[, -(1:2)] <- round(T4[, -(1:2)], 3)
write_table(T4, "table04_loco_per_campaign.csv")

saveRDS(list(res = res), file.path(CFG$out_rds, "validation_schemes.rds"))
cat("\n== Pooled, like-for-like (PLSR only) ==\n")
print(T3[T3$predictor == "PLSR", c("pigment", "scheme", "r2", "R2", "RMSE", "RPD", "R2_within")], row.names = FALSE)
cat("\n== Reference predictors (R2) ==\n")
print(T3[T3$predictor != "PLSR", c("pigment", "scheme", "predictor", "R2", "RMSE")], row.names = FALSE)
cat("\n== Ratio vs naive random (paired cluster bootstrap) ==\n"); print(R3, row.names = FALSE)
