# ==============================================================================
# 05_within_campaign_skill.R -- do the spectra explain chlorophyll BEYOND phenology?
#  (a) campaign-centred PLSR, leave-one-plot-out  (skill relative to "knowing the campaign mean")
#  (b) campaign-centred PLSR, leave-one-campaign-out ("offset-free" transfer to a new campaign)
#  (c) stage-specific (per-campaign) models, plot-grouped 5-fold CV + plot-level bootstrap
#  (d) leave-one-stage-out (LOSO): tillering (60-120 DDC) <-> elongation (150-210 DDC)
# ==============================================================================
source("config.R"); source("R/functions.R")
D <- load_soca(); meta <- D$meta; X <- D$X; out <- list()
rows_a <- list(); rows_c <- list(); rows_d <- list(); pred_c <- list()

for (v in PIGMENTS) {
  y <- meta[[v]]
  # (a) + (b)
  a <- cv_pls_centered(X, y, meta, make_splits(meta, "lopo"), "train", seed = CFG$seed)
  b <- cv_pls_centered(X, y, meta, make_splits(meta, "loco"), "self",  seed = CFG$seed)
  rows_a[[v]] <- rbind(cbind(pigment = PIG_LABEL[v], scheme = "Leave-one-plot-out, campaign-centred", dev_metrics(a$dev_obs, a$dev_pred)),
                       cbind(pigment = PIG_LABEL[v], scheme = "Leave-one-campaign-out, campaign-centred (offset removed)", dev_metrics(b$dev_obs, b$dev_pred)))
  for (cm in levels(meta$campaign)) { i <- meta$campaign == cm
    rows_a[[paste(v, cm)]] <- cbind(pigment = PIG_LABEL[v], scheme = paste0("  LOCO centred, held-out ", DDC_MAP[cm], " DDC"), dev_metrics(b$dev_obs[i], b$dev_pred[i])) }
  out[[v]] <- list(lopo_centered = a, loco_centered = b)

  # (c) stage-specific models: one model per campaign, validated on held-out plots
  for (cm in levels(meta$campaign)) {
    ii <- which(meta$campaign == cm); mm <- meta[ii, ]; yy <- y[ii]; XX <- X[ii, ]
    sp <- make_splits(mm, "lopo"); grp <- split(seq_along(sp), rep(seq_len(CFG$k_random), length.out = length(sp)))   # 5 folds of 4 plots
    sp5 <- lapply(grp, function(g) list(train = setdiff(seq_len(nrow(mm)), unlist(lapply(sp[g], `[[`, "test"))), test = unlist(lapply(sp[g], `[[`, "test"))))
    r <- cv_pls(XX, yy, mm, sp5, seed = CFG$seed)
    bt <- boot_cluster(yy, list(m = r$pred), mm$plot, B = CFG$B_boot); m <- metrics(yy, r$pred)
    rows_c[[paste(v, cm)]] <- cbind(pigment = PIG_LABEL[v], DDC = DDC_MAP[cm], m[, c("n", "r2", "R2", "RMSE", "RPD", "Bias")],
      R2_lo = ci(bt[, "m", "R2"])[1], R2_hi = ci(bt[, "m", "R2"])[2], RPD_lo = ci(bt[, "m", "RPD"])[1], RPD_hi = ci(bt[, "m", "RPD"])[2], mean_LV = mean(r$lv))
    pred_c[[paste(v, cm)]] <- data.frame(idx = ii, pred = r$pred)
  }
  # (d) LOSO
  stage <- ifelse(meta$DDC <= 120, "Tillering", "Elongation")
  for (tr_s in c("Tillering", "Elongation")) { te_s <- setdiff(unique(stage), tr_s); tr <- which(stage == tr_s); te <- which(stage == te_s)
    f <- fit_pls(X[tr, ], y[tr], groups = meta$plot[tr], seed = CFG$seed); p <- predict_pls(f, X[te, ])
    rows_d[[paste(v, tr_s)]] <- cbind(pigment = PIG_LABEL[v], train = tr_s, test = te_s, LV = f$lv, metrics(y[te], p)[, c("n", "r2", "R2", "RMSE", "RPD", "Bias")]) }
}
Ta <- do.call(rbind, rows_a); Tc <- do.call(rbind, rows_c); Td <- do.call(rbind, rows_d)
for (n in c("Ta", "Tc", "Td")) { z <- get(n); z[sapply(z, is.numeric)] <- round(z[sapply(z, is.numeric)], 3); assign(n, z) }
rownames(Ta) <- rownames(Tc) <- rownames(Td) <- NULL
write_table(Ta, "table06_within_campaign_skill.csv"); write_table(Tc, "table07_stage_specific_models.csv"); write_table(Td, "table08_loso_stages.csv")
saveRDS(list(centered = out, stage_pred = pred_c), file.path(CFG$out_rds, "within_campaign.rds"))
cat("\n(a,b) within-campaign skill (R2_dev>0 means better than knowing the campaign mean):\n"); print(Ta, row.names = FALSE)
cat("\n(c) stage-specific models, plot-grouped CV:\n"); print(Tc, row.names = FALSE)
cat("\n(d) LOSO:\n"); print(Td, row.names = FALSE)
