# ==============================================================================
# 04_epo_nested.R -- EPO (external parameter orthogonalisation) under LOCO-CV
#   leaky   : projector estimated from ALL campaigns, as in the original script
#             (the held-out campaign shapes the correction that is applied to it)
#   clean_k : projector estimated from the TRAINING campaigns only, fixed k = 1..4
#   nested  : k in {0,1,2,3} chosen by an inner leave-one-campaign-out on training data
# ==============================================================================
source("config.R"); source("R/functions.R")
D <- load_soca(); meta <- D$meta; X <- D$X
sp <- make_splits(meta, "loco"); rows <- list(); store <- list()
for (v in PIGMENTS) {
  y <- meta[[v]]; preds <- list()
  preds[["No EPO"]] <- cv_pls(X, y, meta, sp, seed = CFG$seed)$pred
  for (k in 1:5) { Q <- epo_Q(X, meta$campaign, k)
    preds[[paste0("EPO leaky (all campaigns), k=", k)]] <- cv_pls(epo_apply(X, Q), y, meta, sp, seed = CFG$seed)$pred }
  for (k in 1:CFG$epo_kmax) preds[[paste0("EPO training-only, k=", k)]] <- cv_pls(X, y, meta, sp, preprocess = pp_epo(k), seed = CFG$seed)$pred
  nes <- cv_pls(X, y, meta, sp, preprocess = pp_epo_nested(kmax = 3, seed = CFG$seed), seed = CFG$seed); preds[["EPO nested-selected k"]] <- nes$pred
  cat(v, "nested k per outer fold:", nes$k, "\n")
  bt <- boot_cluster(y, preds, meta$plot, B = CFG$B_boot)
  for (nm in names(preds)) { m <- metrics(y, preds[[nm]])
    rows[[paste(v, nm)]] <- cbind(pigment = PIG_LABEL[v], variant = nm, m[, c("r2", "R2", "RMSE", "RPD", "Bias")],
      R2_lo = ci(bt[, nm, "R2"])[1], R2_hi = ci(bt[, nm, "R2"])[2], RPD_lo = ci(bt[, nm, "RPD"])[1], RPD_hi = ci(bt[, nm, "RPD"])[2]) }
  store[[v]] <- list(preds = preds, nested_k = nes$k)
}
T5 <- do.call(rbind, rows); rownames(T5) <- NULL; T5[, -(1:2)] <- round(T5[, -(1:2)], 3)
write_table(T5, "table05_epo_loco.csv"); saveRDS(store, file.path(CFG$out_rds, "epo.rds")); print(T5[, 1:8], row.names = FALSE)
