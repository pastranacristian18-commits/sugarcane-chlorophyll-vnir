# ==============================================================================
# 03b_monte_carlo_8020.R -- distribution of the naive random 80/20 metrics
# (repeated random splits stratified by campaign; in every split the latent variables are chosen by an
#  inner 5-fold CV grouped by plot, as in Table 2; manuscript Fig. S8).
# Slow (B_mc x 2 fits). Output: table03d + results/rds/monte_carlo_8020.rds
# ==============================================================================
source("config.R"); source("R/functions.R")
D <- load_soca(); meta <- D$meta; X <- D$X
set.seed(CFG$seed + 100); mc <- list()
for (v in PIGMENTS) { y <- meta[[v]]; out <- matrix(NA_real_, CFG$B_mc, 4, dimnames = list(NULL, c("r2", "R2", "RMSE", "RPD")))
  for (b in seq_len(CFG$B_mc)) {
    ic <- unlist(lapply(levels(meta$campaign), function(cm) { i <- which(meta$campaign == cm); sample(i, floor(0.8 * length(i))) }))
    iv <- setdiff(seq_along(y), ic)
    f <- fit_pls(X[ic, ], y[ic], groups = meta$plot[ic], seed = b); p <- predict_pls(f, X[iv, ])
    mm <- metrics(y[iv], p); out[b, ] <- unlist(mm[, c("r2", "R2", "RMSE", "RPD")]) }
  mc[[v]] <- out }
mct <- do.call(rbind, lapply(PIGMENTS, function(v) data.frame(pigment = PIG_LABEL[v],
  metric = colnames(mc[[v]]), median = apply(mc[[v]], 2, median), lo = apply(mc[[v]], 2, function(z) ci(z)[1]), hi = apply(mc[[v]], 2, function(z) ci(z)[2]))))
mct[, -(1:2)] <- round(mct[, -(1:2)], 3); write_table(mct, "table03d_monte_carlo_random_8020.csv")

saveRDS(mc, file.path(CFG$out_rds, "monte_carlo_8020.rds"))
print(mct, row.names = FALSE)
