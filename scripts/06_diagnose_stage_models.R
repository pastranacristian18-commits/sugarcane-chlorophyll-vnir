# ==============================================================================
# 06_diagnose_stage_models.R -- why the per-campaign ("stage") result changes with the validation design
#   (1) intra-plot correlation (ICC): how similar are the 4 leaves of one plot?
#   (2) naive   : random leaf-level 5-fold, LV chosen by ungrouped CV
#   (3) grouped : plot-grouped 5-fold outer, LV chosen by ungrouped CV
#   (4) grouped : plot-grouped 5-fold outer, LV chosen by plot-grouped CV   (recommended; = table07)
#   (5) manuscript procedure: the 20% most spectrally extreme leaves (max-min "Kennard-Stone" on PC1-5)
#       are used for validation, the rest for calibration (median over inner-CV seeds)
# ==============================================================================
source("config.R"); source("R/functions.R")
D <- load_soca(); meta <- D$meta; X <- D$X

icc <- function(y, g) { g <- factor(g); k <- 4; a <- anova(lm(y ~ g)); msb <- a[1, 3]; msw <- a[2, 3]; (msb - msw) / (msb + (k - 1) * msw) }
ks_extreme <- function(Xs, nv) {                                   # reproduces the original script (selection goes to VALIDATION)
  sc <- prcomp(Xs, center = TRUE, scale. = TRUE)$x[, 1:5]; dm <- as.matrix(dist(sc)); cen <- colMeans(sc)
  ks <- which.max(rowSums(sweep(sc, 2, cen)^2)); rem <- setdiff(seq_len(nrow(Xs)), ks)
  while (length(ks) < nv) { md <- apply(dm[rem, ks, drop = FALSE], 1, min); nx <- rem[which.max(md)]; ks <- c(ks, nx); rem <- setdiff(rem, nx) }
  ks }
grouped5 <- function(mm) { pl <- sort(unique(mm$plot)); set.seed(CFG$seed); g <- sample(rep(1:5, length.out = length(pl)))[match(mm$plot, pl)]
  lapply(1:5, function(f) list(train = which(g != f), test = which(g == f))) }
random5 <- function(n) { set.seed(CFG$seed); f <- sample(rep(1:5, length.out = n)); lapply(1:5, function(k) list(train = which(f != k), test = which(f == k))) }
run_cv <- function(XX, yy, mm, sp, grouped_inner, seed = 1) {
  p <- rep(NA_real_, length(yy))
  for (i in seq_along(sp)) { f <- fit_pls(XX[sp[[i]]$train, ], yy[sp[[i]]$train], groups = if (grouped_inner) mm$plot[sp[[i]]$train] else NULL, seed = seed + i)
    p[sp[[i]]$test] <- predict_pls(f, XX[sp[[i]]$test, ]) }
  p }

rows <- list()
for (v in PIGMENTS) for (cm in levels(meta$campaign)) {
  ii <- which(meta$campaign == cm); mm <- meta[ii, ]; yy <- meta[[v]][ii]; XX <- X[ii, ]
  m2 <- metrics(yy, run_cv(XX, yy, mm, random5(nrow(mm)), FALSE))
  m3 <- metrics(yy, run_cv(XX, yy, mm, grouped5(mm),  FALSE))
  m4 <- metrics(yy, run_cv(XX, yy, mm, grouped5(mm),  TRUE))
  ks <- sapply(1:20, function(s) { set.seed(s); iv <- ks_extreme(XX, max(2, floor(0.2 * nrow(XX)))); ic <- setdiff(seq_len(nrow(XX)), iv)
    f <- fit_pls(XX[ic, ], yy[ic], groups = NULL, max_lv = 12, seed = s); unlist(metrics(yy[iv], predict_pls(f, XX[iv, ]))[c("r2", "R2", "RPD")]) })
  rows[[paste(v, cm)]] <- data.frame(pigment = PIG_LABEL[v], DDC = DDC_MAP[cm], ICC_plot = icc(yy, mm$plot),
    naive_leaf_r2 = m2$r2, naive_leaf_R2 = m2$R2, plotCV_ungroupedLV_r2 = m3$r2, plotCV_ungroupedLV_R2 = m3$R2,
    plotCV_groupedLV_r2 = m4$r2, plotCV_groupedLV_R2 = m4$R2,
    manuscript_KS_r2_median = median(ks["r2", ]), manuscript_KS_R2_median = median(ks["R2", ]), manuscript_KS_RPD_median = median(ks["RPD", ]))
}
T <- do.call(rbind, rows); rownames(T) <- NULL; T[, -(1:2)] <- round(T[, -(1:2)], 3)
write_table(T, "table09_stage_model_diagnosis.csv"); print(T, row.names = FALSE, width = 200)
