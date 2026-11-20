# ==============================================================================
# 13_design_dependence_180ddc.R -- manuscript Table S4 / Fig. S5: at 180 DAC (the campaign with the strongest
# within-campaign correlations) the apparent skill of a campaign-specific model depends on the validation design:
#   (1) spectral-extremes split labelled "Kennard-Stone": the 20 % most extreme spectra (max-min distance on PC1-5)
#       form the validation set. The split is deterministic; only the inner-CV seed (1-20) changes the selected LV.
#   (2) random leaf-level 5-fold CV (LV chosen by ungrouped CV)  -> same as table09 "naive_leaf"
#   (3) whole plots held out: 50 random plot-to-fold assignments (5 folds of 4 plots), LV chosen by plot-grouped inner CV
# ==============================================================================
source("config.R"); source("R/functions.R")
D <- load_soca(); meta <- D$meta; X <- D$X
src <- readLines("scripts/06_diagnose_stage_models.R")            # reuse the same ks_extreme() definition
i0 <- grep("^ks_extreme", src); eval(parse(text = src[i0:(i0 + 4)]))
ii <- which(meta$DDC == 180); mm <- meta[ii, ]; XX <- X[ii, ]; pl <- sort(unique(mm$plot))
q3 <- function(z) c(median(z), quantile(z, c(0.025, 0.975)))
rows <- list()
for (v in PIGMENTS) { yy <- meta[[v]][ii]
  # (1) spectral-extremes split, 20 inner-CV seeds
  iv <- ks_extreme(XX, floor(0.2 * nrow(XX))); ic <- setdiff(seq_len(nrow(XX)), iv)
  ks <- sapply(1:20, function(s) { f <- fit_pls(XX[ic, ], yy[ic], groups = NULL, seed = s); c(unlist(metrics(yy[iv], predict_pls(f, XX[iv, ]))[c("r2", "R2", "RPD")]), LV = unname(f$lv)) })
  lvt <- table(ks["LV", ]); top <- names(lvt)[which.max(lvt)]
  rows[[paste(v, "ks")]] <- data.frame(pigment = PIG_LABEL[v], design = "Spectral-extremes split (Kennard-Stone label), 20 inner-CV seeds",
    r2 = sprintf("%.3f-%.3f", min(ks["r2", ]), max(ks["r2", ])), R2 = sprintf("%.3f-%.3f", min(ks["R2", ]), max(ks["R2", ])),
    RPD = sprintf("%.3f-%.3f", min(ks["RPD", ]), max(ks["RPD", ])), note = sprintf("%s LV in %d/20 seeds; max r2 = %.3f", top, max(lvt), max(ks["r2", ])))
  # (2) random leaf-level 5-fold
  set.seed(CFG$seed); f5 <- sample(rep(1:5, length.out = nrow(XX))); p <- rep(NA_real_, nrow(XX))
  for (k in 1:5) { tr <- which(f5 != k); te <- which(f5 == k); f <- fit_pls(XX[tr, ], yy[tr], groups = NULL, seed = 1 + k); p[te] <- predict_pls(f, XX[te, ]) }
  m <- metrics(yy, p)
  rows[[paste(v, "rand")]] <- data.frame(pigment = PIG_LABEL[v], design = "Random leaf-level 5-fold", r2 = sprintf("%.3f", m$r2), R2 = sprintf("%.3f", m$R2),
    RPD = sprintf("%.3f", m$RPD), note = "")
  # (3) whole plots held out, 50 partitions
  res <- t(sapply(1:50, function(s) { set.seed(CFG$seed + s); g <- sample(rep(1:5, length.out = length(pl)))[match(mm$plot, pl)]
    p <- rep(NA_real_, length(yy))
    for (k in 1:5) { tr <- which(g != k); te <- which(g == k); f <- fit_pls(XX[tr, ], yy[tr], groups = mm$plot[tr], seed = CFG$seed + k); p[te] <- predict_pls(f, XX[te, ]) }
    m <- metrics(yy, p); c(m$r2, m$R2, m$RPD) }))
  fm <- function(z) sprintf("%.3f [%.3f, %.3f]", q3(z)[1], q3(z)[2], q3(z)[3])
  rows[[paste(v, "plots")]] <- data.frame(pigment = PIG_LABEL[v], design = "Whole plots held out, median [2.5th, 97.5th] of 50 partitions",
    r2 = fm(res[, 1]), R2 = fm(res[, 2]), RPD = fm(res[, 3]), note = sprintf("R2 < 0 in %d of 50", sum(res[, 2] < 0))) }
T <- do.call(rbind, rows); rownames(T) <- NULL
write_table(T, "table18_design_dependence_180ddc.csv"); print(T, row.names = FALSE)
