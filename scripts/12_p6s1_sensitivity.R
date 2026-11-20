# ==============================================================================
# 12_p6s1_sensitivity.R -- manuscript Table S9: effect of excluding leaf P6S1 (150 DAC)
# Same five validation schemes and the same (global SG+MSC) matrix with n = 477 and n = 478.
# Adding one leaf changes the random fold assignment, so part of each difference is partition noise.
# ==============================================================================
source("config.R"); source("R/functions.R")
run <- function(drop) {
  cfg <- CFG; cfg$drop_p6s1_150ddc <- drop
  D <- load_soca(cfg = cfg, verbose = FALSE); meta <- D$meta; X <- D$X; rows <- list()
  for (s in c("random", "lopo", "lobo", "loco", "strict")) for (v in PIGMENTS) { y <- meta[[v]]
    m <- metrics(y, cv_pls(X, y, meta, make_splits(meta, s), seed = CFG$seed)$pred)
    rows[[paste(s, v)]] <- data.frame(pigment = PIG_LABEL[v], scheme = s, n = m$n, R2 = m$R2, RMSE = m$RMSE, RPD = m$RPD) }
  do.call(rbind, rows) }
a <- run(TRUE); b <- run(FALSE)
T <- data.frame(a[, c("pigment", "scheme")], R2_477 = a$R2, R2_478 = b$R2, RMSE_477 = a$RMSE, RMSE_478 = b$RMSE,
                RPD_477 = a$RPD, RPD_478 = b$RPD, dRPD = b$RPD - a$RPD)
T[, -(1:2)] <- round(T[, -(1:2)], 3); rownames(T) <- NULL
write_table(T, "table17_p6s1_sensitivity.csv"); print(T, row.names = FALSE)
