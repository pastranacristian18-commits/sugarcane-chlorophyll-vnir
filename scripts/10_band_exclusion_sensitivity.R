# ==============================================================================
# 10_band_exclusion_sensitivity.R -- (a) check the "reflectance < 1" reason for Day 4 / P6S1 (n = 478)
#                                    (b) do the results depend on the ~760 nm spike / last channel? (manuscript Table S10, n = 477)
# ==============================================================================
source("config.R"); source("R/functions.R"); suppressPackageStartupMessages(library(readxl))
cfg478 <- CFG; cfg478$drop_p6s1_150ddc <- FALSE                       # part (a) needs all 478 leaves
D <- load_soca(cfg = cfg478, verbose = FALSE); X <- D$X
raw <- as.data.frame(read_excel(CFG$data_raw)); rk <- paste(raw$Dia, raw$Unidad); Xr <- as.matrix(raw[, -(1:2)]); rownames(Xr) <- rk
k <- "Dia 4 P6S1"
a <- data.frame(check = c("min raw reflectance (%) of Day 4 / P6S1", "min processed (SG+MSC) reflectance (%) of Day 4 / P6S1",
                          "leaves whose raw minimum is < 1 %", "lowest raw minimum among the 478 leaves",
                          "rank of P6S1 among raw minima (1 = lowest)"),
                value = round(c(min(Xr[k, ]), min(X[k, ]), sum(apply(Xr, 1, min) < 1), min(apply(Xr, 1, min)),
                                rank(apply(Xr, 1, min))[k]), 3))
write_table(a, "table12a_P6S1_reflectance_check.csv"); print(a, row.names = FALSE)

# (b) exclude 755-767 nm and >= 815 nm -- analysed dataset of the manuscript (n = 477, default CFG), Table S10.
# "All 500 bands" = main pipeline (global SG+MSC matrix). With bands removed the MSC reference is re-estimated on the
# remaining bands inside every training set (SG-only matrix + fold-wise MSC), so the removed bands cannot enter through MSC.
D <- load_soca(verbose = FALSE); meta <- D$meta; X <- D$X; wl <- D$wl
sg <- as.data.frame(read_excel(CFG$data_sgonly)); Xs <- as.matrix(sg[, -(1:2)]); rownames(Xs) <- paste(sg$Dia, sg$Unidad); Xs <- Xs[meta$id, ]
drop <- wl >= 755 & wl <= 767 | wl >= 815; cat("\nbands removed:", sum(drop), "of", length(wl), "\n")
rows <- list()
for (s in c("random", "loco")) for (v in PIGMENTS) { y <- meta[[v]]; sp <- make_splits(meta, s)
  for (nm in c("all bands", "without 755-767 nm and >= 815 nm")) {
    pr <- if (nm == "all bands") cv_pls(X, y, meta, sp, seed = CFG$seed)$pred else cv_pls(Xs[, !drop], y, meta, sp, preprocess = pp_msc, seed = CFG$seed)$pred
    m <- metrics(y, pr); rows[[paste(s, v, nm)]] <- data.frame(pigment = unname(PIG_LABEL[v]), scheme = s, bands = nm, r2 = m$r2, R2 = m$R2, RMSE = m$RMSE, RPD = m$RPD) } }
b <- do.call(rbind, rows); rownames(b) <- NULL; b[, -(1:3)] <- round(b[, -(1:3)], 3); write_table(b, "table12b_band_exclusion.csv"); print(b, row.names = FALSE)
