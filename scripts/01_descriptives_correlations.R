# ==============================================================================
# 01_descriptives_correlations.R -- manuscript Table 1 and Pearson profiles (Fig. 2)
# ==============================================================================
source("config.R"); source("R/functions.R")
D <- load_soca(); meta <- D$meta; X <- D$X; wl <- D$wl

# ---- Table 1: descriptive statistics per campaign -----------------------------
desc <- function(x) data.frame(n = length(x), Mean = mean(x), Min = min(x), Max = max(x), SD = sd(x), CV_pct = 100 * sd(x) / mean(x))
T1 <- do.call(rbind, lapply(PIGMENTS, function(v) do.call(rbind, lapply(levels(meta$campaign), function(cm)
  cbind(pigment = PIG_LABEL[v], DAC = DDC_MAP[cm], desc(meta[[v]][meta$campaign == cm]))))))
T1[, -(1:2)] <- round(T1[, -(1:2)], 3); write_table(T1, "table01_descriptives_by_campaign.csv")

# ---- Pearson r(reflectance, pigment) per campaign, LOESS-smoothed -------------
span <- 25 / diff(range(wl)); prof <- list()
for (cm in levels(meta$campaign)) { i <- meta$campaign == cm
  for (v in PIGMENTS) { r <- as.vector(cor(X[i, ], meta[[v]][i]))
    sm <- suppressWarnings(predict(loess(r ~ wl, span = span)))
    prof[[paste(cm, v)]] <- data.frame(DAC = DDC_MAP[cm], pigment = PIG_LABEL[v], wavelength = wl, r = r, r_smooth = sm) } }
prof <- do.call(rbind, prof); saveRDS(prof, file.path(CFG$out_rds, "pearson_profiles.rds"))
# Critical |r| for p = 0.05 (two-sided): 80 leaves (naive) vs 20 plots (independent units)
crit <- function(n) { t <- qt(.975, n - 2); t / sqrt(n - 2 + t^2) }
cat(sprintf("critical |r|: n=80 leaves -> %.2f | n=20 plots -> %.2f\n", crit(80), crit(20)))
# Peak |r| per campaign (smoothed)
pk <- aggregate(abs(r_smooth) ~ DAC + pigment, prof, max); names(pk)[3] <- "max_abs_r_smoothed"
write_table(transform(pk, max_abs_r_smoothed = round(max_abs_r_smoothed, 3)), "table01b_peak_correlation_by_campaign.csv")
print(T1, row.names = FALSE)
