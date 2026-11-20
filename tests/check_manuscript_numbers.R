# ==============================================================================
# check_manuscript_numbers.R -- compare the tables written by scripts/run_all.R with the values printed
# in the manuscript / Supporting Information. Run AFTER run_all.R:  Rscript tests/check_manuscript_numbers.R
# Tolerance 0.0025 absolute (the manuscript reports 2-3 decimals). Prints OK / DIFF for every value.
# ==============================================================================
rd <- function(f) read.csv(file.path("results/tables", f), check.names = FALSE, stringsAsFactors = FALSE)
res <- list()
chk <- function(label, got, want, tol = 0.0025) {
  ok <- is.finite(got) && abs(got - want) <= tol
  res[[length(res) + 1]] <<- data.frame(item = label, manuscript = want, repo = round(got, 4), status = ifelse(ok, "OK", "DIFF")) }
pick <- function(d, ..., col) { cond <- list(...); i <- rep(TRUE, nrow(d)); for (n in names(cond)) i <- i & d[[n]] == cond[[n]]; d[[col]][i][1] }

# Table 2 (random 80/20)
t2 <- rd("table02_plsr_random_8020.csv")
for (p in c("Chl a", "Chl b")) { v <- if (p == "Chl a") c(0.494, 0.484, 0.331, 1.399, 2.431, -0.031) else c(0.464, 0.445, 0.140, 1.349, 2.081, -0.008)
  for (k in seq_along(c("r2", "R2", "RMSE", "RPD", "RPIQ", "Bias"))) { cn <- c("r2", "R2", "RMSE", "RPD", "RPIQ", "Bias")[k]
    chk(paste("Table 2", p, "validation", cn), pick(t2, pigment = p, set = "Validation 80/20", col = cn), v[k]) }
  chk(paste("Table 2", p, "LV"), pick(t2, pigment = p, set = "Validation 80/20", col = "LV"), 11, 0) }

# Table 3 (pooled schemes) and Table S2
t3 <- rd("table03_validation_schemes.csv"); t3 <- t3[t3$predictor == "PLSR", ]
spec <- list(
  list("Chl a", "Random 5-fold (leaves)", 0.458, 0.353, 1.360, 2.173), list("Chl a", "Leave-one-plot-out", 0.431, 0.361, 1.327, 2.120),
  list("Chl a", "Leave-one-campaign-out", 0.029, 0.472, 1.016, 1.623), list("Chl b", "Random 5-fold (leaves)", 0.395, 0.160, 1.287, 1.799),
  list("Chl b", "Leave-one-plot-out", 0.381, 0.162, 1.272, 1.778), list("Chl b", "Leave-one-campaign-out", -0.040, 0.210, 0.982, 1.372),
  list("Chl a", "Leave-one-block-out", 0.397, 0.372, 1.289, 2.060), list("Chl a", "New campaign + new plots", -0.027, 0.485, 0.988, 1.578),
  list("Chl b", "Leave-one-block-out", 0.341, 0.167, 1.233, 1.724), list("Chl b", "New campaign + new plots", -0.059, 0.211, 0.973, 1.359))
for (s in spec) for (k in 1:4) { cn <- c("R2", "RMSE", "RPD", "RPIQ")[k]
  chk(paste("Table 3/S2", s[[1]], s[[2]], cn), pick(t3, pigment = s[[1]], scheme = s[[2]], col = cn), s[[2 + k]]) }
rr <- rd("table03b_scheme_ratios_vs_random.csv")
chk("Table 3 RMSE ratio LOCO Chl a", pick(rr, pigment = "Chl a", scheme = "Leave-one-campaign-out", col = "RMSE_ratio_vs_random"), 1.338, 0.0025)
chk("Table 3 RMSE ratio LOCO Chl b", pick(rr, pigment = "Chl b", scheme = "Leave-one-campaign-out", col = "RMSE_ratio_vs_random"), 1.311)
# campaign-mean benchmark (text: RPD 1.59, R2 0.60, Chl a random)
chk("Campaign-mean RPD Chl a random", pick(rd("table03_validation_schemes.csv"), pigment = "Chl a", scheme = "Random 5-fold (leaves)", predictor = "Campaign mean", col = "RPD"), 1.59, 0.005)
chk("Training-mean RMSE LOCO Chl a", pick(rd("table03_validation_schemes.csv"), pigment = "Chl a", scheme = "Leave-one-campaign-out", predictor = "Training mean", col = "RMSE"), 0.540, 0.001)

# Table S3 (LOCO by campaign; Chl a RPD)
s3 <- rd("table04_loco_per_campaign.csv")
for (i in seq_len(6)) chk(paste("Table S3 Chl a RPD", c(60, 90, 120, 150, 180, 210)[i]), s3$RPD[s3$DDC == c(60, 90, 120, 150, 180, 210)[i] & s3$pigment == "Chl a"], c(0.547, 0.780, 0.668, 0.680, 0.795, 0.483)[i])


chk_str <- function(label, got, want) res[[length(res) + 1]] <<- data.frame(item = label, manuscript = want, repo = got, status = ifelse(identical(got, want), "OK", "DIFF"))

# Table 1 (means and SD by campaign)
t1 <- rd("table01_descriptives_by_campaign.csv"); dd <- c(60, 90, 120, 150, 180, 210)
for (i in 1:6) { chk(paste("Table 1 Chl a mean", dd[i]), t1$Mean[t1$pigment == "Chl a" & t1$DAC == dd[i]], c(1.743, 1.407, 0.976, 0.697, 0.877, 0.767)[i], 0.0006)
  chk(paste("Table 1 Chl a SD", dd[i]), t1$SD[t1$pigment == "Chl a" & t1$DAC == dd[i]], c(0.320, 0.253, 0.296, 0.345, 0.326, 0.250)[i], 0.0006)
  chk(paste("Table 1 Chl b mean", dd[i]), t1$Mean[t1$pigment == "Chl b" & t1$DAC == dd[i]], c(0.694, 0.520, 0.350, 0.239, 0.309, 0.284)[i], 0.0006) }

# 1000 repeated random 80/20 splits (text of Results; Fig. S8)
mc <- rd("table03d_monte_carlo_random_8020.csv")
chk("Monte-Carlo median RPD Chl a", pick(mc, pigment = "Chl a", metric = "RPD", col = "median"), 1.36, 0.0055)
chk("Monte-Carlo 2.5% RPD Chl a", pick(mc, pigment = "Chl a", metric = "RPD", col = "lo"), 1.21, 0.0055)
chk("Monte-Carlo 97.5% RPD Chl a", pick(mc, pigment = "Chl a", metric = "RPD", col = "hi"), 1.54, 0.0055)
chk("Monte-Carlo median RPD Chl b", pick(mc, pigment = "Chl b", metric = "RPD", col = "median"), 1.30, 0.0055)
chk("Monte-Carlo 2.5% RPD Chl b", pick(mc, pigment = "Chl b", metric = "RPD", col = "lo"), 1.17, 0.0055)
chk("Monte-Carlo 97.5% RPD Chl b", pick(mc, pigment = "Chl b", metric = "RPD", col = "hi"), 1.46, 0.0055)

# Table 4 (EPO)
e5 <- rd("table05_epo_loco.csv")
ev <- list(list("Chl a", "No EPO", 0.134, 0.029, 0.472, 1.016), list("Chl a", "EPO leaky (all campaigns), k=5", 0.416, -0.274, 0.541, 0.887),
  list("Chl a", "EPO training-only, k=4", 0.388, -0.271, 0.540, 0.888), list("Chl a", "EPO nested-selected k", 0.053, -0.132, 0.510, 0.941),
  list("Chl b", "No EPO", 0.081, -0.040, 0.210, 0.982), list("Chl b", "EPO leaky (all campaigns), k=5", 0.456, -0.269, 0.232, 0.889),
  list("Chl b", "EPO training-only, k=4", 0.430, -0.261, 0.231, 0.891), list("Chl b", "EPO nested-selected k", 0.027, -0.176, 0.223, 0.923))
for (e in ev) for (k in 1:4) { cn <- c("r2", "R2", "RMSE", "RPD")[k]; chk(paste("Table 4", e[[1]], e[[2]], cn), pick(e5, pigment = e[[1]], variant = e[[2]], col = cn), e[[2 + k]]) }

# Tables S5 (ICC), S6 (LOSO, tolerance 0.005), S8, S9, S10
t9 <- rd("table09_stage_model_diagnosis.csv")
for (i in 1:6) { chk(paste("Table S5 ICC Chl a", dd[i]), t9$ICC_plot[t9$pigment == "Chl a" & t9$DDC == dd[i]], c(0.394, 0.015, 0.298, 0.041, 0.535, 0.395)[i], 0.0006)
  chk(paste("Table S5 ICC Chl b", dd[i]), t9$ICC_plot[t9$pigment == "Chl b" & t9$DDC == dd[i]], c(0.395, 0.104, 0.349, 0.014, 0.505, 0.395)[i], 0.0006) }
t8 <- rd("table08_loso_stages.csv")
chk("Table S6 R2 Chl a Till->Elong", pick(t8, pigment = "Chl a", train = "Tillering", col = "R2"), -1.324, 0.005); chk("Table S6 R2 Chl a Elong->Till", pick(t8, pigment = "Chl a", train = "Elongation", col = "R2"), -4.772, 0.005)
chk("Table S6 R2 Chl b Till->Elong", pick(t8, pigment = "Chl b", train = "Tillering", col = "R2"), -1.668, 0.005); chk("Table S6 R2 Chl b Elong->Till", pick(t8, pigment = "Chl b", train = "Elongation", col = "R2"), -2.957, 0.005)
m8 <- rd("table11_msc_leakage_sensitivity.csv")
chk("Table S8 dRMSE% LOCO Chl a", pick(m8, pigment = "Chl a", scheme = "loco", col = "delta_RMSE_pct"), 0.006, 0.0006); chk("Table S8 dRMSE% LOCO Chl b", pick(m8, pigment = "Chl b", scheme = "loco", col = "delta_RMSE_pct"), 0.004, 0.0006)
t17 <- rd("table17_p6s1_sensitivity.csv")
for (r in list(list("Chl a", "random", 1.360, 1.357), list("Chl a", "lopo", 1.327, 1.309), list("Chl a", "loco", 1.016, 1.014), list("Chl b", "loco", 0.982, 0.996), list("Chl b", "strict", 0.973, 0.986))) {
  chk(paste("Table S9 RPD 477", r[[1]], r[[2]]), pick(t17, pigment = r[[1]], scheme = r[[2]], col = "RPD_477"), r[[3]]); chk(paste("Table S9 RPD 478", r[[1]], r[[2]]), pick(t17, pigment = r[[1]], scheme = r[[2]], col = "RPD_478"), r[[4]]) }
t12 <- rd("table12b_band_exclusion.csv"); bs <- "without 755-767 nm and >= 815 nm"
for (r in list(list("Chl a", "random", 1.374), list("Chl b", "random", 1.338), list("Chl a", "loco", 1.026), list("Chl b", "loco", 1.000)))
  chk(paste("Table S10 RPD bands removed", r[[1]], r[[2]]), pick(t12, pigment = r[[1]], scheme = r[[2]], bands = bs, col = "RPD"), r[[3]])

# Table S4 (180 DAC designs; text strings)
t18 <- rd("table18_design_dependence_180ddc.csv")
chk_str("Table S4 Chl a whole plots RPD", pick(t18, pigment = "Chl a", design = "Whole plots held out, median [2.5th, 97.5th] of 50 partitions", col = "RPD"), "0.955 [0.891, 1.063]")
chk_str("Table S4 Chl a whole plots r2", pick(t18, pigment = "Chl a", design = "Whole plots held out, median [2.5th, 97.5th] of 50 partitions", col = "r2"), "0.086 [0.025, 0.200]")
chk_str("Table S4 Chl a spectral-extremes r2 range", pick(t18, pigment = "Chl a", design = "Spectral-extremes split (Kennard-Stone label), 20 inner-CV seeds", col = "r2"), "0.283-0.635")
chk_str("Table S4 Chl b spectral-extremes r2 range", pick(t18, pigment = "Chl b", design = "Spectral-extremes split (Kennard-Stone label), 20 inner-CV seeds", col = "r2"), "0.295-0.709")

out <- do.call(rbind, res); print(out, row.names = FALSE)
cat(sprintf("\n%d values checked: %d OK, %d DIFF\n", nrow(out), sum(out$status == "OK"), sum(out$status == "DIFF")))
if (any(out$status == "DIFF")) quit(status = 1)
