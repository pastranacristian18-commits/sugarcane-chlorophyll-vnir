# ==============================================================================
# run_all.R -- run the complete analysis from the repository root:
#     Rscript scripts/run_all.R
# Approximate run time on one CPU core: ~35 min (the Monte-Carlo step dominates).
# ==============================================================================
steps <- c("00_data_audit", "01_descriptives_correlations", "02_random_split_plsr",
           "03_validation_schemes", "03b_monte_carlo_8020", "04_epo_nested",
           "05_within_campaign_skill", "06_diagnose_stage_models", "07_variable_importance",
           "08_msc_sensitivity", "10_band_exclusion_sensitivity", "11_analisis_adicionales",
           "12_p6s1_sensitivity", "13_design_dependence_180ddc", "09_figures", "14_toc_graphic")
dir.create("results/rds", recursive = TRUE, showWarnings = FALSE)
for (s in steps) { message("\n======== ", s, " ========"); t0 <- Sys.time()
  source(file.path("scripts", paste0(s, ".R")), local = new.env())
  message(sprintf("---- %s finished in %.0f s", s, as.numeric(Sys.time() - t0, units = "secs"))) }
writeLines(capture.output(sessionInfo()), "results/sessionInfo.txt")
message("\nAll steps finished. Tables: results/tables | Figures: results/figures")
