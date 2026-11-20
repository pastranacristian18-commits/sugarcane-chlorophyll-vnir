# Reproducibility notes

All numbers in the manuscript and Supporting Information (SI) are produced by `scripts/run_all.R` (R 4.3.3, seed 42) from the files in `data/`.
`tests/check_manuscript_numbers.R` compares 156 values printed in the manuscript/SI with the tables written by the scripts
(tolerance 0.0025 unless stated; all 156 agree in the version of this release).

## Manuscript item -> script -> table

| Item | Script | Table written (`results/tables/`) |
|---|---|---|
| Table 1 | `01_descriptives_correlations.R` | `table01_descriptives_by_campaign.csv` |
| Table 2 | `02_random_split_plsr.R` | `table02_plsr_random_8020.csv` |
| Table 3, Tables S2-S3 | `03_validation_schemes.R` | `table03_*`, `table03b_*`, `table04_loco_per_campaign.csv` |
| 1000 repeated random 80/20 splits (text, Fig. S8) | `03b_monte_carlo_8020.R` | `table03d_monte_carlo_random_8020.csv` |
| Table 4 | `04_epo_nested.R` | `table05_epo_loco.csv` |
| Within-campaign skill (text, Fig. S5-S6), Table S6 | `05_within_campaign_skill.R` | `table06_*`, `table07_*`, `table08_loso_stages.csv` |
| Table S5 (ICC) | `06_diagnose_stage_models.R` | `table09_stage_model_diagnosis.csv` |
| Table S4 (180 DAC designs) | `13_design_dependence_180ddc.R` | `table18_design_dependence_180ddc.csv` |
| Table S1 (importance), PCA-RF r2 | `07_variable_importance.R` | `table10_*`, `table10b_*` |
| Table S8 (MSC reference) | `08_msc_sensitivity.R` | `table11_msc_leakage_sensitivity.csv` |
| Table S9 (leaf P6S1) | `12_p6s1_sensitivity.R` | `table17_p6s1_sensitivity.csv` |
| Table S10 (band exclusion) and P6S1 reflectance check | `10_band_exclusion_sensitivity.R` | `table12a_*`, `table12b_*` |
| Table S7 (FDR, TCA, vegetation indices) | `11_analisis_adicionales.R` | `table13_*`, `table15_*`, `table16_*` |
| Figures | `09_figures.R` | `results/figures/` (see below) |
| TOC graphic (3.25 x 1.75 in, submission file) | `14_toc_graphic.R` | `results/figures/toc_graphic.png`, `.pdf` |

## Data processing

* The data file has 479 rows; the second copy of the duplicated record Dia 6 / P20S4 is removed (478 leaves).
* Leaf P6S1 at 150 DAC is excluded (`drop_p6s1_150ddc = TRUE`, n = 477). Its Chl a (1.48 mg g-1 FW) differs from the other leaves of the plot
  (0.43-0.59) although its spectrum is normal; this does not establish analytical error. `12_p6s1_sensitivity.R` (Table S9) repeats all
  validation schemes with n = 478 (RPD changes by at most 0.02). Its minimum raw reflectance is 10.7 %, so the earlier "reflectance < 1 %" reason does not apply.
* Two leaves at 150 DAC (P1S3, P4S4) have no valid spectra.
* `Bloque_Num` was corrected against the field map (randomized complete block design, 5 blocks x 4 N rates, each treatment once per block).
  Chlorophyll values, spectra and treatments were not changed.
* Preprocessing: Savitzky-Golay (11 points, 2nd order) and MSC. `data/Base_de_datos_socaaaa.txt` holds the MSC-corrected matrix whose MSC
  reference is the mean of all 478 spectra (used for display and for most schemes); the cross-validation tables that state "reference estimated
  within each calibration set" use the SG-only matrix and re-estimate MSC in every training set (`pp_msc`). Table S8 shows that the two
  implementations differ by <= 0.006 % in RMSE.

## Known differences between this code and the printed manuscript

Values that agree to rounding are not listed. These differ slightly:

* Table S3: R2 and RMSE of single campaigns differ by up to 0.002 (e.g. Chl a 210 DAC R2 -3.332 here, -3.334 printed).
* Table S6: Chl a Elongation -> Tillering R2 = -4.768 here, -4.772 printed.
* Table S1: mean PCA-RF importance differs by up to 0.005 (random forest; PLSR VIP values agree).
* Table S7: TCA with m = 30 components gives r2 0.013 / R2 -0.311 / RMSE 0.5485 / RPD 0.874 (Chl a) and 0.009 / -0.385 / 0.2419 / 0.851 (Chl b) here,
  versus 0.010 / -0.335 / 0.5533 / 0.867 and 0.007 / -0.405 / 0.2436 / 0.845 printed. m = 10 and m = 20 agree. The conclusion (no useful skill) is unchanged.
  Vegetation-index rows differ by <= 0.002.
* Figures: `09_figures.R` draws the same quantities as the manuscript figures but with a different layout and numbering (below);
  they are not pixel-identical to the figures in the submitted files.

## Figure map (manuscript -> file written by `09_figures.R`)

| Manuscript | File |
|---|---|
| Fig. 1 (correlation profiles) | `fig02_pearson_profiles` |
| Fig. 2 (observed vs predicted, 80/20) | `fig04_obs_vs_pred_chla`, `fig05_obs_vs_pred_chlb` |
| Fig. 3 (RPD by scheme) | `fig07_validation_schemes` |
| Fig. 4 (importance) | `fig06_spectral_importance` |
| Fig. S2 (spectra) | `fig03_spectra_raw_vs_processed` |
| Fig. S3 (LOCO predictions) | `fig11_loco_predictions_by_campaign` |
| Fig. S4 (EPO) | `fig10_epo_correlation_vs_error` |
| Fig. S5 (180 DAC designs) | `fig08_stage_models_by_validation_design` |
| Fig. S6 (campaign-centred models) | `fig12_within_campaign_skill` |
| Fig. S7 (N-rate trajectories) | `fig09_chlorophyll_trajectory_by_N_dose` |
| Fig. S8 (repeated splits vs LOCO) | `fig07s_monte_carlo_vs_loco_per_fold` |
| Fig. S1 (site map) | not generated by code |
