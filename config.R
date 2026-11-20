# ==============================================================================
# config.R -- central configuration (all scripts source this file)
# Paths are relative to the repository root. Run scripts from the root:
#   Rscript scripts/run_all.R
# ==============================================================================
CFG <- list(
  # ---- input data ----------------------------------------------------------
  data_txt    = "data/Base_de_datos_socaaaa.txt",                       # chemistry + SG+MSC spectra
  data_raw    = "data/matlab_export/DATOS_SOCA_MATRIZ_CRUDA.xlsx",      # raw spectra (MATLAB export)
  data_sgonly = "data/matlab_export/DATOS_SOCA_SG_SIN_MSC.xlsx",        # SG only, no MSC (MATLAB export)

  # ---- output --------------------------------------------------------------
  out_tables  = "results/tables",
  out_figures = "results/figures",
  out_rds     = "results/rds",

  # ---- data cleaning -------------------------------------------------------
  # The raw file contains one duplicated key (Day 6 / P20S4); the second copy is
  # always dropped (this is what the original code did with row 455).
  # The original code ALSO dropped row 468 = Day 4 (150 DDC) / P6S1, a valid,
  # non-outlying observation (z = 2.19). Set FALSE to keep all 478 leaves.
  # TRUE reproduces n = 477 of the manuscript. THE REASON MUST BE STATED IN THE PAPER.
  drop_p6s1_150ddc = TRUE,

  # ---- modelling -----------------------------------------------------------
  seed      = 42,
  max_lv    = 15,     # max latent variables explored
  inner_k   = 5,      # inner CV folds for LV selection (grouped by plot)
  k_random  = 5,      # folds of the "random 80/20" pooled scheme
  B_mc      = 1000,   # Monte-Carlo random 80/20 repeats (as in the manuscript)
  B_boot    = 1000,   # cluster-bootstrap replicates (resampling PLOTS, no refits)
  epo_kmax  = 4,      # fixed-k EPO curve (k = 1..4; max rank with 5 training campaigns)

  # ---- figure style (journal print sizes, inches) --------------------------
  fig_width_full = 7.2,    # ~ two-column width
  fig_width_half = 3.5,    # ~ one-column width
  dpi            = 600
)

# Crear las carpetas de salida automaticamente (si no existen), sin importar
# que script se corra primero o en que orden -- evita errores de "No such
# file or directory" al escribir tablas, figuras u objetos intermedios.
invisible(lapply(c(CFG$out_tables, CFG$out_figures, CFG$out_rds),
                  dir.create, recursive = TRUE, showWarnings = FALSE))
