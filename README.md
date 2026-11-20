# Random cross-validation overstates transferability of leaf chlorophyll estimates from in situ VIS-NIR reflectance in sugarcane

Code, data, tables and figures for the manuscript

> Pastrana Palacios, C. C.; Reyes-Trujillo, A.; Benavides-Bolaños, J. A. *Random Cross-Validation Overstates Transferability of
> Leaf Chlorophyll Estimates from In Situ VIS-NIR Reflectance in Sugarcane (Saccharum spp.)*. Universidad del Valle, Cali, Colombia.
> (submitted; citation and DOI to be added after release)

PLSR for leaf chlorophyll a and b from 477 sugarcane leaves (CC 01-1940, first ratoon, LASA station, 20 plots, six campaigns,
60-210 days after cutting) under random, leave-one-plot-out, leave-one-block-out and leave-one-campaign-out validation, with
campaign-mean benchmarks, EPO, campaign-specific and stage-transfer models, and sensitivity analyses.

## Quick start

```bash
# R >= 4.3 with: pls, randomForest, ggplot2, patchwork, readxl
Rscript tests/test_leakage.R            # ~5 s: validation schemes never share plots / blocks / campaigns
Rscript scripts/run_all.R               # ~35 min on one core (the Monte-Carlo step dominates)
Rscript tests/check_manuscript_numbers.R  # compares the tables written by run_all.R with the manuscript values
```

Run everything from the repository root (in RStudio: `setwd()` to the repository folder, then `source("scripts/run_all.R")`).
Scripts can also be run one by one. Outputs: `results/tables/*.csv`, `results/figures/*.png|pdf`, `results/rds/` (intermediate
objects, git-ignored), `results/sessionInfo.txt`. Every random step uses a fixed seed (`CFG$seed = 42`), so repeated runs give the same numbers.

## Layout

```
config.R                       paths, seeds and options (all scripts source it)
R/functions.R                  data loading, metrics, MSC / EPO, PLSR, validation schemes, cluster bootstrap
R/plot_style.R                 figure theme and export at print size
scripts/00-13_*.R              analyses; run_all.R runs them in order
matlab/preprocess_spectra.m    raw spectra -> 500 bands, Savitzky-Golay, MSC (MATLAB; not run here)
tests/                         validation-scheme invariants; manuscript-number check
data/                          leaf chemistry and spectra (see data/README.md)
REPRODUCIBILITY.md             manuscript item -> script -> output, and data-processing notes
```

## Manuscript item -> script

See `REPRODUCIBILITY.md` for the full table (every manuscript/SI table and figure -> script -> output file), the data-processing
notes and the few known small differences. `scripts/00-13_*.R` are run in order by `run_all.R`.

## Metrics

`r2` squared Pearson correlation between observed and predicted; `R2` predictive coefficient of determination `1 - SSE/SST`
(SST around the mean of the evaluated observations; negative when the error exceeds that reference); `RMSE` in mg g-1 fresh weight;
`RPD = SD(obs)/RMSE`; `RPIQ = IQR(obs)/RMSE`. See `R/functions.R::metrics`.

## Validation schemes

All schemes give out-of-fold predictions for the same 477 leaves: `random` (5-fold, stratified by campaign, ignores plots), `lopo`
(leave-one-plot-out), `lobo` (leave-one-block-out), `loco` (leave-one-campaign-out), `strict` (campaign and plot group held out
together). In every training set the MSC reference is re-estimated from the training spectra and the number of latent variables
(1-15) is chosen by an inner 5-fold CV grouped by plot. Uncertainty comes from resampling the 20 plots (cluster bootstrap, B = 1000).

## Data and licence

The code is released under the MIT licence (`LICENSE`). The data in `data/` are released under the licence stated in
`data/README.md`. Please cite the manuscript (and the archived release, DOI to be added) if you use them.
