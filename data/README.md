# Data

| File | Content |
|---|---|
| `Base_de_datos_socaaaa.txt` | One row per leaf (tab-separated): `Unidad` (P = plot 1-20, S = leaf 1-4), `Bloque_Num`, `Tratamiento` (T1-T4 = 0/80/160/240 kg N ha-1), `Chla`, `Chlb` (mg g-1 fresh weight), `Dia` (Dia 1-6 = 60-210 days after cutting) and 500 reflectance columns (450-820 nm), already Savitzky-Golay smoothed and MSC-corrected. 479 rows, 478 unique leaves. |
| `matlab_export/DATOS_SOCA_MATRIZ_CRUDA.xlsx` | Raw reflectance resampled to the 500 bands (478 leaves). |
| `matlab_export/DATOS_SOCA_SG_SIN_MSC.xlsx` | Savitzky-Golay only (no MSC): input for the MSC reference estimated inside each calibration set. |

Cleaning applied by `R/functions.R::load_soca()` (details in `REPRODUCIBILITY.md`):
the second copy of the duplicated record Dia 6 / P20S4 is removed (479 -> 478 leaves), and leaf P6S1 at 150 days after cutting is
excluded by default (`drop_p6s1_150ddc = TRUE` in `config.R`; 478 -> 477 leaves; Table S9 gives the sensitivity). Two leaves at
150 days (P1S3, P4S4) have no valid spectra and are absent from the file.

`Bloque_Num` was corrected against the field map (randomized complete block design, 5 blocks x 4 N rates). Chlorophyll values,
spectra and `Tratamiento` were not altered.

## Licence

Data licence: **to be confirmed by the authors and the institution before this repository is made public** (CC BY 4.0 is the usual
choice). The MIT licence of the repository covers the code only.
