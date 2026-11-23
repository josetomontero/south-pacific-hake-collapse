# Understanding the Collapse of the South Pacific Hake (*Merluccius gayi gayi*)

Data and code of the manuscript EAP25-0808 (*Ecological Applications*): Figures 1 to 4, Table 1,
all the candidate models of each period and their selection, and the fits of the candidate models to
the assessment cases and to the acoustic-survey biomass. Population-dynamics models
(Verhulst-Royama reproductive curve with lateral and vertical exogenous effects, and a prey model for
jumbo squid predation) are fitted to the IFOP stock-assessment biomass of South Pacific hake,
1968-2020.

Repository: <https://github.com/josetomontero/south-pacific-hake-collapse>

Archive: Zenodo, <https://doi.org/10.5281/zenodo.24471896> (all versions). Version 1.0.0,
<https://doi.org/10.5281/zenodo.24471897>.

License: MIT for the code and CC BY 4.0 for the data and outputs produced by the authors. The
illustration of the hake is copyrighted (see "License and copyright" below).

## How to run

```sh
# from the repository root
Rscript run_all.R
PYTHON=/path/to/python Rscript run_all.R   # to choose the Python interpreter of Figure 1
```

`run_all.R` runs the R steps in one R session, then the analyses of `appendix_analyses/` (each in
its own R process, from its folder), then the Python map. All outputs are written to `results/`,
`data/derived/` and `appendix_analyses/*/tables` and `figs`. The complete run takes about 40 s. Each
R script can also be run alone from the repository root (`Rscript R/<script>.R`) once the earlier
steps have been run. The scripts of `appendix_analyses/` are run from their own folder (for example
`cd appendix_analyses/A_assessment_cases && Rscript R/A1_fit_P3_cases.R`) and need the folders
`tables/` and `figs/` to exist.

Figures are written as PNG at 600 dpi. The figure number is not printed inside the files.

## Software

- R 4.5.1. The packages loaded by the steps and their versions are written to
  `results/sessionInfo.txt`.
- Optional step `R/00_extract_SST_HadISST.R` also needs the packages raster, xts, caTools, fields and
  ncdf4.
- Python 3.9 or later with matplotlib, numpy, pyproj and shapely 2 (Figure 1 only).

AICc in the model tables is `qpcR::AICc`, where k is the number of estimated coefficients.

Random numbers. Every step that draws random numbers sets its seed in the script
(`set.seed(122)` before the candidate models, `set.seed(123)` in the bootstrap, seed 6 inside
`simula3()`).

## Steps and outputs

| Step | Script | Output | Item |
|---|---|---|---|
| 00 | `R/00_study_area_polygons.R` | `data/derived/fishery_unit_polygon.gpkg`, `.geojson`, `LME_Humboldt.geojson` | Polygons of Figure 1 and of the SST extraction |
| 00 (optional) | `R/00_extract_SST_HadISST.R` | `data/derived/SSTchileStocks_2022.txt`, `SSTchileStocks_raw_2022.txt` | SST predictor (column `CHU`). Skipped when the HadISST file is absent |
| 01 | `R/01_smooth_climate_indices.R` | `data/derived/enso34_smoothed.txt`, `oni_smoothed.txt`, `results/figs/climate_smoothing.png` | ENSO predictor, figure of the raw and smoothed series |
| 02 | `R/02_build_dataset.R` | `data/derived/hake_dataset.txt` | Analysis dataset |
| 03 | `R/03_breakpoints.R` | `results/tables/breakpoints_hake_biomass.txt`, `results/figs/breakpoints.png` | Breakpoints 1974, 1983, 1992, 2001 |
| 04 | `R/04_fit_models_P1_P2.R` | `results/tables/models_1968_2020.txt`, `models_1975_2000.txt`, `results/models/best_P1.rds`, `best_P2.rds`, `sim_P1.txt`, `sim_P2.txt` | The 27 models of 1968-2020 and of 1975-2000 at the lags of the selected models, Table 1 |
| 05 | `R/05_fit_models_P3.R` | `results/tables/models_2001_2020.txt`, `results/models/best_P3.rds`, `sim_P3.txt` | The eight models of 2001-2020, Table 1 |
| 06 | `R/06_bootstrap_R.R` | `results/tables/bootstrap_R_P1.txt`, `_P2.txt`, `_P3.txt` | Residual bootstrap of the reproductive rate (Figure 3d-f) |
| 07 | `R/07_Fig2.R` | `results/figs/Fig_2.png`, `results/tables/exploitation_rate_anomalies.txt` | Figure 2 |
| 08 | `R/08_Fig3.R` | `results/figs/Fig_3.png` | Figure 3 |
| 09 | `R/09_Fig4.R` | `results/figs/Fig_4.png`, `results/tables/Fig_4_equilibria.txt` | Figure 4, drawn from the models of steps 04 and 05 |
| 10 | `R/10_lag_fits.R` | `results/tables/lag_fits_1968_2020.csv`, `lag_fits_1975_2000.csv`, `lag_fits_2001_2020.csv`, `lag_fits_2001_2020_27_structures.csv`, `lag_fits_summary.txt` | Lag tables: every model structure fitted at every lag combination of SST and ENSO (lags 0 to 3), from which the best models are selected |
| 20 | `R/20_candidate_set_selection.R` | `results/tables/candidate_models_1968_2020.csv`, `candidate_models_1975_2000.csv`, `candidate_models_2001_2020.csv`, `model_structures_1968_2020.csv`, `model_structures_1975_2000.csv`, `model_structures_2001_2020.csv`, `candidate_models_selection.txt` | Model selection over the full candidate set (243, 243 and 20 models): weights of Table 1 and the model tables and lag tables of the appendix |
| 19 | `R/19_ENSO_episodes.R` | `results/figs/ENSO_episodes.png`, `results/tables/enso_episodes_1968_2020.txt`, `enso_episodes_by_period.txt` | El Niño and La Niña episodes, 1968-2020 |
| A1 | `appendix_analyses/A_assessment_cases/R/A1_fit_P3_cases.R` | `tables/A1_P3_summary.txt`, `A1_P3_models_by_case.txt`, `A1_runs.rds` | Models of 2001-2020 fitted to the four assessment cases |
| A2 | `appendix_analyses/A_assessment_cases/R/A2_figures.R` | `figs/A2_cases_biomass.png`, `A2_cases_dAICc.png` | Figures of the assessment cases |
| A3 | `appendix_analyses/A_assessment_cases/R/A3_U_plus_squid.R` | `tables/A3_U_plus_squid_summary.txt`, `A3_U_plus_squid_all_models.txt` | Models with exploitation rate and squid predation together |
| B1 | `appendix_analyses/B_acoustic_survey/R/B1_survey_tests.R` | `tables/B1_validation.txt`, `B2_models_survey.txt`, `survey_biomass_series.txt`, `B_models.rds`, `figs/B1_*.png`, `B2_*.png` | Comparison with the acoustic-survey biomass |
| Python | `python/fig1_map.py` | `results/figs/Fig_1.png` | Figure 1 |

`run_all.R` also writes the R version and the packages of the run to `results/sessionInfo.txt`.

The columns of every table are described in `METADATA.md`.

Model tables. `candidate_models_*.csv` hold all the candidate models of a period, with their AICc,
Akaike weights and RMSD. `model_structures_*.csv` hold one row per model structure, at its valid lag
version with the lowest AICc, and are the model tables of the appendix of the manuscript. The model
selected in each period is in `candidate_models_selection.txt`. `models_*.txt` (steps 04 and 05) hold
the models of one lag combination, with weights computed within that combination. Steps 06, 08 and 09
read the fitted models of those steps.

RMSD. The column `RMSD` of the `candidate_models` and `model_structures` files compares every
predicted year of the complete trajectory with the observed log biomass of the same year, and these
are the values reported in the manuscript. The column `RMSD` of `models_*.txt` uses the expression of
the submitted analysis, which counted the first year (the initial value given to the model) and
paired the last observed year with the prediction of the following year. That expression is kept in
the column `RMSD_steps_04_05` of the other files.

Figure 4. The H0 curve of each regime is model 1 of `models_1975_2000.txt` and
`models_2001_2020.txt`. The dashed curves are the best models (`best_P2.rds`, `best_P3.rds`) with
ENSO at the upper quartile of 1975-2000 and SST at its mean, and with the exploitation rate at its
maximum. Every equilibrium of the figure lies within the biomass observed in its regime
(`Fig_4_equilibria.txt`).

Each folder of `appendix_analyses/` has its own README with the description of the analysis. The
stock-recruitment models of the appendix of the manuscript are a complementary analysis and are not
part of this repository.

Shared code. `R/config.R` (project root, periods, break years), `R/lib/auto_R_Verhulst_8.R`
(`pop_models()`, builds and fits the candidate models listed in `data/model_templates/`),
`R/lib/simula_v3.R` (`simula3()`, complete-trajectory simulation), `R/lib/model_selector.R`,
`R/lib/fig_theme.R`.

## Outputs that are not regenerated

- `data/derived/SSTchileStocks_2022.txt` and `SSTchileStocks_raw_2022.txt` are distributed. They are
  regenerated only when the HadISST1 file `HadISST_sst.nc` (about 485 MB, not distributed) is placed
  in `data/external/` or the environment variable `HADISST_NC` points to it. Source
  <https://www.metoffice.gov.uk/hadobs/hadisst/>.
- The CSV files extracted or digitized from IFOP reports are distributed
  (`appendix_analyses/*/data/*.csv`). The Python scripts that produced them
  (`appendix_analyses/*/python/`) need the report PDFs, which are not distributed, and poppler
  (`pdfimages`, `pdftotext`).

## Data

`data/fisheries/`
- `ch_ts_sst.txt`. Yearly hake catch (SERNAPESCA), hake total biomass (IFOP statistical
  catch-at-age assessment), exploitation rate, reproductive rate and jumbo squid biomass (IFOP
  assessment, 2001-2022).
- `chile_stocks_admb_outs_AES2.txt`. IFOP assessment outputs of South Pacific hake, 1960-2020 (total
  biomass, spawning biomass and fishing mortality), used by step 03.

`data/climate/`
- `enso34.txt`. Monthly Niño 3.4 SST from HadISST1 (1870-2021), with seasonal and annual means.
- `ONI.txt`. Monthly Oceanic Niño Index (NOAA CPC) and annual mean.
- `ONI_NOAA_CPC_seasonal_20261111.txt`. Seasonal Oceanic Niño Index (NOAA CPC, ERSSTv5), downloaded
  on 2026-11-11 from <https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt>. Used by step 19.

`data/external/`
- `LME_Humboldt.gpkg`. Humboldt Current Large Marine Ecosystem polygon (NOAA LME program).
- `natural_earth_paises_10m.geojson`. Natural Earth 1:10m countries.

`data/model_templates/`. Predictor combinations expanded by `pop_models()` into the 27 model
structures.

`appendix_analyses/A_assessment_cases/data/IFOP2022_cases0-3_1997_2021.csv`. Tables 2-5 of the IFOP
2022 hake status report (Tascheri 2022). `appendix_analyses/B_acoustic_survey/data/`. Acoustic-survey
biomass from Table 10 of the IFOP 2024 survey report (Molina et al. 2025) and values digitized from
Figure 3 of Tascheri (2022).

`figures/hake_illustration_Salinas_Nass.png`. Illustration by Pablo Salinas Nass, used in Figure 1.
© Pablo Salinas Nass, reproduced with permission (see `figures/COPYRIGHT.md`).

## License and copyright

- Code (`R/`, `python/`, `appendix_analyses/*/R/`, `appendix_analyses/*/python/`, `run_all.R`). MIT
  License, see `LICENSE`.
- Data and outputs produced by the authors (`data/derived/`, `results/`, the tables and figures of
  `appendix_analyses/`, and the files digitized or extracted by the authors). Creative Commons
  Attribution 4.0 International (CC BY 4.0), <https://creativecommons.org/licenses/by/4.0/>.
- Source data of other institutions. The catch series of SERNAPESCA, the stock-assessment and
  survey series of IFOP, the HadISST1 sea surface temperature of the Met Office Hadley Centre, the
  Oceanic Niño Index of the NOAA Climate Prediction Center, the Large Marine Ecosystem boundary of
  NOAA and the Natural Earth coastline are public data of those institutions. They are included to
  reproduce the analysis and remain under the terms of their providers, who must be cited as their
  source.
- Illustration of the hake (`figures/hake_illustration_Salinas_Nass.png`, also part of
  `results/figs/Fig_1.png`). © Pablo Salinas Nass, all rights reserved. Reproduced with the written
  permission of its author. It is not covered by the licenses above and may not be reused without
  the permission of its author (`figures/COPYRIGHT.md`).
