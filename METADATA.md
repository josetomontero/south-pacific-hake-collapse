# Metadata of the tables

Files with the extension `.csv` are comma separated. Files with the extension `.txt` in `data/` and
`results/` are tab separated, with a header row. Missing values are `NA`. Years are calendar years.

Notation of the models. `R` is the realized per capita rate of change of the hake biomass,
R_t = ln(B_t+1 / B_t). `Xt` is the natural logarithm of the hake biomass in t. The first number of
a model is the maximum per capita rate of change (Rm), fixed at the highest observed R_t of the
period. `U` is the exploitation rate, `SST` the smoothed sea surface temperature over the fishery
management unit and `ENSO` the smoothed sea surface temperature of the Niño 3.4 region. `Y` is the
natural logarithm of the jumbo squid biomass. The year of each variable is given in brackets
(`[t]`, `[t-1]`). Variables inside the parentheses that multiply `exp(Xt)` act on intraspecific
competition (lateral effect) and variables outside act on the reproductive rate (vertical effect).

## data/derived/hake_dataset.txt

Analysis dataset, one row per year (1968-2021).

| Column | Description | Unit |
|---|---|---|
| `year` | Year | |
| `catch` | Hake catch | t |
| `biomass` | Hake total biomass of the IFOP stock assessment | t |
| `logB` | Natural logarithm of `biomass` | |
| `R` | Reproductive rate, ln(B_t+1 / B_t) | per year |
| `U` | Exploitation rate, catch / biomass | |
| `squid_B` | Jumbo squid total biomass of the IFOP assessment (2001-2022) | thousand t |
| `log_squid_B` | Natural logarithm of `squid_B` | |
| `SST0`, `SST1`, `SST2` | Smoothed annual mean sea surface temperature over the fishery management unit, at lags of 0, 1 and 2 years | °C |
| `ENSO34`, `ENSO341`, `ENSO342` | Smoothed annual mean sea surface temperature of the Niño 3.4 region, at lags of 0, 1 and 2 years | °C |
| `ONI0`, `ONI1` | Smoothed annual mean Oceanic Niño Index, at lags of 0 and 1 year. Not used in the models | °C |

## results/tables/candidate_models_1968_2020.csv, candidate_models_1975_2000.csv, candidate_models_2001_2020.csv

All the candidate models of a period (243, 243 and 20 rows), one row per model, ordered by validity
and AICc.

| Column | Description |
|---|---|
| `rank` | Order of the model in the file |
| `Models` | Model formula |
| `structure_ID` | Number of the model structure (1 to 27, or 1 to 8 in 2001-2020) |
| `sst_lag`, `enso_lag` | Lag of SST and of ENSO in the model, in years. `NA` when the model does not include the variable |
| `N.param` | Number of estimated parameters |
| `c` | Intraspecific competition parameter |
| `b`, `e`, `d` | Coefficients of the predictors. In 1968-2020 and 1975-2000, `b` is the coefficient of U, `e` of SST and `d` of ENSO. In 2001-2020, `d` is the coefficient of a variable with a vertical effect, `e` of a variable with a lateral effect and `b` of the jumbo squid biomass |
| `AICc` | Akaike information criterion corrected for small samples |
| `R2` | Squared correlation between the observed and the fitted reproductive rates |
| `RMSD` | Root mean squared deviation between the observed and the simulated log biomass, every predicted year against the observation of the same year |
| `RMSD_steps_04_05` | RMSD with the expression of the submitted analysis (see README) |
| `min_k` | Lowest value of the competition term of the model over the years of the period |
| `valid` | `TRUE` when the competition term is positive in every year |
| `converged` | `TRUE` when the fit converged |
| `n_versions` | Number of valid lag versions of the structure of the model |
| `DAICc` | Difference in AICc from the valid model with the lowest AICc |
| `wi` | Akaike weight over the valid models, with prior weights equal among structures and, within a structure, equal among its valid lag versions |
| `wi_equal_priors` | Akaike weight over the valid models with the same prior weight for every model |
| `W_structure` | Weight of the structure, sum of `wi` over its lag versions |
| `W_structure_equal_priors` | Sum of `wi_equal_priors` over the lag versions of the structure |

## results/tables/model_structures_1968_2020.csv, model_structures_1975_2000.csv, model_structures_2001_2020.csv

One row per model structure, at its valid lag version with the lowest AICc. `ID` is the number of the
structure and `lag_versions` the number of valid lag versions out of the versions fitted. The other
columns are those of the `candidate_models` files.

## results/tables/lag_fits_1968_2020.csv, lag_fits_1975_2000.csv, lag_fits_2001_2020.csv, lag_fits_2001_2020_27_structures.csv

Every model structure fitted at every lag combination of SST and ENSO (lags 0 to 3), one row per
model. `lag_fits_2001_2020.csv` holds the eight structures of 2001-2020 and
`lag_fits_2001_2020_27_structures.csv` the 27 structures of the other periods fitted to 2001-2020,
which are not candidate models of that period. Columns `rank` to `AICc`, `R2`, `min_k`, `valid` and
`converged` are those of the `candidate_models` files.

| Column | Description |
|---|---|
| `DAICc`, `wi` | Difference in AICc and Akaike weight over all the models of the file |
| `DAICc_valid`, `wi_valid` | Difference in AICc and Akaike weight over the valid models, with the same prior weight for every model |

## results/tables/models_1968_2020.txt, models_1975_2000.txt, models_2001_2020.txt

Models of one lag combination (SST and ENSO at a lag of one year in 1968-2020 and without lag in the
other periods). `ID` is the number of the structure. `DAICc` and `wi` are computed within the file.
`RMSD` uses the expression of the submitted analysis (see README). The other columns are those of the
`candidate_models` files.

## Other tables of results/tables

| File | Columns |
|---|---|
| `candidate_models_selection.txt`, `lag_fits_summary.txt` | Text summaries of the selection and of the lag fits |
| `bootstrap_R_P1.txt`, `_P2.txt`, `_P3.txt` | `Year`, observed reproductive rate (`R_obs`), mean of the bootstrap predictions (`R_pred`), 2.5 and 97.5 percentiles (`ci_low`, `ci_high`) and number of bootstrap fits that converged (`n_ok`). P1 is 1968-2020, P2 1975-2000 and P3 2001-2020 |
| `breakpoints_hake_biomass.txt` | Number of breaks (`n_breaks`), residual sum of squares (`RSS`), `BIC` and years of the breaks (`break_years`) |
| `exploitation_rate_anomalies.txt` | `year`, exploitation rate (`U`) and its `anomaly`, the residual of the linear regression of U on year |
| `Fig_4_equilibria.txt` | Curves of Figure 4: `key` and `regime` of the curve, `model` number, `form` of the effect, `variable` and its `value`, variables held `fixed`, `Rm`, effective competition parameter (`c_eff`), equilibrium in log biomass (`Xstar`) and in thousand t (`Bstar_kt`), and the range of biomass observed in the regime (`B_min_kt`, `B_max_kt`) |
| `enso_episodes_1968_2020.txt` | One row per El Niño or La Niña episode: `phase`, year and season of onset and end, number of `seasons` and `peak` of the Oceanic Niño Index (°C) |
| `enso_episodes_by_period.txt` | Episodes by `period`: number of `years`, onsets of El Niño (`EN_onsets`) and La Niña (`LN_onsets`), onsets per decade, percentage of seasons within episodes, mean index (`mean_ONI`, °C) and number of episodes with a peak of at least 1.5 °C (`EN_peak_ge_1.5`) or at most −1.5 °C (`LN_peak_le_m1.5`) |

## results/models

`best_P1.rds`, `best_P2.rds` and `best_P3.rds` are R objects with the selected model of steps 04 and
05, its number, the simulation formulas and the data of the period. `sim_P1.txt`, `sim_P2.txt` and
`sim_P3.txt` hold the simulations of Figure 3, in thousand t: `years`, observed biomass (`Nt`),
one-step-ahead prediction (`Sim1`), complete-trajectory prediction (`Sim2`), their 95% confidence
intervals from parameter uncertainty (`.CI.L`, `.CI.U`) and the 95% prediction interval with process
error (`Sim2.PI.L`, `Sim2.PI.U`, median `Sim2.PI.median`).

## appendix_analyses

The tables of `appendix_analyses/A_assessment_cases/tables` and `B_acoustic_survey/tables` are
described in the README of each folder.
