# A. Sensitivity of the 2001-2020 results to the assessment case

Question. Are the results of 2001-2020 the same when the hake biomass comes from assessment runs that
attribute more or less natural mortality to jumbo squid?

Data. `data/IFOP2022_cases0-3_1997_2021.csv` holds Tables 2-5 of the IFOP 2022 hake status report
(Tascheri 2022) for the four base cases of the assessment, 1997-2021 (100 rows, 25 per case).

| Column | Description | Unit |
|---|---|---|
| `case` | Case 0, official landings. Case 1, official landings with Francis (2011) weighting of the age compositions. Case 2, catch series of the Scientific Technical Committee (corrected for discards and artisanal under-reporting). Case 3, that catch series with Francis weighting | |
| `year` | Year | |
| `total_biomass_kt` | Total biomass (age 2+) | thousand t |
| `ssb_kt` | Spawning biomass | thousand t |
| `exploitable_biomass_kt` | Exploitable biomass | thousand t |
| `recruits_1e6` | Recruits | millions |
| `F` | Fishing mortality | per year |
| `u` | Exploitation rate printed in the case table | |

`python/extract_case_tables.py` extracted the file from the report PDF, which is not distributed.

Models. The eight models of 2001-2020 (`R/05_fit_models_P3.R` of the repository), with the same
formulas, starting values, fitting routine and rule for Rm, fitted to the total biomass of each case.
The exploitation rate is the catch of the series of the paper divided by the biomass of the case. As
a sensitivity, the exploitation rate printed in each case table is used instead. SST, ENSO and squid
biomass are those of `data/derived/hake_dataset.txt`. Fitting the series of the paper reproduces
`results/tables/models_2001_2020.txt`, and the script stops otherwise.

Run, from this folder, after the steps of the repository root.

```sh
Rscript R/A1_fit_P3_cases.R    # eight models per series
Rscript R/A2_figures.R         # figures
Rscript R/A3_U_plus_squid.R    # models with exploitation rate and squid predation together
```

Outputs.

- `tables/A1_P3_models_by_case.txt`. One row per series and model: number of years (`n`), `Rm`,
  parameters (`c`, `d`, `e`, `b`), `converged`, `AICc`, `R2`, RMSD of the simulated log biomass
  (`RMSD_logB`), Shapiro-Wilk P of the residuals (`shapiro_p`), residual autocorrelation at lag 1
  (`acf1`), difference in AICc (`dAICc`) and Akaike weight (`wi`) within the series.
- `tables/A1_P3_summary.txt`. One row per series: `Rm`, best model (`best`), its weight (`w_best`),
  R² (`R2_best`) and added parameter (`best_param`), difference in AICc, weight and parameter of the
  squid predation model (`dAICc_squid`, `w_squid`, `b_squid`) and difference in AICc of the model with
  the exploitation rate on intraspecific competition (`dAICc_U_lateral`).
- `tables/A3_U_plus_squid_all_models.txt`. The eight models plus model 9 (U lateral + squid predation)
  and model 10 (U vertical + squid predation), with the columns of `A1_P3_models_by_case.txt` and the
  standard error of the predation parameter (`b_SE`).
- `tables/A3_U_plus_squid_summary.txt`. One row per series: best of the ten models, differences in
  AICc of models 5, 8, 9 and 10, weights of models 9 and 10, predation parameter and its standard
  error in models 9 and 10, and P of the likelihood-ratio test of b = 0 (`LR_p_9_vs_5`,
  `LR_p_10_vs_2`).
- `figs/A2_cases_biomass.png`, `figs/A2_cases_dAICc.png`.
