# B. Comparison with the acoustic-survey biomass

Question. Do the models of 2001-2020 agree with the biomass of the acoustic surveys, which does not
depend on the assessment model?

Data.

- `data/survey_biomass_table10_2024.csv`. Biomass of the winter acoustic surveys of 2002 and
  2004-2024, from Table 10 of the IFOP 2024 survey report (Molina et al. 2025). Columns `year` and
  `biomass_t` (t).
- `data/survey_biomass_fig3_IFOP2022_digitised.csv`. Survey biomass digitized from Figure 3 of the
  IFOP 2022 status report (Tascheri 2022), used for 1995, 1997 and 1999-2001. Columns `year`,
  `biomass_kt_digitised` (thousand t) and `biomass_kt_table10_2024` (the value of Table 10 for the
  same year, thousand t, for checking).
- The Python scripts of `python/` produced both files from the report PDFs, which are not distributed.

There were no surveys in 1996, 1998 and 2003. The value of 2002 was flagged as inconsistent in the
IFOP 2022 report and is not used in any fit.

Analyses (`R/B1_survey_tests.R`, run from this folder after the steps of the repository root).

1. The eight models of 2001-2020 fitted to the assessment biomass are run as deterministic complete
   trajectories from the biomass of 2001 and compared with the survey biomass of 2004-2020 (n = 17)
   on the log scale, after removing the mean log difference (a catchability q).
2. The eight model structures are fitted directly to the survey biomass, with R_t = ln(B_t+1 / B_t)
   for consecutive surveys of 2004-2020 and U = catch / survey biomass.

Outputs.

- `tables/survey_biomass_series.txt`. `year`, survey biomass (`survey_t`, t) and `source`.
- `tables/B1_validation.txt`. One row per model and one for the assessment series: number of years
  (`n`), correlation with the survey on log biomass (`r_logB`), mean log difference (`log_q`) and
  catchability (`q`), RMSD on log biomass after scaling by q (`RMSD_logB_after_q`), and AICc and
  difference in AICc of the model fitted to the assessment biomass.
- `tables/B2_models_survey.txt`. Models fitted to the survey biomass: `n`, `Rm`, parameters (`c`,
  `d`, `e`, `b`), `converged`, `AICc`, `R2`, Shapiro-Wilk P of the residuals (`shapiro_p`), residual
  autocorrelation at lag 1 (`acf1`), difference in AICc (`dAICc`) and Akaike weight (`wi`).
- `figs/B1_survey_vs_assessment.png`, `figs/B1_predicted_vs_survey.png`,
  `figs/B2_survey_models_dAICc.png`.
