# Run the complete analysis. From the repository root:  Rscript run_all.R
# Step 00_extract_SST_HadISST needs the HadISST netCDF (see R/config.R); it is skipped when the file is
# absent, because its outputs (data/derived/SSTchileStocks_2022.txt and SSTchileStocks_raw_2022.txt) are
# distributed with the repository.
source(file.path(rprojroot::find_root(rprojroot::has_file("R/config.R")), "R", "config.R"))
steps <- c("00_study_area_polygons", "01_smooth_climate_indices", "02_build_dataset", "03_breakpoints",
           "04_fit_models_P1_P2", "05_fit_models_P3", "06_bootstrap_R",
           "07_Fig2", "08_Fig3", "09_Fig4", "10_lag_fits", "20_candidate_set_selection",
           "19_ENSO_episodes")
# The SST extraction reads the management unit polygon, so it runs after 00_study_area_polygons.
if (file.exists(HADISST_NC)) steps <- append(steps, "00_extract_SST_HadISST", after = 1) else
  message("HadISST_sst.nc not found: skipping 00_extract_SST_HadISST.R (using distributed SST series)")
for (s in steps) {
  cat("\n==========", s, "==========\n")
  .step <- s; .t0 <- proc.time()[["elapsed"]]   # the scripts run in this session and define variables such as s
  source(file.path(PROJ, "R", paste0(s, ".R")), echo = FALSE)
  cat(sprintf("[%s: %.1f s]\n", .step, proc.time()[["elapsed"]] - .t0))
}

# Candidate models fitted to other biomass series (appendix_analyses/, Appendix S1: Sections S3 and S4). Each script
# finds its own folder from the working directory and defines its own PROJ and RES, so it is run in a separate R
# process from its folder. They read data/derived/hake_dataset.txt and results/tables/models_2001_2020.txt written above.
rscript <- file.path(R.home("bin"), "Rscript")
run_appendix <- function(folder, script) {
  cat("\n========== appendix_analyses/", folder, "/R/", script, " ==========\n", sep = "")
  dir <- file.path(PROJ, "appendix_analyses", folder)
  for (p in c("tables", "figs")) dir.create(file.path(dir, p), showWarnings = FALSE)
  old <- setwd(dir); on.exit(setwd(old))
  t0 <- proc.time()[["elapsed"]]
  st <- system2(rscript, shQuote(file.path("R", script)))
  if (st != 0) stop("appendix step failed (exit ", st, "): ", folder, "/R/", script)
  cat(sprintf("[%s: %.1f s]\n", script, proc.time()[["elapsed"]] - t0))
}
run_appendix("A_assessment_cases", "A1_fit_P3_cases.R")
run_appendix("A_assessment_cases", "A2_figures.R")
run_appendix("A_assessment_cases", "A3_U_plus_squid.R")
run_appendix("B_acoustic_survey", "B1_survey_tests.R")

# Figure 1 (map) is drawn in Python from the polygons written by 00_study_area_polygons.
# Needs a python3 (>= 3.9) with matplotlib, numpy, pyproj and shapely. Set the environment variable PYTHON
# to choose the interpreter.
py <- Sys.getenv("PYTHON", "python3")
run_py <- function(args) {
  cat("\n==========", paste(args, collapse = " "), "==========\n")
  st <- system2(py, shQuote(args))   # system2 quotes the command but not the arguments; paths may contain spaces
  if (st != 0) message("Python step failed (exit ", st, "); run it by hand: ", py, " ", paste(args, collapse = " "))
}
run_py(c(file.path(PROJ, "python", "fig1_map.py"), "--out", file.path(RES, "figs")))

# R version and packages used in this run
writeLines(capture.output(sessionInfo()), file.path(RES, "sessionInfo.txt"))
cat("sessionInfo written to results/sessionInfo.txt\n")
