# Project configuration. Every script sources this file first.
PROJ <- rprojroot::find_root(rprojroot::has_file("R/config.R"))

# HadISST netCDF (~485 MB) is not distributed with the repository.
# Download from https://www.metoffice.gov.uk/hadobs/hadisst/ and place it in
# data/external/, or point the environment variable HADISST_NC to it.
HADISST_NC <- Sys.getenv("HADISST_NC",
  unset = file.path(PROJ, "data", "external", "HadISST_sst.nc"))

# Environmental predictors are annual means: ENSO = annual mean of monthly Nino 3.4 SST
# (HadISST1); SST = annual mean of the HadISST1 cells inside the hake fishery management unit
# polygon (00_study_area_polygons.R), both smoothed (01_smooth_climate_indices.R). All outputs
# are written to results/.
RES <- file.path(PROJ, "results")

# Output folders. lapply (not a for loop) so no stray variable is left in the
# global environment: the model scripts use parameter names such as b, c, d, e.
invisible(lapply(c(file.path(PROJ, "data", "derived"), file.path(RES, c("tables", "models", "figs"))),
  function(p) dir.create(p, recursive = TRUE, showWarnings = FALSE)))

# Shared helpers
source(file.path(PROJ, "R", "lib", "auto_R_Verhulst_8.R"))   # pop_models(): builds and fits the 27 candidate models
source(file.path(PROJ, "R", "lib", "simula_v3.R"))           # simula3(): complete-trajectory simulation
source(file.path(PROJ, "R", "lib", "model_selector.R"))      # modelSelector(): picks one model + its simulation strings

# Breaks in the hake biomass series (strucchange, see 03_breakpoints.R)
BREAKS <- c(1974, 1983, 1992, 2001)

# Study periods (years of the reproductive rate R_t = log(B_{t+1}/B_t))
PERIODS <- list(
  P1 = list(label = "Whole period 1968-2020",             from = 1968, to = 2020, rm = 0.31, lag = 1),
  P2 = list(label = "Low exploitation period 1975-2000",  from = 1975, to = 2000, rm = 0.21, lag = 0),
  P3 = list(label = "Collapsed period 2001-2020",         from = 2001, to = 2020, rm = 0.31, lag = 0)
)

tsv_write <- function(x, path, rownames = FALSE) {
  write.table(x, path, sep = "\t", row.names = rownames, quote = FALSE)
  cat("wrote", sub(paste0(PROJ, "/"), "", path), "\n")
}
tsv_read <- function(path) read.delim(path, check.names = FALSE)
