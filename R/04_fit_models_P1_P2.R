# Fit the 27 candidate Verhulst-Royama models (equations 3-4) for
#   P1: whole period 1968-2020, predictors U_t, SST_{t-1}, ENSO_{t-1}, Rm = 0.31
#   P2: low exploitation period 1975-2000, predictors U_t, SST_t, ENSO_t, Rm = 0.21
# Rm is fixed at the highest observed R_t of each period (rounded to 2 decimals).
# Outputs: results/tables/models_1968_2020.txt, models_1975_2000.txt (the 27 models, AICc, wi, R2, RMSD),
#          results/models/best_P1.rds, best_P2.rds, sim_P1.txt, sim_P2.txt (Fig. 3a,b).
if (!exists("PROJ")) source(file.path(rprojroot::find_root(rprojroot::has_file("R/config.R")), "R", "config.R"))
suppressPackageStartupMessages({ library(Metrics); library(doParallel); library(foreach); library(stringr) })

hd <- tsv_read(file.path(PROJ, "data", "derived", "hake_dataset.txt"))
registerDoParallel(max(1, parallel::detectCores() - 1))

fit_period <- function(p, id) {
  per <- PERIODS[[p]]
  data <- data.frame(Year = hd$year, Xt = hd$logB, R = hd$R, V1 = hd$U,
                     V2 = if (per$lag == 1) hd$SST1 else hd$SST0,
                     V3 = if (per$lag == 1) hd$ENSO341 else hd$ENSO34)
  data <- na.omit(subset(data, Year >= per$from & Year <= per$to))
  rm <- round(max(data$R), 2)
  stopifnot(rm == per$rm)
  cat(sprintf("\n%s: n = %d (R_t for %d-%d), Rm = %.2f\n", per$label, nrow(data), min(data$Year), max(data$Year), rm))

  set.seed(122)
  allmodels <- pop_models(predictors = 3, Rmax = "obs", Rm = rm, modelpars = "yes", model = "Verhulst", data = data)

  # Complete-trajectory simulation of every candidate (simula3, 200 parameter draws) and RMSD on log biomass
  sim <- foreach(i = seq_along(allmodels[[1]][[1]]), .export = c("modelSelector", "simula3")) %dopar% {
    m <- modelSelector(allmodels, obj = i)
    x <- simula3(data$Xt[1], m[[1]][[1]], m[[2]], data$Year, data$Xt, data$V1, data$V2, data$V3,
                 iter = 200, log = "no", process_error = FALSE)
    list(x, rmse(log(x$Nt[-length(x$Nt)]), log(x$Sim2[-length(data$Xt)])))
  }
  tab <- allmodels[[1]]
  tab$RMSD <- sapply(sim, `[[`, 2)
  lab <- if (per$lag == 1) c("SST[t-1]", "ENSO[t-1]") else c("SST[t]", "ENSO[t]")
  tab$Models <- gsub("V3", lab[2], gsub("V2", lab[1], gsub("V1", "U[t]", tab$Models)))
  tab <- cbind(ID = seq_len(nrow(tab)), tab)
  tsv_write(tab, file.path(RES, "tables", sprintf("models_%d_%d.txt", per$from, per$to)))

  best <- which.min(tab$AICc)
  cat(sprintf("best model ID %d: %s  (wi = %.2f, R2 = %.2f, RMSD = %.3f)\n",
              best, tab$Models[best], tab$wi[best], tab$R2[best], tab$RMSD[best]))
  m <- modelSelector(allmodels, obj = best)
  saveRDS(list(period = p, ID = best, model = m[[1]][[1]], sim_formulas = m[[2]], data = data),
          file.path(RES, "models", paste0("best_", p, ".rds")))

  # Final prediction of the best model (5000 draws) for Figure 3: parameter-uncertainty band
  # (Sim2.CI) and parameter + process-error prediction band (Sim2.PI)
  simm <- simula3(data$Xt[1], m[[1]][[1]], m[[2]], data$Year, data$Xt, data$V1, data$V2, data$V3,
                  iter = 5000, log = "no", process_error = TRUE)
  tsv_write(simm, file.path(RES, "models", paste0("sim_", p, ".txt")))
  invisible(tab)
}

t1 <- fit_period("P1", 1)
t2 <- fit_period("P2", 2)
