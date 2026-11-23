# Residual bootstrap (5000 resamples) of the best model of each period to obtain the
# predicted reproductive rate R_t with 95% intervals (Figure 3d-f).
# Outputs: results/tables/bootstrap_R_P1.txt, _P2.txt, _P3.txt
if (!exists("PROJ")) source(file.path(rprojroot::find_root(rprojroot::has_file("R/config.R")), "R", "config.R"))
library(boot)

boot_R <- function(p, nboot = 5000) {
  b <- readRDS(file.path(RES, "models", paste0("best_", p, ".rds")))
  model <- b$model; data <- b$data
  data$rm <- PERIODS[[p]]$rm   # Rm appears as a fixed constant in the model formulas
  f <- formula(model); start <- coef(model)
  fitted_values <- fitted(model); res <- residuals(model)
  set.seed(123)
  bf <- function(dat, idx) {
    dat$R <- fitted_values + sample(res, length(res), replace = TRUE)
    refit <- try(nls(f, data = dat, start = start), silent = TRUE)
    if (inherits(refit, "try-error")) rep(NA, length(start)) else coef(refit)
  }
  br <- boot(data, bf, R = nboot)
  pred <- t(apply(br$t, 1, function(cf) {
    if (anyNA(cf)) return(rep(NA, nrow(data)))
    e <- as.list(data); e[names(start)] <- cf; e$rm <- PERIODS[[p]]$rm
    eval(f[[3]], e)
  }))
  out <- data.frame(Year = data$Year, R_obs = data$R,
                    R_pred = colMeans(pred, na.rm = TRUE),
                    ci_low = apply(pred, 2, quantile, 0.025, na.rm = TRUE),
                    ci_high = apply(pred, 2, quantile, 0.975, na.rm = TRUE),
                    n_ok = colSums(!is.na(pred)))
  cat(sprintf("%s: model ID %d, %d/%d bootstrap refits converged\n", p, b$ID, min(out$n_ok), nboot))
  tsv_write(out, file.path(RES, "tables", paste0("bootstrap_R_", p, ".txt")))
}
for (p in names(PERIODS)) boot_R(p)
