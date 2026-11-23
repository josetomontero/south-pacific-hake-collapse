# Collapsed period 2001-2020: eight candidate models with at most two estimated
# parameters (c plus one of U, SST, ENSO in vertical or lateral position, or the
# predation parameter b of the prey model, equation 6). Predictors are not lagged.
# Jumbo squid biomass (Y_t) is available from 2001, which defines the period.
# Outputs: results/tables/models_2001_2020.txt, results/models/best_P3.rds, sim_P3.txt (Fig. 3c).
if (!exists("PROJ")) source(file.path(rprojroot::find_root(rprojroot::has_file("R/config.R")), "R", "config.R"))
suppressPackageStartupMessages({ library(minpack.lm); library(qpcR); library(Metrics) })

hd  <- tsv_read(file.path(PROJ, "data", "derived", "hake_dataset.txt"))
per <- PERIODS$P3
data3 <- data.frame(Year = hd$year, Yt = hd$log_squid_B, Xt = hd$logB, R = hd$R,
                    V1 = hd$U, V2 = hd$SST0, V3 = hd$ENSO34)
data3 <- na.omit(subset(data3, Year >= per$from & Year <= per$to))
rm <- round(max(data3$R), 2)
stopifnot(rm == per$rm)
cat(sprintf("%s: n = %d (R_t for %d-%d), Rm = %.2f\n", per$label, nrow(data3), min(data3$Year), max(data3$Year), rm))

# Model set in the order of Table S3 (d = vertical effect, e = lateral effect, b = predation)
cand <- list(
  list(f = R ~ rm - c * exp(Xt),                 s = c(c = 1),             sim = "RM-(c[j])*exp(N)"),
  list(f = R ~ rm - c * exp(Xt) + d * V1,        s = c(c = 1, d = 1),      sim = "RM-(c[j])*exp(N)+d[j]*V1[t]"),
  list(f = R ~ rm - c * exp(Xt) + d * V2,        s = c(c = 1, d = 1),      sim = "RM-(c[j])*exp(N)+d[j]*V2[t]"),
  list(f = R ~ rm - c * exp(Xt) + d * V3,        s = c(c = 1, d = 1),      sim = "RM-(c[j])*exp(N)+d[j]*V3[t]"),
  list(f = R ~ rm - (c + e * V1) * exp(Xt),      s = c(c = 1, e = 1),      sim = "RM-(c[j]+e[j]*V1[t])*exp(N)"),
  list(f = R ~ rm - (c + e * V2) * exp(Xt),      s = c(c = 1, e = 1),      sim = "RM-(c[j]+e[j]*V2[t])*exp(N)"),
  list(f = R ~ rm - (c + e * V3) * exp(Xt),      s = c(c = 1, e = 1),      sim = "RM-(c[j]+e[j]*V3[t])*exp(N)"),
  list(f = R ~ rm - b * exp(Yt) - c * exp(Xt),   s = c(c = 1, b = 1e-4),   sim = "RM-b[j]*exp(Y[t])-(c[j])*exp(N)")
)
labels <- c("U[t]", "SST[t]", "ENSO[t]", "Y[t]")
# Write the numeric Rm into the formulas so the saved model objects are self-contained
cand <- lapply(cand, function(k) { k$f <- as.formula(gsub("\\brm\\b", sprintf("%.2f", rm), deparse(k$f))); k })

models <- lapply(cand, function(k) nlsLM(k$f, data = data3, start = k$s,
                                         control = nls.lm.control(maxiter = 1000), algorithm = "BFGS"))
aic <- sapply(models, qpcR::AICc); w <- akaike.weights(aic)
pars <- plyr::rbind.fill(lapply(models, function(m) as.data.frame(t(coef(m)))))
tab <- data.frame(ID = seq_along(models),
                  Models = sapply(models, function(m) {
                    f <- deparse(formula(m)); f <- sub("rm", sprintf("%.2f", rm), f)
                    f <- gsub("V1", labels[1], gsub("V2", labels[2], gsub("V3", labels[3], f))); gsub("Yt", labels[4], f) }),
                  N.param = sapply(models, function(m) length(coef(m))),
                  pars[, c("c", "d", "e", "b")],
                  AICc = aic, DAICc = w$deltaAIC, wi = w$weights,
                  R2 = sapply(models, function(m) cor(data3$R, predict(m))^2))

# Complete-trajectory simulation (the prey model needs the squid biomass; the simulation
# strings reference Y, which is made available in the calling frame)
sim_strings <- function(k) c(gsub("N", "N1a", gsub("RM", sprintf("%.2f", rm), k$sim)),
                             gsub("N", "N1b", gsub("RM", sprintf("%.2f", rm), k$sim)))
Y <- data3$Yt
sims <- lapply(seq_along(models), function(i)
  simula3(data3$Xt[1], models[[i]], sim_strings(cand[[i]]), data3$Year, data3$Xt, data3$V1, data3$V2, data3$V3,
          iter = 200, log = "no", process_error = FALSE))
tab$RMSD <- sapply(sims, function(x) rmse(log(x$Nt[-length(x$Nt)]), log(x$Sim2[-length(data3$Xt)])))
tsv_write(tab, file.path(RES, "tables", "models_2001_2020.txt"))

best <- which.min(tab$AICc)
cat(sprintf("best model ID %d: %s  (wi = %.2f, R2 = %.2f, RMSD = %.3f)\n",
            best, tab$Models[best], tab$wi[best], tab$R2[best], tab$RMSD[best]))
cat(sprintf("evidence ratio best vs second: %.2f; predation model dAICc = %.2f, wi = %.3f\n",
            tab$wi[best] / sort(tab$wi, decreasing = TRUE)[2], tab$DAICc[8], tab$wi[8]))
saveRDS(list(period = "P3", ID = best, model = models[[best]], sim_formulas = sim_strings(cand[[best]]), data = data3),
        file.path(RES, "models", "best_P3.rds"))
simm <- simula3(data3$Xt[1], models[[best]], sim_strings(cand[[best]]), data3$Year, data3$Xt, data3$V1, data3$V2, data3$V3,
                iter = 5000, log = "no", process_error = TRUE)
tsv_write(simm, file.path(RES, "models", "sim_P3.txt"))
