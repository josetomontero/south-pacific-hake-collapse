# Complement to A1: does jumbo squid predation add anything once the
# exploitation rate is in the model? The paper's candidate set for 2001-2020 has at most
# two parameters per model (c plus one predictor), so U and squid never appear together. Here two
# three-parameter models are added to the eight of the paper and fitted to the paper series and to
# the four IFOP 2022 assessment cases (same data and settings as A1_fit_P3_cases.R):
#   9  U lateral + squid predation   R_t = Rm - (c + e U_t) exp(X_t) - b exp(Y_t)
#   10 U vertical + squid predation  R_t = Rm - c exp(X_t) - b exp(Y_t) + d U_t
# For each series: AICc, dAICc and Akaike weights over the ten models, the predation parameter b with
# its SE, and a likelihood-ratio test of b = 0 (model 9 vs 5; model 10 vs 2; anova on nested nls fits).
# Outputs: tables/A3_U_plus_squid_all_models.txt (ten models per series), tables/A3_U_plus_squid_summary.txt.
suppressPackageStartupMessages({ library(minpack.lm); library(qpcR) })
A <- rprojroot::find_root(rprojroot::has_file("R/A1_fit_P3_cases.R"))
HAKE <- normalizePath(file.path(A, "..", ".."))   # repository root 

hd <- read.delim(file.path(HAKE, "data", "derived", "hake_dataset.txt"))
cs <- read.csv(file.path(A, "data", "IFOP2022_cases0-3_1997_2021.csv"))
s3 <- read.delim(file.path(HAKE, "results", "tables", "models_2001_2020.txt"))

cand <- list(
  m1  = list(f = "R ~ RM - c * exp(Xt)",                              s = c(c = 1)),
  m2  = list(f = "R ~ RM - c * exp(Xt) + d * V1",                     s = c(c = 1, d = 1)),
  m3  = list(f = "R ~ RM - c * exp(Xt) + d * V2",                     s = c(c = 1, d = 1)),
  m4  = list(f = "R ~ RM - c * exp(Xt) + d * V3",                     s = c(c = 1, d = 1)),
  m5  = list(f = "R ~ RM - (c + e * V1) * exp(Xt)",                   s = c(c = 1, e = 1)),
  m6  = list(f = "R ~ RM - (c + e * V2) * exp(Xt)",                   s = c(c = 1, e = 1)),
  m7  = list(f = "R ~ RM - (c + e * V3) * exp(Xt)",                   s = c(c = 1, e = 1)),
  m8  = list(f = "R ~ RM - b * exp(Yt) - c * exp(Xt)",                s = c(c = 1, b = 1e-4)),
  m9  = list(f = "R ~ RM - (c + e * V1) * exp(Xt) - b * exp(Yt)",     s = c(c = 1, e = 1, b = 1e-4)),
  m10 = list(f = "R ~ RM - c * exp(Xt) - b * exp(Yt) + d * V1",       s = c(c = 1, d = 1, b = 1e-4)))
mlab <- c("1 basic", "2 U vertical", "3 SST vertical", "4 ENSO vertical", "5 U lateral", "6 SST lateral",
          "7 ENSO lateral", "8 squid predation", "9 U lateral + squid", "10 U vertical + squid")
`%||%` <- function(a, b) if (is.null(a)) b else a

fit_set <- function(B, years, U, label) {
  d <- data.frame(Year = years, Xt = log(B), R = c(log(B[-1] / B[-length(B)]), NA), V1 = U,
                  V2 = hd$SST0[match(years, hd$year)], V3 = hd$ENSO34[match(years, hd$year)],
                  Yt = hd$log_squid_B[match(years, hd$year)])
  d <- na.omit(subset(d, Year >= 2001 & Year <= 2020))
  rm <- round(max(d$R), 2)
  fits <- lapply(names(cand), function(k) {
    f <- as.formula(gsub("RM", sprintf("%.2f", rm), cand[[k]]$f))
    nlsLM(f, data = d, start = cand[[k]]$s, control = nls.lm.control(maxiter = 1000))
  })
  names(fits) <- names(cand)
  tab <- do.call(rbind, lapply(names(fits), function(k) {
    m <- fits[[k]]; cf <- as.list(coef(m)); se <- sqrt(diag(vcov(m)))
    data.frame(series = label, model = mlab[match(k, names(cand))], n = nrow(d), Rm = rm,
               c = cf$c, d = cf$d %||% NA, e = cf$e %||% NA, b = cf$b %||% NA, b_SE = if ("b" %in% names(se)) se[["b"]] else NA,
               converged = m$convInfo$isConv, AICc = qpcR::AICc(m), R2 = cor(d$R, fitted(m))^2,
               shapiro_p = shapiro.test(resid(m))$p.value, acf1 = acf(resid(m), plot = FALSE)$acf[2])
  }))
  tab$dAICc <- tab$AICc - min(tab$AICc)
  tab$wi <- exp(-tab$dAICc / 2) / sum(exp(-tab$dAICc / 2))
  lr9  <- anova(fits$m5, fits$m9); lr10 <- anova(fits$m2, fits$m10)
  summ <- data.frame(series = label, n = nrow(d), Rm = rm,
                     best = tab$model[which.min(tab$AICc)], w_best = tab$wi[which.min(tab$AICc)],
                     dAICc_5_U_lateral = tab$dAICc[5], dAICc_8_squid = tab$dAICc[8],
                     dAICc_9_Ulat_squid = tab$dAICc[9], w_9 = tab$wi[9], b_9 = tab$b[9], b_SE_9 = tab$b_SE[9],
                     LR_p_9_vs_5 = lr9$`Pr(>F)`[2],
                     dAICc_10_Uvert_squid = tab$dAICc[10], w_10 = tab$wi[10], b_10 = tab$b[10], b_SE_10 = tab$b_SE[10],
                     LR_p_10_vs_2 = lr10$`Pr(>F)`[2])
  list(tab = tab, summ = summ)
}

catch <- hd$catch[match(1997:2021, hd$year)]
runs <- list()
ok <- !is.na(hd$biomass)
runs[["Paper series"]] <- fit_set(hd$biomass[ok], hd$year[ok], hd$U[ok], "Paper series")
case_lab <- c("Case 0: official landings", "Case 1: official landings, Francis weighting",
              "Case 2: CCT catches", "Case 3: CCT catches, Francis weighting")
for (k in 0:3) {
  x <- subset(cs, case == k); x <- x[order(x$year), ]
  B <- x$total_biomass_kt * 1000
  runs[[case_lab[k + 1]]] <- fit_set(B, x$year, catch[match(x$year, 1997:2021)] / B, case_lab[k + 1])
}

# check: models 1-8 of the series of the paper reproduce models_2001_2020.txt
p <- runs[["Paper series"]]$tab[1:8, ]
chk <- max(abs(p$AICc - s3$AICc)); cat(sprintf("check vs models_2001_2020.txt (models 1-8): max |dAICc| = %.2e\n", chk))
stopifnot(chk < 1e-3)

all_tab <- do.call(rbind, lapply(runs, `[[`, "tab")); rownames(all_tab) <- NULL
rnd <- function(t) { for (v in c("AICc", "dAICc")) t[[v]] <- round(t[[v]], 2)
  for (v in c("wi", "R2", "shapiro_p")) t[[v]] <- round(t[[v]], 3); t$acf1 <- round(t$acf1, 2); t }
write.table(rnd(all_tab), file.path(A, "tables", "A3_U_plus_squid_all_models.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
summ <- do.call(rbind, lapply(runs, `[[`, "summ")); rownames(summ) <- NULL
num <- sapply(summ, is.numeric); summ[num] <- lapply(summ[num], signif, 3)
write.table(summ, file.path(A, "tables", "A3_U_plus_squid_summary.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
print(summ[, c("series", "best", "w_best", "dAICc_8_squid", "dAICc_9_Ulat_squid", "w_9", "b_9", "b_SE_9", "LR_p_9_vs_5",
               "dAICc_10_Uvert_squid", "w_10", "LR_p_10_vs_2")], row.names = FALSE)
