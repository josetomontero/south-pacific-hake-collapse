# Are the 2001-2020 results the same when the hake biomass comes from
# assessment runs that attribute more or less natural mortality to jumbo squid?
#
# Biomass: total biomass 2+ of Cases 0-3 of the IFOP 2022 assessment (Tables 2-5, 1997-2021;
# data/IFOP2022_cases0-3_1997_2021.csv, extracted by python/extract_case_tables.py). All cases include
# the squid-proportional M from 1999; its estimate is largest in Cases 0-1 (official landings) and
# "not significant" in Cases 2-3 (CCT catches with discards and under-reporting) (report printed pp. 49 and 59).
# Reference: the series of the paper (data/derived/hake_dataset.txt of the repository, close to Case 2).
#
# Models: the 8 candidate models of the paper for 2001-2020 (R/05_fit_models_P3.R of the repository), same
# formulas, starting values, nlsLM call, Rm rule (Rm = round(max R_t, 2)) and predictors:
#   U   = catch of the paper series (hake_dataset; 1.6-1.9 x SERNAPESCA official landings in 2015-2020,
#         consistent with CCT-corrected catches) / total biomass of the case (paper definition);
#         sensitivity: the exploitation rate u printed in each case table (case-specific catches)
#   SST, ENSO = smoothed annual means (hake_dataset); Y = log jumbo squid biomass (hake_dataset)
# R_t = log(B_t+1 / B_t), t = 2001-2020.
# Outputs: tables/A1_P3_models_by_case.txt (all models), tables/A1_P3_summary.txt (one row per case),
# figures in A2_figures.R.
suppressPackageStartupMessages({ library(minpack.lm); library(qpcR) })
A <- rprojroot::find_root(rprojroot::has_file("R/A1_fit_P3_cases.R"))
HAKE <- normalizePath(file.path(A, "..", ".."))   # repository root 

hd <- read.delim(file.path(HAKE, "data", "derived", "hake_dataset.txt"))
cs <- read.csv(file.path(A, "data", "IFOP2022_cases0-3_1997_2021.csv"))

# Check with the series of the paper: reproduces results/tables/models_2001_2020.txt (stops otherwise)
s3 <- read.delim(file.path(HAKE, "results", "tables", "models_2001_2020.txt"))

cand <- list(
  m1 = list(f = "R ~ RM - c * exp(Xt)",               s = c(c = 1)),
  m2 = list(f = "R ~ RM - c * exp(Xt) + d * V1",      s = c(c = 1, d = 1)),
  m3 = list(f = "R ~ RM - c * exp(Xt) + d * V2",      s = c(c = 1, d = 1)),
  m4 = list(f = "R ~ RM - c * exp(Xt) + d * V3",      s = c(c = 1, d = 1)),
  m5 = list(f = "R ~ RM - (c + e * V1) * exp(Xt)",    s = c(c = 1, e = 1)),
  m6 = list(f = "R ~ RM - (c + e * V2) * exp(Xt)",    s = c(c = 1, e = 1)),
  m7 = list(f = "R ~ RM - (c + e * V3) * exp(Xt)",    s = c(c = 1, e = 1)),
  m8 = list(f = "R ~ RM - b * exp(Yt) - c * exp(Xt)", s = c(c = 1, b = 1e-4)))
mlab <- c("1 basic", "2 U vertical", "3 SST vertical", "4 ENSO vertical", "5 U lateral",
          "6 SST lateral", "7 ENSO lateral", "8 squid predation")

`%||%` <- function(a, b) if (is.null(a)) b else a

fit_set <- function(B, years, U, label) {
  # B: total biomass (t) for 'years'; R_t needs B_t+1, so R is defined for years[-length(years)]
  d <- data.frame(Year = years, Xt = log(B), R = c(log(B[-1] / B[-length(B)]), NA), V1 = U,
                  V2 = hd$SST0[match(years, hd$year)], V3 = hd$ENSO34[match(years, hd$year)],
                  Yt = hd$log_squid_B[match(years, hd$year)])
  d <- na.omit(subset(d, Year >= 2001 & Year <= 2020))
  rm <- round(max(d$R), 2)
  out <- lapply(names(cand), function(k) {
    f <- as.formula(gsub("RM", sprintf("%.2f", rm), cand[[k]]$f))
    m <- nlsLM(f, data = d, start = cand[[k]]$s, control = nls.lm.control(maxiter = 1000))
    # complete trajectory
    cf <- as.list(coef(m)); X <- d$Xt[1]
    for (t in seq_len(nrow(d) - 1)) {
      env <- c(cf, list(Xt = X[t], V1 = d$V1[t], V2 = d$V2[t], V3 = d$V3[t], Yt = d$Yt[t]))
      X[t + 1] <- X[t] + eval(f[[3]], env)
    }
    r <- resid(m)
    data.frame(series = label, model = mlab[match(k, names(cand))], n = nrow(d), Rm = rm,
               c = cf$c, d = cf$d %||% NA, e = cf$e %||% NA, b = cf$b %||% NA,
               converged = m$convInfo$isConv, AICc = qpcR::AICc(m), R2 = cor(d$R, fitted(m))^2,
               RMSD_logB = sqrt(mean((X - d$Xt)^2)),
               shapiro_p = shapiro.test(r)$p.value, acf1 = acf(r, plot = FALSE)$acf[2])
  })
  tab <- do.call(rbind, out)
  tab$dAICc <- tab$AICc - min(tab$AICc)
  tab$wi <- exp(-tab$dAICc / 2) / sum(exp(-tab$dAICc / 2))
  list(tab = tab, data = d)
}

catch <- hd$catch[match(1997:2021, hd$year)]
runs <- list()
yrs_paper <- hd$year[!is.na(hd$biomass)]
runs[["Paper series"]] <- fit_set(hd$biomass[!is.na(hd$biomass)], yrs_paper, hd$U[!is.na(hd$biomass)], "Paper series")
case_lab <- c("Case 0: official landings", "Case 1: official landings, Francis weighting",
              "Case 2: CCT catches", "Case 3: CCT catches, Francis weighting")
for (k in 0:3) {
  x <- subset(cs, case == k); x <- x[order(x$year), ]
  B <- x$total_biomass_kt * 1000
  runs[[case_lab[k + 1]]] <- fit_set(B, x$year, catch[match(x$year, 1997:2021)] / B, case_lab[k + 1])
  lab_u <- paste0(case_lab[k + 1], " [u of the case table]")
  runs[[lab_u]] <- fit_set(B, x$year, x$u, lab_u)
}

# check: the series of the paper reproduces models_2001_2020.txt (AICc and parameters)
p <- runs[["Paper series"]]$tab
chk <- max(abs(p$AICc - s3$AICc)); chk_c <- max(abs(p$c / s3$c - 1))
cat(sprintf("check vs models_2001_2020.txt: max |dAICc| = %.2e, max rel. diff in c = %.2e\n", chk, chk_c))
stopifnot(chk < 1e-3, chk_c < 1e-3)

all_tab <- do.call(rbind, lapply(runs, `[[`, "tab")); rownames(all_tab) <- NULL
dir.create(file.path(A, "tables"), showWarnings = FALSE)
write.table(transform(all_tab, AICc = round(AICc, 2), dAICc = round(dAICc, 2), wi = round(wi, 3), R2 = round(R2, 3),
                      RMSD_logB = round(RMSD_logB, 3), shapiro_p = round(shapiro_p, 3), acf1 = round(acf1, 2)),
            file.path(A, "tables", "A1_P3_models_by_case.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

summ <- do.call(rbind, lapply(runs, function(r) { t <- r$tab; b <- which.min(t$AICc); s <- t[t$model == "8 squid predation", ]
  data.frame(series = t$series[1], Rm = t$Rm[1], best = t$model[b], w_best = round(t$wi[b], 3), R2_best = round(t$R2[b], 2),
             best_param = signif(c(t$d[b], t$e[b], t$b[b])[!is.na(c(t$d[b], t$e[b], t$b[b]))][1], 3),
             dAICc_squid = round(s$dAICc, 2), w_squid = round(s$wi, 3), b_squid = signif(s$b, 3),
             dAICc_U_lateral = round(t$dAICc[t$model == "5 U lateral"], 2)) }))
rownames(summ) <- NULL
write.table(summ, file.path(A, "tables", "A1_P3_summary.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
print(summ)
saveRDS(runs, file.path(A, "tables", "A1_runs.rds"))
