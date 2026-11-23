# Tests of the 2001-2020 results against the acoustic-survey biomass, which is
# estimated by echo-integration and does not depend on the assessment model or its squid term.
#
# Survey biomass (winter surveys, thousand t -> t): Table 10 of the 2024 survey report (2002, 2004-2024;
# data/survey_biomass_table10_2024.csv) and, for 1995, 1997, 1999-2001, the values digitised from
# Figure 3 of the IFOP 2022 status report (data/survey_biomass_fig3_IFOP2022_digitised.csv; check against
# Table 10: mean absolute error 4 kt). No surveys in 1994, 1996, 1998, 2003.
#
# B1 (as Lima et al. 2020, Front. Mar. Sci. 7:123): the 2001-2020 models fitted to the assessment biomass
#   of the paper predict the biomass trajectory (deterministic complete trajectory from the 2001 biomass,
#   observed predictors); the prediction is compared with the survey biomass of 2004-2020 (the paper
#   series ends in 2020). Because the
#   survey measures biomass on its own scale (catchability q), agreement is measured on log biomass by
#   (i) the correlation and (ii) the RMSD after removing the mean log difference (log q).
# B2: the same 8 candidate models fitted directly to the survey biomass, R_t = log(B_t+1/B_t) for
#   t = 2004-2020 (consecutive surveys; the 2001-2002 and 2002-2004 steps are not usable because of the
#   2003 gap). U = official catch (hake_dataset) / survey biomass. Rm = round(max R_t, 2) as in the paper.
# Outputs: tables/B1_validation.txt, tables/B2_models_survey.txt, tables/survey_biomass_series.txt,
#          figs/B1_*.png, figs/B2_*.png
suppressPackageStartupMessages({ library(minpack.lm); library(qpcR) })
B <- rprojroot::find_root(rprojroot::has_file("R/B1_survey_tests.R"))
PROJ <- B; RES <- B
HAKE <- normalizePath(file.path(B, "..", ".."))   # repository root 
source(file.path(B, "R", "lib", "fig_theme.R"))

hd <- read.delim(file.path(HAKE, "data", "derived", "hake_dataset.txt"))
t10 <- read.csv(file.path(B, "data", "survey_biomass_table10_2024.csv"))
dg  <- read.csv(file.path(B, "data", "survey_biomass_fig3_IFOP2022_digitised.csv"))
sv <- rbind(data.frame(year = t10$year, survey_t = t10$biomass_t, source = "Table 10, IFOP survey report 2024"),
            data.frame(year = dg$year[dg$year < 2002], survey_t = dg$biomass_kt_digitised[dg$year < 2002] * 1000,
                       source = "digitised, Fig. 3 IFOP 2022 status report"))
sv <- sv[order(sv$year), ]
write.table(sv, file.path(B, "tables", "survey_biomass_series.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

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
pred_vars <- function(years, Bt) data.frame(V1 = hd$catch[match(years, hd$year)] / Bt, V2 = hd$SST0[match(years, hd$year)],
                                            V3 = hd$ENSO34[match(years, hd$year)], Yt = hd$log_squid_B[match(years, hd$year)])
fit_all <- function(d) {
  rm <- round(max(d$R), 2)
  lapply(setNames(names(cand), mlab), function(k) {
    f <- as.formula(gsub("RM", sprintf("%.2f", rm), cand[[k]]$f))
    nlsLM(f, data = d, start = cand[[k]]$s, control = nls.lm.control(maxiter = 1000))
  })
}
traj <- function(m, X0, P) {             # deterministic complete trajectory of log biomass
  f <- formula(m); cf <- as.list(coef(m)); X <- X0
  for (t in seq_len(nrow(P) - 1)) X[t + 1] <- X[t] + eval(f[[3]], c(cf, list(Xt = X[t]), as.list(P[t, ])))
  X
}

# ---------------- B1: models fitted to the assessment biomass of the paper, checked against the survey
yrs <- 2001:2021
Bp <- hd$biomass[match(yrs, hd$year)]
dp <- data.frame(Year = yrs, Xt = log(Bp), R = c(log(Bp[-1] / Bp[-length(Bp)]), NA), pred_vars(yrs, Bp))
dp$V1 <- hd$U[match(yrs, hd$year)]                       # the paper's U (identical definition)
mods_p <- fit_all(na.omit(subset(dp, Year <= 2020)))
Pp <- dp[, c("V1", "V2", "V3", "Yt")]; Pp$Yt[is.na(Pp$Yt)] <- 0
pred <- sapply(mods_p, function(m) traj(m, dp$Xt[1], Pp))     # 2001-2021
s_y <- sv$year[sv$year %in% 2004:2020]   # the paper series ends in 2020
s_x <- log(sv$survey_t[match(s_y, sv$year)])
val <- do.call(rbind, lapply(c(colnames(pred), "Assessment biomass (paper)"), function(k) {
  x <- if (k == "Assessment biomass (paper)") dp$Xt[match(s_y, yrs)] else pred[match(s_y, yrs), k]
  lq <- mean(s_x - x)
  data.frame(series = k, n = length(s_y), r_logB = cor(x, s_x), log_q = lq, q = exp(lq),
             RMSD_logB_after_q = sqrt(mean((s_x - x - lq)^2)), AICc_on_assessment =
               if (k %in% names(mods_p)) qpcR::AICc(mods_p[[k]]) else NA)
}))
val$dAICc_on_assessment <- val$AICc_on_assessment - min(val$AICc_on_assessment, na.rm = TRUE)
write.table(transform(val, r_logB = round(r_logB, 3), log_q = round(log_q, 3), q = round(q, 3),
                      RMSD_logB_after_q = round(RMSD_logB_after_q, 3), AICc_on_assessment = round(AICc_on_assessment, 2),
                      dAICc_on_assessment = round(dAICc_on_assessment, 2)),
            file.path(B, "tables", "B1_validation.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nB1: assessment-fitted models vs survey biomass 2004-2020\n"); print(val, digits = 3)

# ---------------- B2: models fitted to the survey biomass
yrs2 <- 2004:2021
Bs <- sv$survey_t[match(yrs2, sv$year)]
ds <- data.frame(Year = yrs2, Xt = log(Bs), R = c(log(Bs[-1] / Bs[-length(Bs)]), NA), pred_vars(yrs2, Bs))
ds <- na.omit(subset(ds, Year <= 2020))
mods_s <- fit_all(ds)
tab2 <- do.call(rbind, lapply(names(mods_s), function(k) { m <- mods_s[[k]]; cf <- as.list(coef(m)); r <- resid(m)
  data.frame(model = k, n = nrow(ds), Rm = round(max(ds$R), 2), c = cf$c, d = cf$d %||% NA, e = cf$e %||% NA, b = cf$b %||% NA,
             converged = m$convInfo$isConv, AICc = qpcR::AICc(m), R2 = cor(ds$R, fitted(m))^2,
             shapiro_p = shapiro.test(r)$p.value, acf1 = acf(r, plot = FALSE)$acf[2]) }))
tab2$dAICc <- tab2$AICc - min(tab2$AICc); tab2$wi <- exp(-tab2$dAICc / 2) / sum(exp(-tab2$dAICc / 2))
write.table(transform(tab2, AICc = round(AICc, 2), dAICc = round(dAICc, 2), wi = round(wi, 3), R2 = round(R2, 3),
                      shapiro_p = round(shapiro_p, 3), acf1 = round(acf1, 2)),
            file.path(B, "tables", "B2_models_survey.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nB2: models fitted to the survey biomass, R_t 2004-2020\n"); print(tab2[, c("model", "n", "Rm", "AICc", "dAICc", "wi", "R2", "acf1")], digits = 3)

# ---------------- figures
ab <- data.frame(year = hd$year, B = hd$biomass / 1000); ab <- subset(ab, year >= 1993 & year <= 2021)
f1 <- ggplot() +
  geom_line(data = ab, aes(year, B, colour = "Assessment (paper series)"), linewidth = 0.8) +
  geom_point(data = sv, aes(year, survey_t / 1000, colour = "Acoustic survey", shape = grepl("digitised", source)), size = 1.8) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), labels = c("Survey, Table 10 (2024 report)", "Survey, digitised (IFOP 2022 Fig. 3)")) +
  scale_colour_manual(values = c("Assessment (paper series)" = "black", "Acoustic survey" = "#0072B2"), guide = "none") +
  scale_x_continuous(breaks = seq(1993, 2024, 4)) + scale_y_continuous(limits = c(0, NA), expand = expansion(c(0, 0.05))) +
  annotate("text", x = 1998.3, y = 1250, label = "Assessment biomass\n(paper series)", hjust = 0, size = 2.8) +
  labs(x = "Year", y = "Total biomass (kt)") + theme_pub() + theme(legend.position = c(0.75, 0.85))
save_fig(f1, "B1_survey_vs_assessment", 12, 8, dir = file.path(B, "figs"))

show <- c("1 basic", "5 U lateral", "8 squid predation")
pl <- do.call(rbind, lapply(show, function(k) data.frame(year = yrs, model = k, B = exp(pred[, k] + val$log_q[val$series == k]) / 1000)))
f2 <- ggplot() +
  geom_line(data = pl, aes(year, B, colour = model), linewidth = 0.7) +
  geom_point(data = subset(sv, year >= 2001), aes(year, survey_t / 1000), size = 1.6) +
  scale_colour_manual(values = c("1 basic" = "grey50", "5 U lateral" = "#0072B2", "8 squid predation" = "#C0392B")) +
  coord_cartesian(ylim = c(0, 1000)) +
  labs(x = "Year", y = "Biomass (kt; model prediction scaled by q)",
       subtitle = "Lines: models fitted to the assessment biomass 2001-2020, complete trajectory from 2001\nPoints: acoustic survey (2002 = 1555 kt, off scale)") +
  theme_pub() + theme(legend.position = c(0.75, 0.8), plot.subtitle = element_text(size = 7))
save_fig(f2, "B1_predicted_vs_survey", 12, 9, dir = file.path(B, "figs"))

tab2$model <- factor(tab2$model, rev(mlab))
f3 <- ggplot(tab2, aes(dAICc, model, fill = model == "8 squid predation")) + geom_col(width = 0.7) +
  geom_vline(xintercept = 2, linetype = "dashed", linewidth = 0.3) +
  scale_fill_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "grey60"), guide = "none") +
  labs(x = expression(Delta * "AICc, models fitted to the survey biomass (R"[t] * " 2004-2020)"), y = NULL) + theme_pub()
save_fig(f3, "B2_survey_models_dAICc", 11, 7, dir = file.path(B, "figs"))
saveRDS(list(assessment = mods_p, survey = mods_s, ds = ds), file.path(B, "tables", "B_models.rds"))
