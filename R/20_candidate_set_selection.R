# Model selection over the full candidate set.
# The candidate set of each period is every model that was fitted: the model structures at every
# lag combination of SST and ENSO (lags 0-3), 243 models in 1968-2020 and in 1975-2000 and 20 in
# 2001-2020. The structures, starting values, fitting routine, seed, Rm, validity check and
# simulation are those of steps 04, 05 and 10.
#   - The RMSD of the complete trajectory is computed for every model of the set, with every
#     predicted year paired with the observation of the same year.
#   - AICc differences and Akaike weights are computed over the valid models of the set.
#   - The lag versions of one structure are models of one type. Following Burnham & Anderson
#     (2002, sections 2.9.2 and 4.6, model redundancy), the weights use prior weights that are
#     equal among structures and, within a structure, equal among its valid lag versions,
#     w_i proportional to exp(-Delta_i / 2) * tau_i with tau_i = 1 / (S * n_s). The weight of a
#     structure is the sum of the weights of its lag versions. The weights with equal priors
#     over all the models (tau_i = 1 / R) are also written, for comparison.
#   - Selection rule: models within two AICc units of the lowest AICc have substantial support
#     (Burnham & Anderson 2002, section 2.6). Among them the model with fewest parameters is
#     taken and, if they have the same number, the one with the lowest RMSD.
#
# Outputs: results/tables/candidate_models_1968_2020.csv, candidate_models_1975_2000.csv and
#          candidate_models_2001_2020.csv (every candidate model),
#          results/tables/model_structures_1968_2020.csv, model_structures_1975_2000.csv and
#          model_structures_2001_2020.csv (one row per structure, at its best valid lag version),
#          results/tables/candidate_models_selection.txt
if (!exists("PROJ")) source(file.path(rprojroot::find_root(rprojroot::has_file("R/config.R")), "R", "config.R"))
suppressPackageStartupMessages({ library(minpack.lm); library(qpcR); library(plyr); library(Metrics); library(stringr) })

hd   <- tsv_read(file.path(PROJ, "data", "derived", "hake_dataset.txt"))
sst  <- tsv_read(file.path(PROJ, "data", "derived", "SSTchileStocks_2022.txt"))
enso <- tsv_read(file.path(PROJ, "data", "derived", "enso34_smoothed.txt"))
LAGS <- 0:3
lagged <- function(x, years, target, k) x[match(target - k, years)]
SST  <- sapply(LAGS, function(k) lagged(sst$CHU, sst$Years, hd$year, k))
ENSO <- sapply(LAGS, function(k) lagged(enso$enso34, enso$Year, hd$year, k))
stopifnot(max(abs(SST[, 1:3] - as.matrix(hd[, c("SST0", "SST1", "SST2")])), na.rm = TRUE) < 1e-12,
          max(abs(ENSO[, 1:3] - as.matrix(hd[, c("ENSO34", "ENSO341", "ENSO342")])), na.rm = TRUE) < 1e-12)
lag_grid <- expand.grid(sst_lag = LAGS, enso_lag = LAGS)
lab <- function(v, k) if (k == 0) sprintf("%s[t]", v) else sprintf("%s[t-%d]", v, k)

# competition term of a fitted model in each year (as in step 10)
k_min <- function(m, f, d) {
  lat <- if (grepl("\\) \\* exp\\(Xt\\)", f)) sub("^.*- \\((.*)\\) \\* exp\\(Xt\\).*$", "\\1", f)
         else sub("^.*- (\\w+) \\* exp\\(Xt\\).*$", "\\1", f)
  min(rep_len(eval(parse(text = lat), c(as.list(coef(m)), as.list(d))), nrow(d)))
}
# RMSD of the complete trajectory on log biomass.
# simula3() starts from the observed biomass of the first year and predicts the biomass of the
# following Tn years, the last of which is the year after the last row of the data. The observed
# biomass of that year is X[Tn] + R[Tn], the value that defines the last R of the fit.
#   RMSD            every predicted year against the observed biomass of the same year (Tn pairs).
#                   The first year is the initial value given to the model, not a prediction, and
#                   is not counted.
#   RMSD_steps_04_05  the expression of steps 04 and 05, kept only to check the reproduction of
#                   their tables. It counts the first year, whose deviation is zero, and pairs the
#                   last observed year with the prediction of the following year.
rmsd2 <- function(m, strings, d) {
  x <- simula3(d$Xt[1], m, strings, d$Year, d$Xt, d$V1, d$V2, d$V3, iter = 200, log = "no", process_error = FALSE)
  Tn <- nrow(d)
  obs <- c(d$Xt, d$Xt[Tn] + d$R[Tn])                       # observed log biomass of years 1 to Tn + 1
  pred <- log(x$Sim2 * 1000)                               # simula3 returns thousand t
  stopifnot(length(pred) == Tn + 1, x$years == c(d$Year, max(d$Year) + 1), abs(pred[1] - obs[1]) < 1e-9)
  c(corrected = sqrt(mean((obs[-1] - pred[-1])^2)),
    steps_04_05 = rmse(log(x$Nt[-length(x$Nt)]), log(x$Sim2[-Tn])))
}

pooled <- function(models, raw, labels, tag, datas, strings) {
  keep <- !duplicated(labels)
  models <- models[keep]; raw <- raw[keep]; labels <- labels[keep]; tag <- tag[keep, ]; datas <- datas[keep]; strings <- strings[keep]
  pars <- rbind.fill(lapply(models, function(m) as.data.frame(t(coef(m)))))
  for (q in c("c", "b", "e", "d")) if (is.null(pars[[q]])) pars[[q]] <- NA
  kmin <- mapply(k_min, models, raw, datas)
  rm2 <- t(mapply(rmsd2, models, strings, datas))
  tab <- data.frame(Models = labels, structure_ID = tag$structure_ID,
                    sst_lag = ifelse(grepl("SST", labels), tag$sst_lag, NA),
                    enso_lag = ifelse(grepl("ENSO", labels), tag$enso_lag, NA),
                    N.param = sapply(models, function(m) length(coef(m))), pars[, c("c", "b", "e", "d")],
                    AICc = sapply(models, qpcR::AICc),
                    R2 = mapply(function(m, d) cor(d$R, predict(m))^2, models, datas),
                    RMSD = rm2[, "corrected"], RMSD_steps_04_05 = rm2[, "steps_04_05"],
                    min_k = kmin, valid = kmin > 0, converged = sapply(models, function(m) m$convInfo$isConv))
  # weights over the valid models
  v <- tab$valid
  tab$n_versions <- as.integer(ave(v, tab$structure_ID, FUN = sum))           # valid lag versions of the structure
  S <- length(unique(tab$structure_ID[v]))
  tab$DAICc <- NA; tab$wi <- NA; tab$wi_equal_priors <- NA
  D <- tab$AICc[v] - min(tab$AICc[v]); L <- exp(-D / 2)
  tau <- 1 / (S * tab$n_versions[v])
  tab$DAICc[v] <- D; tab$wi[v] <- L * tau / sum(L * tau); tab$wi_equal_priors[v] <- L / sum(L)
  tab$W_structure <- ave(ifelse(v, tab$wi, 0), tab$structure_ID, FUN = sum)
  tab$W_structure_equal_priors <- ave(ifelse(v, tab$wi_equal_priors, 0), tab$structure_ID, FUN = sum)
  o <- order(!tab$valid, tab$AICc)
  cbind(rank = seq_along(o), tab[o, ])
}

fit27 <- function(p) {
  per <- PERIODS[[p]]
  rows <- hd$year >= per$from & hd$year <= per$to & !is.na(hd$R) & !is.na(hd$logB) & complete.cases(SST, ENSO)
  base <- data.frame(Year = hd$year, Xt = hd$logB, R = hd$R, V1 = hd$U)[rows, ]
  rm <- round(max(base$R), 2); stopifnot(rm == per$rm)
  models <- list(); raw <- c(); labels <- c(); tag <- data.frame(); datas <- list(); strings <- list()
  for (k in seq_len(nrow(lag_grid))) {
    ls <- lag_grid$sst_lag[k]; le <- lag_grid$enso_lag[k]
    d <- base; d$V2 <- SST[rows, ls + 1]; d$V3 <- ENSO[rows, le + 1]
    set.seed(122)
    am <- pop_models(predictors = 3, Rmax = "obs", Rm = rm, modelpars = "yes", model = "Verhulst", data = d)
    f <- as.character(am[[1]]$Models)
    models <- c(models, am[[2]]); raw <- c(raw, f)
    labels <- c(labels, gsub("V3", lab("ENSO", le), gsub("V2", lab("SST", ls), gsub("V1", "U[t]", f)), fixed = TRUE))
    tag <- rbind(tag, data.frame(structure_ID = seq_along(f), sst_lag = ls, enso_lag = le))
    datas <- c(datas, rep(list(d), length(f)))
    strings <- c(strings, lapply(seq_along(f), function(i) c(am[[3]][[1]][i], am[[3]][[2]][i])))
  }
  list(tab = pooled(models, raw, labels, tag, datas, strings), n = nrow(base), per = per)
}

fit_P3 <- function() {
  per <- PERIODS$P3
  rows <- hd$year >= per$from & hd$year <= per$to & !is.na(hd$R) & !is.na(hd$logB) & !is.na(hd$log_squid_B) & complete.cases(SST, ENSO)
  base <- data.frame(Year = hd$year, Yt = hd$log_squid_B, Xt = hd$logB, R = hd$R, V1 = hd$U)[rows, ]
  rm <- round(max(base$R), 2); stopifnot(rm == per$rm)
  assign("Y", base$Yt, envir = globalenv())                 # the simulation string of the predation model reads Y[t]
  cand <- list(                                             # order, starting values and simulation strings of 05_fit_models_P3.R
    list(f = "R ~ RM - c * exp(Xt)",               s = c(c = 1),           sim = "RM-(c[j])*exp(N)"),
    list(f = "R ~ RM - c * exp(Xt) + d * V1",      s = c(c = 1, d = 1),    sim = "RM-(c[j])*exp(N)+d[j]*V1[t]"),
    list(f = "R ~ RM - c * exp(Xt) + d * V2",      s = c(c = 1, d = 1),    sim = "RM-(c[j])*exp(N)+d[j]*V2[t]"),
    list(f = "R ~ RM - c * exp(Xt) + d * V3",      s = c(c = 1, d = 1),    sim = "RM-(c[j])*exp(N)+d[j]*V3[t]"),
    list(f = "R ~ RM - (c + e * V1) * exp(Xt)",    s = c(c = 1, e = 1),    sim = "RM-(c[j]+e[j]*V1[t])*exp(N)"),
    list(f = "R ~ RM - (c + e * V2) * exp(Xt)",    s = c(c = 1, e = 1),    sim = "RM-(c[j]+e[j]*V2[t])*exp(N)"),
    list(f = "R ~ RM - (c + e * V3) * exp(Xt)",    s = c(c = 1, e = 1),    sim = "RM-(c[j]+e[j]*V3[t])*exp(N)"),
    list(f = "R ~ RM - b * exp(Yt) - c * exp(Xt)", s = c(c = 1, b = 1e-4), sim = "RM-b[j]*exp(Y[t])-(c[j])*exp(N)"))
  sim_strings <- function(k) c(gsub("N", "N1a", gsub("RM", sprintf("%.2f", rm), k$sim)),
                               gsub("N", "N1b", gsub("RM", sprintf("%.2f", rm), k$sim)))
  models <- list(); raw <- c(); labels <- c(); tag <- data.frame(); datas <- list(); strings <- list()
  for (k in seq_len(nrow(lag_grid))) {
    ls <- lag_grid$sst_lag[k]; le <- lag_grid$enso_lag[k]
    d <- base; d$V2 <- SST[rows, ls + 1]; d$V3 <- ENSO[rows, le + 1]
    for (i in seq_along(cand)) {
      f <- sub("RM", sprintf("%.2f", rm), cand[[i]]$f)
      models <- c(models, list(nlsLM(as.formula(f), data = d, start = cand[[i]]$s,
                                     control = nls.lm.control(maxiter = 1000), algorithm = "BFGS")))
      raw <- c(raw, f)
      labels <- c(labels, gsub("Yt", "Y[t]", gsub("V3", lab("ENSO", le), gsub("V2", lab("SST", ls), gsub("V1", "U[t]", f)), fixed = TRUE)))
      tag <- rbind(tag, data.frame(structure_ID = i, sst_lag = ls, enso_lag = le)); datas <- c(datas, list(d))
      strings <- c(strings, list(sim_strings(cand[[i]])))
    }
  }
  list(tab = pooled(models, raw, labels, tag, datas, strings), n = nrow(base), per = per)
}

# one row per structure: its valid lag version with the lowest AICc
by_structure <- function(tab) {
  v <- subset(tab, valid)
  best <- do.call(rbind, lapply(split(v, v$structure_ID), function(x) x[which.min(x$AICc), ]))
  total <- table(tab$structure_ID)[as.character(best$structure_ID)]
  out <- data.frame(ID = best$structure_ID, Models = best$Models, lag_versions = sprintf("%d of %d", best$n_versions, as.integer(total)),
                    best[, c("N.param", "c", "b", "e", "d", "AICc", "DAICc", "wi", "W_structure", "R2", "RMSD", "RMSD_steps_04_05")])
  missing <- setdiff(unique(tab$structure_ID), out$ID)
  if (length(missing)) cat("structures without a valid lag version:", paste(missing, collapse = ", "), "\n")
  out[order(out$ID), ]
}

select_model <- function(tab) {
  near <- subset(tab, valid & DAICc < 2)
  near <- near[order(near$N.param, near$RMSD), ]
  near[1, ]
}

summarise <- function(title, x, old_file, old_id) {
  tab <- x$tab; v <- subset(tab, valid); near <- subset(v, DAICc < 2); sel <- select_model(tab)
  st <- unique(v[order(-v$W_structure), c("structure_ID", "n_versions", "W_structure", "W_structure_equal_priors")])
  old <- tsv_read(file.path(RES, "tables", old_file))                        # check against the tables of steps 04 and 05
  same <- v[v$Models == gsub("~", "~", old$Models[old$ID == old_id], fixed = TRUE), ]
  c(sprintf("== %s, n = %d", title, x$n),
    sprintf("models fitted: %d; converged: %d; valid: %d; structures: %d; models with DAICc < 2: %d",
            nrow(tab), sum(tab$converged), nrow(v), length(unique(tab$structure_ID)), nrow(near)),
    "models with DAICc < 2 (evidence ratio = w of the lowest AICc model / w of the model):",
    sprintf("   %-72s K %d  DAICc %.2f  ER %.2f  wi %.3f  R2 %.3f  RMSD %.4f",
            near$Models, near$N.param, near$DAICc, exp(near$DAICc / 2), near$wi, near$R2, near$RMSD),
    sprintf("selected (fewest parameters, then lowest RMSD): %s  [structure %d, DAICc %.2f, wi %.3f, W structure %.3f, R2 %.3f, RMSD %.4f]",
            sel$Models, sel$structure_ID, sel$DAICc, sel$wi, sel$W_structure, sel$R2, sel$RMSD),
    sprintf("model of %s (ID %d) in this set: AICc %.3f (table %.3f), RMSD with the expression of that step %.4f (table %.4f)", old_file, old_id,
            same$AICc, old$AICc[old$ID == old_id], same$RMSD_steps_04_05, old$RMSD[old$ID == old_id]),
    "weight of the structures (priors of section 4.6 | equal priors over all models), first five:",
    sprintf("   structure %2d  lag versions %2d  W %.3f | %.3f", st$structure_ID, st$n_versions, st$W_structure, st$W_structure_equal_priors)[1:min(5, nrow(st))],
    sprintf("95%% confidence set of models (section 4.5): %d models", sum(cumsum(sort(v$wi, decreasing = TRUE)) < 0.95) + 1), "")
}

out <- c("Model selection over the full candidate set (R/20_candidate_set_selection.R)",
         "RMSD: every predicted year against the observed log biomass of the same year", "")
for (p in c("P1", "P2", "P3")) {
  x <- if (p == "P3") fit_P3() else fit27(p)
  i <- match(p, c("P1", "P2", "P3"))
  write.csv(x$tab, file.path(RES, "tables", sprintf("candidate_models_%d_%d.csv", x$per$from, x$per$to)), row.names = FALSE)
  write.csv(by_structure(x$tab), file.path(RES, "tables", sprintf("model_structures_%d_%d.csv", x$per$from, x$per$to)), row.names = FALSE)
  lines <- summarise(x$per$label, x, sprintf("models_%d_%d.txt", x$per$from, x$per$to), c(18, 7, 5)[i])
  writeLines(lines); out <- c(out, lines)
}
writeLines(out, file.path(RES, "tables", "candidate_models_selection.txt"))
