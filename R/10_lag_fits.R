# Lag tables: every model structure fitted at every lag combination of SST and ENSO.
#   - SST and ENSO at lags of 0, 1, 2 and 3 years (16 lag combinations);
#   - the three periods: 1968-2020, 1975-2000 and 2001-2020;
#   - every fitted model is checked for the sign of its competition term, k_t = the lateral
#     term of the model evaluated in each year (c, c + e SST_t + d ENSO_t, ...). A model with
#     k_t <= 0 in any year of its period has no equilibrium in that year and is flagged invalid;
#   - Akaike weights are given over all unique models and over the valid ones.
# The candidate structures, starting values, fitting routine, seed and Rm are those of the main
# analysis (pop_models() for the 27 structures; the eight structures of 05_fit_models_P3.R).
# For P3 two sets are fitted: (a) the eight structures of the paper (at most two parameters,
# with the predation model), SST and ENSO at lags 0-3; (b) the 27 structures of P1 and P2, which
# is outside the rule of the paper for this period (n = 19) and is reported for completeness.
# U is never lagged. All lag combinations use the same years, so AICc values are comparable.
#
# Outputs: results/tables/lag_fits_1968_2020.csv, lag_fits_1975_2000.csv, lag_fits_2001_2020.csv,
#          lag_fits_2001_2020_27_structures.csv, lag_fits_summary.txt
if (!exists("PROJ")) source(file.path(rprojroot::find_root(rprojroot::has_file("R/config.R")), "R", "config.R"))
suppressPackageStartupMessages({ library(minpack.lm); library(qpcR); library(plyr) })

hd   <- tsv_read(file.path(PROJ, "data", "derived", "hake_dataset.txt"))
sst  <- tsv_read(file.path(PROJ, "data", "derived", "SSTchileStocks_2022.txt"))
enso <- tsv_read(file.path(PROJ, "data", "derived", "enso34_smoothed.txt"))
LAGS <- 0:3
lagged <- function(x, years, target, k) x[match(target - k, years)]
SST  <- sapply(LAGS, function(k) lagged(sst$CHU, sst$Years, hd$year, k))
ENSO <- sapply(LAGS, function(k) lagged(enso$enso34, enso$Year, hd$year, k))
# lags 0-2 must be the columns of the analysis dataset
stopifnot(max(abs(SST[, 1:3] - as.matrix(hd[, c("SST0", "SST1", "SST2")])), na.rm = TRUE) < 1e-12,
          max(abs(ENSO[, 1:3] - as.matrix(hd[, c("ENSO34", "ENSO341", "ENSO342")])), na.rm = TRUE) < 1e-12)
lag_grid <- expand.grid(sst_lag = LAGS, enso_lag = LAGS)
lab <- function(v, k) if (k == 0) sprintf("%s[t]", v) else sprintf("%s[t-%d]", v, k)

# competition term of a fitted model in each year: text between "- (" and ") * exp(Xt)", or the
# single coefficient that multiplies exp(Xt)
k_min <- function(m, f, d) {
  lat <- if (grepl("\\) \\* exp\\(Xt\\)", f)) sub("^.*- \\((.*)\\) \\* exp\\(Xt\\).*$", "\\1", f)
         else sub("^.*- (\\w+) \\* exp\\(Xt\\).*$", "\\1", f)
  min(rep_len(eval(parse(text = lat), c(as.list(coef(m)), as.list(d))), nrow(d)))
}
resid_diag <- function(m) {
  r <- residuals(m)
  sprintf("converged %s, Shapiro-Wilk P = %.2f, residual lag-1 autocorrelation = %.2f",
          m$convInfo$isConv, shapiro.test(r)$p.value, acf(r, plot = FALSE)$acf[2])
}

# pooled table of one period from lists of models, raw formulas, labels and data
pooled <- function(models, raw, labels, tag, datas) {
  keep <- !duplicated(labels)
  models <- models[keep]; raw <- raw[keep]; labels <- labels[keep]; tag <- tag[keep, ]; datas <- datas[keep]
  aic <- sapply(models, qpcR::AICc); w <- akaike.weights(aic)
  pars <- rbind.fill(lapply(models, function(m) as.data.frame(t(coef(m)))))
  for (q in c("c", "b", "e", "d")) if (is.null(pars[[q]])) pars[[q]] <- NA
  kmin <- mapply(k_min, models, raw, datas)
  tab <- data.frame(Models = labels, structure_ID = tag$structure_ID,
                    sst_lag = ifelse(grepl("SST", labels), tag$sst_lag, NA),
                    enso_lag = ifelse(grepl("ENSO", labels), tag$enso_lag, NA),
                    N.param = sapply(models, function(m) length(coef(m))), pars[, c("c", "b", "e", "d")],
                    AICc = aic, DAICc = w$deltaAIC, wi = w$weights,
                    R2 = mapply(function(m, d) cor(d$R, predict(m))^2, models, datas),
                    min_k = kmin, valid = kmin > 0, converged = sapply(models, function(m) m$convInfo$isConv))
  tab$DAICc_valid <- NA; tab$wi_valid <- NA
  v <- which(tab$valid); wv <- akaike.weights(tab$AICc[v])
  tab$DAICc_valid[v] <- wv$deltaAIC; tab$wi_valid[v] <- wv$weights
  o <- order(tab$AICc)
  list(tab = cbind(rank = seq_along(o), tab[o, ]), models = models[o])
}

summarise <- function(p, title, res, n, main_id) {
  tab <- res$tab; v <- subset(tab, valid); best <- v[1, ]
  near <- subset(v, DAICc_valid < 2)
  wl <- function(col) paste(sprintf("%.2f", sapply(LAGS, function(k) sum(v$wi_valid[v[[col]] %in% k]))), collapse = " / ")
  pre <- subset(v, structure_ID == main_id)
  pre <- pre[order(pre$AICc), c("Models", "AICc", "DAICc_valid", "wi_valid", "R2")]
  c(sprintf("== %s (%s), n = %d", p, title, n),
    sprintf("unique models: %d; converged: %d; valid (k_t > 0 in every year): %d; valid with dAICc < 2: %d",
            nrow(tab), sum(tab$converged), nrow(v), nrow(near)),
    sprintf("best of all models: %s (structure %d, AICc %.2f, valid %s)", tab$Models[1], tab$structure_ID[1], tab$AICc[1], tab$valid[1]),
    sprintf("best valid model:   %s (structure %d, AICc %.2f, wi %.3f, R2 %.2f)", best$Models, best$structure_ID, best$AICc, best$wi_valid, best$R2),
    sprintf("   diagnostics of the best valid model: %s", resid_diag(res$models[[which(tab$valid)[1]]])),
    "valid models with dAICc < 2:",
    sprintf("   %-78s dAICc %.2f  wi %.3f  R2 %.2f", near$Models, near$DAICc_valid, near$wi_valid, near$R2),
    sprintf("summed wi (valid models) by SST lag  0/1/2/3: %s", wl("sst_lag")),
    sprintf("summed wi (valid models) by ENSO lag 0/1/2/3: %s", wl("enso_lag")),
    sprintf("structure of the paper (model %d) by lag, best five:", main_id),
    capture.output(print(head(pre, 5), row.names = FALSE, digits = 4)), "")
}

fit27 <- function(p) {
  per <- PERIODS[[p]]
  rows <- hd$year >= per$from & hd$year <= per$to & !is.na(hd$R) & !is.na(hd$logB) & complete.cases(SST, ENSO)
  if (p == "P3") rows <- rows & !is.na(hd$log_squid_B)          # same years as 05_fit_models_P3.R
  base <- data.frame(Year = hd$year, Xt = hd$logB, R = hd$R, V1 = hd$U)[rows, ]
  rm <- round(max(base$R), 2); stopifnot(rm == per$rm)
  models <- list(); raw <- c(); labels <- c(); tag <- data.frame(); datas <- list()
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
  }
  list(res = pooled(models, raw, labels, tag, datas), n = nrow(base), per = per)
}

fit_P3_paper <- function() {
  per <- PERIODS$P3
  rows <- hd$year >= per$from & hd$year <= per$to & !is.na(hd$R) & !is.na(hd$logB) & !is.na(hd$log_squid_B) & complete.cases(SST, ENSO)
  base <- data.frame(Year = hd$year, Yt = hd$log_squid_B, Xt = hd$logB, R = hd$R, V1 = hd$U)[rows, ]
  rm <- round(max(base$R), 2); stopifnot(rm == per$rm)
  cand <- list(                                   # order and starting values of 05_fit_models_P3.R
    list(f = "R ~ RM - c * exp(Xt)",               s = c(c = 1)),
    list(f = "R ~ RM - c * exp(Xt) + d * V1",      s = c(c = 1, d = 1)),
    list(f = "R ~ RM - c * exp(Xt) + d * V2",      s = c(c = 1, d = 1)),
    list(f = "R ~ RM - c * exp(Xt) + d * V3",      s = c(c = 1, d = 1)),
    list(f = "R ~ RM - (c + e * V1) * exp(Xt)",    s = c(c = 1, e = 1)),
    list(f = "R ~ RM - (c + e * V2) * exp(Xt)",    s = c(c = 1, e = 1)),
    list(f = "R ~ RM - (c + e * V3) * exp(Xt)",    s = c(c = 1, e = 1)),
    list(f = "R ~ RM - b * exp(Yt) - c * exp(Xt)", s = c(c = 1, b = 1e-4)))
  models <- list(); raw <- c(); labels <- c(); tag <- data.frame(); datas <- list()
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
    }
  }
  list(res = pooled(models, raw, labels, tag, datas), n = nrow(base), per = per)
}

out <- c("Lag fits: SST and ENSO at lags 0-3, three periods (R/10_lag_fits.R)", "")
for (p in c("P1", "P2", "P3")) {
  x <- fit27(p)
  main_id <- if (p == "P3") 2 else readRDS(file.path(RES, "models", paste0("best_", p, ".rds")))$ID
  name <- if (p == "P3") "lag_fits_2001_2020_27_structures.csv" else sprintf("lag_fits_%d_%d.csv", x$per$from, x$per$to)
  write.csv(x$res$tab, file.path(RES, "tables", name), row.names = FALSE)
  title <- if (p == "P3") paste(x$per$label, "- 27 structures, outside the rule of the paper for this period") else x$per$label
  lines <- summarise(p, title, x$res, x$n, main_id)
  if (p == "P3") lines[grepl("^structure of the paper", lines)] <- "structure 2 (U lateral, c + b U_t; no SST or ENSO, a single model):"
  writeLines(lines); out <- c(out, lines)
}
x <- fit_P3_paper()
write.csv(x$res$tab, file.path(RES, "tables", "lag_fits_2001_2020.csv"), row.names = FALSE)
lines <- summarise("P3", paste(x$per$label, "- the eight structures of the paper"), x$res, x$n, 5)
writeLines(lines); out <- c(out, lines)
writeLines(out, file.path(RES, "tables", "lag_fits_summary.txt"))
