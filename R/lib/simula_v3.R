# simula3(): complete-trajectory simulation of a fitted R-function model with propagation of
# uncertainty.
#   - Every parameter draw j (multivariate normal with the estimates and their covariance matrix)
#     runs its own trajectory from N0, so the interval widens with time.
#   - With process_error = TRUE a second set of trajectories adds a N(0, sigma^2) innovation to R_t
#     each year (sigma = residual standard error of the nls fit), which gives a prediction interval.
# Arguments: initial value No, fitted model Mobj, the two simulation strings formN, years, observed
#   log biomass Nt, predictors V1 to V3, number of draws iter, process_error (TRUE/FALSE) and seed.
# Returns a data frame with years, Nt, Sim1, Sim2, Sim1.CI.L/U, Sim2.CI.L/U and, when
# process_error = TRUE, Sim2.PI.L / Sim2.PI.U (parameter + process error).
# Sim1 = one-step-ahead prediction from the observed value; Sim2 = complete trajectory.
simula3 <- function(No, Mobj, formN, years, Nt, V1, V2, V3, iter = 5000,
                    ylab = "", xlab = "", plot = "no", log = "no",
                    process_error = TRUE, seed = 6) {
  require(MASS)
  set.seed(seed)
  param <- c(Mobj$m$getAllPars())
  pr <- mvrnorm(iter, param, vcov(Mobj))
  pr <- matrix(pr, ncol = length(param), dimnames = list(NULL, names(param)))
  sigma <- summary(Mobj)$sigma
  Tn <- length(Nt)
  j <- seq_len(iter)                       # the formula strings index parameters as c[j], b[j], ...
  env <- list2env(as.data.frame(pr), parent = environment())

  N2  <- matrix(NA_real_, Tn, iter)        # complete trajectory, parameter uncertainty only
  N2p <- matrix(NA_real_, Tn, iter)        # complete trajectory, parameter + process error
  N3  <- matrix(NA_real_, Tn, iter)        # one-step-ahead from the observed value
  prev  <- rep(No, iter)
  prevp <- rep(No, iter)
  for (t in seq_len(Tn)) {
    assign("N1a", prev, envir = env)
    N2[t, ]  <- prev + eval(parse(text = formN[1]), envir = env)
    assign("N1a", prevp, envir = env)
    eps <- if (process_error) rnorm(iter, 0, sigma) else 0
    N2p[t, ] <- prevp + eval(parse(text = formN[1]), envir = env) + eps
    assign("N1b", Nt[t], envir = env)
    N3[t, ]  <- Nt[t] + eval(parse(text = formN[2]), envir = env)
    prev  <- N2[t, ]
    prevp <- N2p[t, ]
  }
  q <- function(m) t(apply(m, 1, quantile, c(0.025, 0.5, 0.975), na.rm = TRUE))
  qqN2 <- q(N2); qqN2p <- q(N2p); qqN3 <- q(N3)

  f <- if (log == "yes") function(x) x else function(x) exp(x) / 1000
  out <- data.frame(years = seq(min(years), max(years) + 1, 1),
                    Nt = c(f(Nt), NA)[1:(Tn + 1)],
                    Sim1 = f(c(No, qqN3[, 2])), Sim2 = f(c(No, qqN2[, 2])),
                    Sim1.CI.L = f(c(No, qqN3[, 1])), Sim1.CI.U = f(c(No, qqN3[, 3])),
                    Sim2.CI.L = f(c(No, qqN2[, 1])), Sim2.CI.U = f(c(No, qqN2[, 3])))
  if (process_error) {
    out$Sim2.PI.L <- f(c(No, qqN2p[, 1])); out$Sim2.PI.U <- f(c(No, qqN2p[, 3]))
    out$Sim2.PI.median <- f(c(No, qqN2p[, 2]))
  }
  attr(out, "sigma") <- sigma
  if (plot == "yes") {
    plot(out$years, out$Nt, pch = 19, cex = .6, xlab = xlab, ylab = ylab,
         ylim = range(c(out$Nt, out$Sim2.CI.L, out$Sim2.CI.U, out$Sim2.PI.L, out$Sim2.PI.U), na.rm = TRUE))
    if (process_error) polygon(c(out$years, rev(out$years)), c(out$Sim2.PI.L, rev(out$Sim2.PI.U)), col = rgb(0, 0, 1, .1), border = NA)
    polygon(c(out$years, rev(out$years)), c(out$Sim2.CI.L, rev(out$Sim2.CI.U)), col = rgb(0, 0, 1, .2), border = NA)
    lines(out$years, out$Sim2, col = "blue", lwd = 2)
  }
  out
}
