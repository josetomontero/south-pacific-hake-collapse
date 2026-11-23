# Figure 3. Left (a-c): observed total biomass and the complete-trajectory prediction of the best
# model of each period (simula3: median of 5000 parameter draws; dark band = 95% interval from
# parameter uncertainty; light band = 95% prediction interval adding the residual process error).
# Right (d-f): observed reproductive rate R_t and the bootstrap prediction of the same models.
if (!exists("PROJ")) source(file.path(rprojroot::find_root(rprojroot::has_file("R/config.R")), "R", "config.R"))
source(file.path(PROJ, "R", "lib", "fig_theme.R"))

d <- subset(tsv_read(file.path(PROJ, "data", "derived", "hake_dataset.txt")), year <= 2020)
d$biomass_kt <- d$biomass / 1000
xlim <- c(1966, 2021)
labels <- c(P1 = "1968–2020", P2 = "1975–2000", P3 = "2001–2020")
PT_IN <- 2.0; PT_OUT <- 1.6          # point sizes: period analysed / other years
PT_STROKE <- 0.35                    # thin white rim on the black points so that neighbours do not merge (as in Figure 2)
BOOT_COL <- "#0072B2"                # bootstrap 95% CI band of R_t (panels d-f); blue so that it differs from the parameter CI of a-c

# Axis weight: thin (0.3) looks lighter; the previous heavier setting is kept for comparison.
AXIS_LW <- 0.3      # thin axes
# AXIS_LW <- 0.5    # original heavier axes
theme_nat <- function(show_x, show_y) theme_pub(base_size = 9) +
  theme(axis.ticks.length = unit(0.12, "cm"), axis.line = element_line(linewidth = AXIS_LW),
        axis.ticks = element_line(linewidth = AXIS_LW),
        plot.margin = margin(6, 8, 2, 4),
        # every panel keeps both axis titles (invisible where not needed) so all panels have the same size
        axis.title.x = element_text(colour = if (show_x) "black" else "transparent"),
        axis.title.y = element_text(colour = if (show_y) "black" else "transparent"))
breaks_layer <- geom_vline(xintercept = c(1974, 2001), linetype = "22", colour = "grey55", linewidth = 0.4)
# Inputs and output name. 08c_Fig3_lag_selected.R sets FIG3_TAG (suffix of the input files of a period) and
# FIG3_NAME before sourcing this script; without them the figure is the one of steps 04-06 (Fig_3).
if (!exists("FIG3_TAG")) FIG3_TAG <- c(P1 = "", P2 = "", P3 = "")
if (!exists("FIG3_NAME")) FIG3_NAME <- "Fig_3"
sims <- setNames(lapply(names(PERIODS), function(p) tsv_read(file.path(RES, "models", paste0("sim_", p, FIG3_TAG[[p]], ".txt")))), names(PERIODS))
YB <- c(0, max(sapply(sims, function(s) max(s$Sim2.PI.U))))   # biomass axis ends at the maximum 95% PI (+ process) of the three models
YR <- c(-0.6, 0.5)          # R_t axis: covers observations and bootstrap bands
fy <- function(lim, f) lim[1] + f * diff(lim)      # position as a fraction of the axis range (same height in all panels)

# letter top flush with the upper end of the y axis (vjust = 1 at the axis maximum)
panel_letter <- function(letter, lim) annotate("text", x = xlim[1] + 0.4, y = lim[2], label = paste0(letter, ")"),
                                               hjust = 0, vjust = 1.05, size = 3.4)

panels <- list()
for (i in seq_along(PERIODS)) {
  p <- names(PERIODS)[i]
  s <- sims[[p]]
  # simula3() leaves the observed biomass of the last predicted year empty. When that year belongs to the period
  # (2020 in 1968-2020 and in 2001-2020) its biomass is observed and modeled, so it is drawn as a modeled observation.
  # In 1975-2000 the last predicted year is 2001, outside the period, and it stays as an observation of another year.
  last <- is.na(s$Nt); stopifnot(sum(last) == 1, which(last) == nrow(s))
  if (s$years[last] <= PERIODS[[p]]$to) s$Nt[last] <- d$biomass_kt[match(s$years[last], d$year)]
  b <- tsv_read(file.path(RES, "tables", paste0("bootstrap_R_", p, FIG3_TAG[[p]], ".txt")))

  panels[[2 * i - 1]] <- ggplot() + breaks_layer +
    geom_ribbon(data = s, aes(years, ymin = Sim2.PI.L, ymax = Sim2.PI.U), fill = PRED_COL, alpha = 0.12) +
    geom_ribbon(data = s, aes(years, ymin = Sim2.CI.L, ymax = Sim2.CI.U), fill = PRED_COL, alpha = 0.30) +
    geom_line(data = s, aes(years, Sim2), colour = PRED_COL, linewidth = 0.9) +
    geom_point(data = d, aes(year, biomass_kt), shape = 21, fill = "white", colour = "grey45", size = PT_OUT, stroke = 0.4) +
    geom_point(data = s, aes(years, Nt), shape = 21, fill = OBS_COL, colour = "white", size = PT_IN, stroke = PT_STROKE) +
    panel_letter(letters[i], YB) +
    scale_y_continuous(breaks = seq(0, 2500, 500), expand = c(0, 0)) +
    scale_x_continuous(breaks = seq(1970, 2020, 10), expand = c(0, 0)) +
    coord_cartesian(xlim = xlim, ylim = YB, expand = FALSE) +
    labs(x = "Year", y = "Hake total biomass (1000 t)") + theme_nat(i == 3, i == 2)

  panels[[2 * i]] <- ggplot() + breaks_layer +
    geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.3) +
    geom_ribbon(data = b, aes(Year, ymin = ci_low, ymax = ci_high), fill = BOOT_COL, alpha = 0.25) +
    geom_line(data = b, aes(Year, R_pred), colour = PRED_COL, linewidth = 0.9) +
    geom_point(data = d, aes(year, R), shape = 21, fill = "white", colour = "grey45", size = PT_OUT, stroke = 0.4) +
    geom_point(data = b, aes(Year, R_obs), shape = 21, fill = OBS_COL, colour = "white", size = PT_IN, stroke = PT_STROKE) +
    panel_letter(letters[i + 3], YR) +
    scale_y_continuous(breaks = seq(-0.4, 0.4, 0.2), expand = c(0, 0)) +
    scale_x_continuous(breaks = seq(1970, 2020, 10), expand = c(0, 0)) +
    coord_cartesian(xlim = xlim, ylim = YR, expand = FALSE) +
    labs(x = "Year", y = expression("Reproductive rate, " * italic(R)[t])) + theme_nat(i == 3, i == 2)
}
# legend, drawn once in panel a: top aligned with the panel letter, in the upper-right corner (right after
# the 2001 break line), clear of the observations and of the prediction bands
LEG_X <- 2001.8; LEG_W <- 1.6; LEG_T <- 2.3; LEG_S <- 2.2   # key start, key length (years), text offset, text size
LEG_F <- c(0.975, 0.925, 0.875, 0.825, 0.775, 0.725)     # row positions as fractions of the y range (first row level with the letter)
leg_rows <- function(lim, keys, labs) {
  ys <- fy(lim, LEG_F[seq_along(keys)]); dy <- 0.013 * diff(lim); out <- list()
  for (k in seq_along(keys)) out <- c(out, keys[[k]](LEG_X, ys[k], dy))
  c(out, list(annotate("text", x = LEG_X + LEG_T, y = ys, hjust = 0, vjust = 0.5, size = LEG_S, colour = "grey20", label = labs)))
}
k_obs  <- function(x, y, dy) list(annotate("point", x = x + LEG_W / 2, y = y, shape = 21, fill = OBS_COL, colour = "white", size = PT_IN, stroke = PT_STROKE))
k_oth  <- function(x, y, dy) list(annotate("point", x = x + LEG_W / 2, y = y, shape = 21, fill = "white", colour = "grey45", size = PT_OUT))
k_line <- function(x, y, dy) list(annotate("segment", x = x, xend = x + LEG_W, y = y, yend = y, colour = PRED_COL, linewidth = 0.9))
k_ci   <- function(x, y, dy) list(annotate("rect", xmin = x, xmax = x + LEG_W, ymin = y - dy, ymax = y + dy, fill = PRED_COL, alpha = 0.30))
k_pi   <- function(x, y, dy) list(annotate("rect", xmin = x, xmax = x + LEG_W, ymin = y - dy, ymax = y + dy, fill = PRED_COL, alpha = 0.12))
k_boot <- function(x, y, dy) list(annotate("rect", xmin = x, xmax = x + LEG_W, ymin = y - dy, ymax = y + dy, fill = BOOT_COL, alpha = 0.25))
panels[[1]] <- panels[[1]] + leg_rows(YB, list(k_obs, k_oth, k_line, k_ci, k_pi, k_boot),
  c("Observed, modeled", "Observed, other years", "Predicted", "95% CI (parameters)", "95% PI (+ process)", "95% CI (bootstrap)"))
# The key of panel d was removed on 2026-11-03: one key in panel a, which includes the bootstrap band of d-f.

fig <- (panels[[1]] | panels[[2]]) / (panels[[3]] | panels[[4]]) / (panels[[5]] | panels[[6]])
save_fig(fig, FIG3_NAME, width_cm = 18, height_cm = 20)
