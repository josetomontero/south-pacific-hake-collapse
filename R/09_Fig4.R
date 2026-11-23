# Figure 4. Phase portrait (Berryman 1999, Royama 1992) of the two regimes on one set of axes.
# Everything comes from the models of the paper: the best model of each period (results/models/best_P2.rds and
# best_P3.rds, steps 04 and 05) and model 1 of each period, without exogenous variables (H0, models_1975_2000.txt and models_2001_2020.txt).
#   - observed (B_t, R_t) of 1975-2019 joined in time order; points colored by regime (regime I in the blue of
#     Figure 2, regime II black), the collapse years 2001-2003 in their own color;
#   - complete H0 reproduction curve of each regime, R = Rm - c exp(X), drawn over the whole axis;
#   - one partial curve per regime from its best model, dashed. Regime I, model 7, with ENSO at the upper quartile
#     of 1975-2000 and SST at its period mean. Regime II, model 5, with U at its observed maximum. Both effects are
#     lateral (on c). ENSO is not set at its maximum because the equilibrium would then lie far outside the observed
#     biomass. Every equilibrium of the figure lies within the biomass observed in its
#     regime (checked below);
#   - equilibria X* = ln(Rm / c_eff) as open circles on R = 0; the collapse shift between the two H0 equilibria as a
#     curved dashed arrow; the effect of each variable as a solid arrow between equilibria, labeled with the variable.
# Outputs: results/figs/Fig_4.png, results/tables/Fig_4_equilibria.txt
if (!exists("PROJ")) source(file.path(rprojroot::find_root(rprojroot::has_file("R/config.R")), "R", "config.R"))
source(file.path(PROJ, "R", "lib", "fig_theme.R"))
COL <- c("Regime I" = "#0072B2",    # the blue of Figure 2
         "Regime II" = "#000000",   # black, as in Figures 2 and 3
         "Shift" = "#9E2A2B")       # collapse years 2001-2003
# H0: model 1 of each period (no exogenous variable) in the tables of steps 04 and 05
h0 <- function(file, rm) {
  t <- tsv_read(file.path(RES, "tables", file)); stopifnot(t$ID[1] == 1, t$N.param[1] == 1, grepl(sprintf("^R ~ %.2f - \\(?c\\)? \\* exp\\(Xt\\)$", rm), t$Models[1]))
  t$c[1]
}

# ---- equilibria -------------------------------------------------------------------------------------------
bI <- readRDS(file.path(RES, "models", "best_P2.rds")); bII <- readRDS(file.path(RES, "models", "best_P3.rds"))
RmI <- PERIODS$P2$rm; RmII <- PERIODS$P3$rm
c0I <- h0("models_1975_2000.txt", RmI); c0II <- h0("models_2001_2020.txt", RmII)
dI <- bI$data; dII <- bII$data
stopifnot(bI$ID == 7, bII$ID == 5)
cfI <- coef(bI$model); ensoHigh <- unname(quantile(dI$V3, 0.75))   # regime I: model 7 with ENSO at its upper quartile, SST at its mean
ceI <- cfI[["c"]] + cfI[["e"]] * mean(dI$V2) + cfI[["d"]] * ensoHigh; stopifnot(ceI > 0)
cfII <- coef(bII$model); Umax <- max(dII$V1)                      # regime II: model 5 at U max (lateral)
ceII <- cfII[["c"]] + cfII[["e"]] * Umax; stopifnot(ceII > 0)

eq <- data.frame(key = c("I_H0", "II_H0", "I_partial", "II_partial"), regime = c("I", "II", "I", "II"),
                 model = c(1, 1, bI$ID, bII$ID), form = c("H0", "H0", "lateral", "lateral"),
                 variable = c(NA, NA, "ENSO", "U"), value = c(NA, NA, round(ensoHigh, 3), round(Umax, 3)),
                 fixed = c(NA, NA, sprintf("SST = %.2f", mean(dI$V2)), NA),
                 Rm = c(RmI, RmII, RmI, RmII), c_eff = signif(c(c0I, c0II, ceI, ceII), 5))
eq$Xstar <- round(log(eq$Rm / eq$c_eff), 3); eq$Bstar_kt <- round(eq$Rm / eq$c_eff / 1000, 1)
RNG <- list(I = range(exp(dI$Xt)) / 1000, II = range(exp(dII$Xt)) / 1000)          # biomass observed in each regime (thousand t)
eq$B_min_kt <- round(sapply(eq$regime, function(r) RNG[[r]][1]), 1); eq$B_max_kt <- round(sapply(eq$regime, function(r) RNG[[r]][2]), 1)
stopifnot(eq$Bstar_kt >= eq$B_min_kt, eq$Bstar_kt <= eq$B_max_kt)                    # every equilibrium within the observed biomass of its regime
tsv_write(eq, file.path(RES, "tables", "Fig_4_equilibria.txt"))
print(eq[, c("key", "model", "variable", "value", "Rm", "c_eff", "Xstar", "Bstar_kt", "B_min_kt", "B_max_kt")], row.names = FALSE)
g <- function(k, v) eq[[v]][eq$key == k]

# ---- data for the drawing --------------------------------------------------------------------------------
pts <- rbind(data.frame(regime = "Regime I", Year = dI$Year, B_kt = exp(dI$Xt) / 1000, R = dI$R),
             data.frame(regime = "Regime II", Year = dII$Year, B_kt = exp(dII$Xt) / 1000, R = dII$R))
pts <- pts[order(pts$Year), ]; stopifnot(!any(duplicated(pts$Year)), all(diff(pts$Year) == 1))
pts$group <- ifelse(pts$Year %in% 2001:2003, "Shift", pts$regime)
XLIM <- c(120, 3200); xg <- seq(log(XLIM[1] * 1e3), log(XLIM[2] * 1e3), length.out = 800)
curve <- function(key) data.frame(key = key, B_kt = pmin(pmax(exp(xg) / 1000, XLIM[1]), XLIM[2]), R = g(key, "Rm") - g(key, "c_eff") * exp(xg))
cu <- do.call(rbind, lapply(eq$key, curve))
SH <- c(I_H0 = COL[["Regime I"]], I_partial = COL[["Regime I"]],      # partial curves in the colour of their regime,
        II_H0 = COL[["Regime II"]], II_partial = COL[["Regime II"]])   # told apart from H0 by the dashed line
LT_PART <- "42"
PT <- COL[c("Regime I", "Regime II", "Shift")]
XI <- g("I_H0", "Bstar_kt"); XII <- g("II_H0", "Bstar_kt"); XIe <- g("I_partial", "Bstar_kt"); XIIe <- g("II_partial", "Bstar_kt")
EDGE <- 1.035            # circle radius in the log-biomass axis (arrows start and end at the circle edge)
YLIM <- c(-0.5, 0.455)   # the top is the top of the key
# 16 x 12 cm. Text enlarged as far as the layout allows
# without overlaps in the PDF (cairo fonts): key and annotations 7.8 pt (were 6.8), X* labels 9 pt, axis text
# 10 pt and axis titles 11 pt (theme_pub base 11). The long key label of regime I fixes the limit.
FIG_W <- 16; FIG_H <- 12; TXT <- 7.8 / .pt; BASE <- 11

# Semi-bold label: the PDF device (cairo) has no medium weight of Helvetica, so the plain text is drawn
# twice with a horizontal offset of 0.03 mm, which thickens the strokes to about half way to bold.
semibold <- function(x, y, label, colour, size = 7 / .pt, dx = 1.0008)
  list(annotate("text", x = x, y = y, label = label, size = size, colour = colour),
       annotate("text", x = x * dx, y = y, label = label, size = size, colour = colour))

# White box of the size of an arrow label, without border. It is drawn just above the grey trajectory and below
# the curves, arrows and points, so that the trajectory does not run through the text and nothing else is masked.
knockout <- function(x, y, label, size = 7 / .pt)
  annotate("label", x = x, y = y, label = label, size = size, colour = "white", fill = "white", border.colour = NA,
           label.padding = unit(0.3, "pt"), label.r = unit(0, "pt"))
LAB_I <- list(x = XIe * 1.27, y = 0.024, lab = sprintf("ENSO %.1f °C", ensoHigh))   # at the right of the equilibrium, clear of the observed points
LAB_II <- list(x = sqrt(XII * XIIe), y = 0.02, lab = sprintf("U %.2f", Umax))            # centred over the arrow

p <- ggplot() +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey45") +
  geom_path(data = pts, aes(B_kt, R), colour = "grey72", linewidth = 0.3) +
  knockout(LAB_I$x, LAB_I$y, LAB_I$lab) + knockout(LAB_II$x, LAB_II$y, LAB_II$lab) +
  geom_line(data = cu[cu$key %in% c("I_H0", "II_H0"), ], aes(B_kt, R, colour = key), linewidth = 1.0) +
  geom_line(data = cu[cu$key %in% c("I_partial", "II_partial"), ], aes(B_kt, R, colour = key), linewidth = 0.7, linetype = LT_PART) +
  # collapse: curved dashed arrow below R = 0 from the edge of X*_I to the edge of X*_II
  annotate("curve", x = XI / EDGE, xend = XII * EDGE, y = -0.012, yend = -0.012, curvature = -0.28, colour = "grey40",   # grey, so that near X*_II the only black dashed line is the partial curve
           linewidth = 0.5, linetype = "22", arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
  annotate("text", x = sqrt(XI * XII), y = -0.108, label = "collapse 2001–2003", size = TXT, colour = "grey30") +
  # effect arrows on R = 0, circle edge to circle edge, labelled with the variable
  annotate("segment", x = XI * EDGE, xend = XIe / EDGE, y = 0, yend = 0, colour = SH[["I_partial"]], linewidth = 0.6,
           arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
  annotate("segment", x = XII / EDGE, xend = XIIe * EDGE, y = 0, yend = 0, colour = SH[["II_partial"]], linewidth = 0.6,
           arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
  semibold(LAB_I$x, LAB_I$y, LAB_I$lab, SH[["I_partial"]]) +
  semibold(LAB_II$x, LAB_II$y, LAB_II$lab, SH[["II_partial"]]) +
  # X*_I and X*_II: dotted drop line from the circle to a label just above the x axis
  annotate("segment", x = c(XI, XII), xend = c(XI, XII), y = -0.012, yend = YLIM[1] + 0.062, colour = SH[c("I_H0", "II_H0")],
           linewidth = 0.4, linetype = "dotted") +
  annotate("text", x = c(XI, XII), y = YLIM[1] + 0.052, label = c("italic(X)[I]^'*'", "italic(X)[II]^'*'"), parse = TRUE, size = TXT * 1.15,
           colour = SH[c("I_H0", "II_H0")], vjust = 1) +
  geom_point(data = eq, aes(Bstar_kt, 0, colour = key), shape = 21, fill = "white", size = 2.8, stroke = 0.9, show.legend = FALSE) +
  geom_point(data = pts, aes(B_kt, R, fill = group), shape = 21, colour = "white", size = 2.3, stroke = 0.3, show.legend = FALSE) +
  # the equilibrium under high ENSO lies among the observed points of 1993-1996, so its circle is drawn again above them
  geom_point(data = eq[eq$key == "I_partial", ], aes(Bstar_kt, 0, colour = key), shape = 21, fill = "white", size = 2.8, stroke = 0.9, show.legend = FALSE) +
  scale_colour_manual(values = SH, guide = "none") + scale_fill_manual(values = PT, guide = "none") +
  scale_x_log10(limits = XLIM, breaks = c(200, 500, 1000, 2000), labels = scales::label_comma(), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(-0.5, 0.25, 0.25)) +
  coord_cartesian(ylim = YLIM, clip = "on", expand = FALSE) +
  labs(x = expression("Hake total biomass, " * italic(B[t]) * " (1000 t, log scale)"), y = expression("Reproductive rate, " * italic(R[t]))) +
  theme_pub(base_size = BASE) + theme_axes_fig3 + theme(legend.position = "none")

# key drawn in data coordinates: entries = list(list(type = "line"/"point", col, lw, lt, lab)), one row per entry (lt optional, solid by default)
key_rows <- function(p, x0, y0, entries, dy, seglen = 1.3, size = TXT) {
  for (i in seq_along(entries)) {
    e <- entries[[i]]; y <- y0 - (i - 1) * dy
    if (e$type == "line") p <- p + annotate("segment", x = x0, xend = x0 * seglen, y = y, yend = y, colour = e$col, linewidth = e$lw,
                                            linetype = if (is.null(e$lt)) "solid" else e$lt)
    else p <- p + annotate("point", x = x0 * sqrt(seglen), y = y, shape = 21, fill = e$col, colour = "white", size = 2.4, stroke = 0.3)
    p <- p + annotate("text", x = x0 * seglen * 1.06, y = y, label = e$lab, hjust = 0, vjust = 0.5, size = size, colour = "grey15")
  }
  p
}
L <- list(list(type = "line", col = SH[["I_H0"]], lw = 1.0, lab = "Regime I 1975–2000, H0"),
          list(type = "line", col = SH[["I_partial"]], lw = 0.7, lt = LT_PART, lab = "Regime I, ENSO at its upper quartile"),
          list(type = "line", col = SH[["II_H0"]], lw = 1.0, lab = "Regime II 2001–2020, H0"),
          list(type = "line", col = SH[["II_partial"]], lw = 0.7, lt = LT_PART, lab = "Regime II, U at its maximum"))
P <- list(list(type = "point", col = PT[["Regime I"]], lab = "Regime I"), list(type = "point", col = PT[["Regime II"]], lab = "Regime II"),
          list(type = "point", col = PT[["Shift"]], lab = "Collapse 2001–2003"))
p <- key_rows(p, 130, 0.44, P[1], 0.04); p <- key_rows(p, 260, 0.44, P[2], 0.04); p <- key_rows(p, 520, 0.44, P[3], 0.04)
p <- key_rows(p, 130, 0.40, L[1:2], 0.04); p <- key_rows(p, 950, 0.40, L[3:4], 0.04)   # right column at 950 kt: the long left label ends near 800 kt in the PDF
# "Regime II 2001–2020" is the fitting period of the regime II models (PERIODS$P3); R_t ends in 2019 because it needs the biomass of the following year

save_fig(p, "Fig_4", width_cm = FIG_W, height_cm = FIG_H)
