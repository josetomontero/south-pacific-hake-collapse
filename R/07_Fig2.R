# Figure 2. (a) hake total biomass with structural breaks, (b) exploitation-rate anomalies,
# (c) exploitation rate vs biomass, (d) jumbo squid biomass, (e) ENSO 3.4, (f) SST.
# Anomalies in (b): U is linearly detrended (lm(U ~ year), 1968-2019) and the detrended
# values are centred on their mean. Colours in (a) and (c) follow the segments between the
# breaks (1974, 1983, 1992, 2001). Biomass axes start at zero and the squid biomass in (d)
# uses the same scale as the hake biomass in (a).
if (!exists("PROJ")) source(file.path(rprojroot::find_root(rprojroot::has_file("R/config.R")), "R", "config.R"))
source(file.path(PROJ, "R", "lib", "fig_theme.R"))

d <- tsv_read(file.path(PROJ, "data", "derived", "hake_dataset.txt"))
d <- subset(d, year <= 2020)
d$biomass_kt <- d$biomass / 1000
d$segment <- seg_label(d$year)

ua <- subset(d, year <= 2019)
trend <- lm(U ~ year, data = ua)
ua$anomaly <- resid(trend) - mean(resid(trend))
tsv_write(ua[, c("year", "U", "anomaly")], file.path(RES, "tables", "exploitation_rate_anomalies.txt"))

xlim <- c(1965, 2021)
PT_SIZE <- 2.0; PT_STROKE <- 0.35   # points as in Figure 3 (size 2.0), with a thin white rim so that neighbours do not merge
# Panel letters as in Figure 3: plain "a)", size 3.4 mm, inside the panel with the top of the
# letter flush with the top of the y axis
letter_theme <- theme(plot.tag = element_text(size = 3.4 * .pt, face = "plain", hjust = 0, vjust = 1.05),
                      plot.tag.location = "panel", plot.tag.position = c(0.008, 1)) +
  theme_axes_fig3   # same axis line and tick length as Figure 3
breaks_layer <- geom_vline(xintercept = BREAKS, linetype = "dotted", colour = "grey55", linewidth = 0.4)

pa <- ggplot(d, aes(year, biomass_kt)) + breaks_layer +
  geom_path(colour = "grey60", linewidth = 0.4) +
  geom_point(aes(fill = segment), shape = 21, colour = "white", stroke = PT_STROKE, size = PT_SIZE) +
  scale_fill_manual(values = SEG_COLS) +
  scale_y_continuous(limits = c(0, 1500), breaks = seq(0, 1500, 500), expand = c(0, 0)) +
  scale_x_continuous(limits = xlim, breaks = seq(1970, 2020, 10)) +
  labs(x = NULL, y = "Hake total biomass (1000 t)", tag = "a)") +
  theme_pub() + letter_theme +
  theme(legend.position = c(0.98, 0.97), legend.justification = c(1, 1))   # upper right: clear of the series and the break lines

pb <- ggplot(ua, aes(year, anomaly)) + breaks_layer +
  annotate("rect", xmin = 1974, xmax = 2001, ymin = -Inf, ymax = 0, fill = "grey85", alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "grey30", linewidth = 0.4) +
  geom_line(colour = PRED_COL, linewidth = 0.8) +
  scale_y_continuous(limits = c(-0.1, 0.3), breaks = seq(-0.1, 0.2, 0.1), expand = c(0, 0)) +   # headroom for the letter
  scale_x_continuous(limits = xlim, breaks = seq(1970, 2020, 10)) +
  labs(x = NULL, y = "Exploitation rate anomaly", tag = "b)") + theme_pub() + letter_theme

pc <- ggplot(d, aes(U, biomass_kt)) +
  geom_path(colour = "grey60", linewidth = 0.4) +
  geom_point(aes(fill = segment), shape = 21, colour = "white", stroke = PT_STROKE, size = PT_SIZE, show.legend = FALSE) +
  scale_fill_manual(values = SEG_COLS) +
  scale_y_continuous(limits = c(0, 1500), breaks = seq(0, 1500, 500), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(0.05, 0.2, 0.35)) +
  labs(x = "Exploitation rate (U)", y = "Hake total biomass (1000 t)", tag = "c)") + theme_pub() + letter_theme

pd <- ggplot(subset(d, !is.na(squid_B)), aes(year, squid_B)) +
  geom_path(colour = "grey60", linewidth = 0.4) + geom_point(shape = 21, fill = "black", colour = "white", stroke = PT_STROKE, size = PT_SIZE) +
  scale_y_continuous(limits = c(0, 1500), breaks = seq(0, 1500, 500), expand = c(0, 0)) +   # same scale as the hake biomass in (a)
  scale_x_continuous(limits = xlim, breaks = seq(1970, 2020, 10)) +
  labs(x = NULL, y = "Jumbo squid total biomass (1000 t)", tag = "d)") + theme_pub() + letter_theme

pe <- ggplot(d, aes(year, ENSO34)) + geom_line(colour = "#0072B2", linewidth = 0.8) +
  scale_x_continuous(limits = xlim, breaks = seq(1970, 2020, 10)) +
  labs(x = NULL, y = expression("Niño 3.4 SST, smoothed (" * degree * "C)"), tag = "e)") + theme_pub() + letter_theme

pf <- ggplot(d, aes(year, SST0)) + geom_line(colour = "#D55E00", linewidth = 0.8) +
  scale_x_continuous(limits = xlim, breaks = seq(1970, 2020, 10)) +
  labs(x = "Year", y = expression("SST fishing area, smoothed (" * degree * "C)"), tag = "f)") + theme_pub() + letter_theme

fig <- (pa | pd) / (pb | pe) / (pc | pf)
save_fig(figure_label(fig, 2), "Fig_2", width_cm = 18, height_cm = 20)
