# El Nino and La Nina episodes, 1968-2020, against the periods of the population models.
# Descriptive only. Nothing here
# enters the population models, whose ENSO predictor is the smoothed annual mean of Nino 3.4 SST
# (01_smooth_climate_indices.R).
#
# Data: data/climate/ONI_NOAA_CPC_seasonal_20261111.txt, the Oceanic Nino Index of NOAA CPC
#   (3-month running mean of ERSSTv5 Nino 3.4 SST anomalies, columns SEAS, YR, TOTAL, ANOM),
#   downloaded on 2026-11-11 from https://www.cpc.ncep.noaa.gov/data/indices/oni.ascii.txt.
#   It is a different vintage from data/climate/ONI.txt (not touched, used by step 01).
# Episode rule (NOAA CPC): ONI rounded to one decimal at or above +0.5 C (El Nino) or at or
#   below -0.5 C (La Nina) for at least five consecutive overlapping seasons. The rule is applied
#   mechanically to the file. Episode counts depend on the vintage of the index and on this rule
#   (the 1979-80 and 2019-20 warm spells and the 2000-01 and 2016-17 cold spells are borderline).
#
# Outputs: results/figs/ENSO_episodes.png
#          results/tables/enso_episodes_1968_2020.txt   (one row per episode)
#          results/tables/enso_episodes_by_period.txt   (onsets and share of seasons per period)
if (!exists("PROJ")) source(file.path(rprojroot::find_root(rprojroot::has_file("R/config.R")), "R", "config.R"))
source(file.path(PROJ, "R", "lib", "fig_theme.R"))

SEAS <- c("DJF", "JFM", "FMA", "MAM", "AMJ", "MJJ", "JJA", "JAS", "ASO", "SON", "OND", "NDJ")
THR <- 0.5; MIN_SEASONS <- 5; Y0 <- 1968; Y1 <- 2020
EN_COL <- "#C0392B"; LN_COL <- "#0072B2"

oni <- read.table(file.path(PROJ, "data", "climate", "ONI_NOAA_CPC_seasonal_20261111.txt"), header = TRUE, stringsAsFactors = FALSE)
oni$mo <- match(oni$SEAS, SEAS)                       # centre month of the season (DJF = January)
stopifnot(!anyNA(oni$mo), all(diff(oni$YR * 12 + oni$mo) == 1))
oni$t <- oni$YR + (oni$mo - 0.5) / 12
oni$a1 <- round(oni$ANOM, 1)

# Runs of at least MIN_SEASONS consecutive seasons beyond the threshold
episodes <- function(flag, phase) {
  r <- rle(flag); e <- cumsum(r$lengths); s <- e - r$lengths + 1; k <- which(r$values & r$lengths >= MIN_SEASONS)
  data.frame(phase = phase, i0 = s[k], i1 = e[k], seasons = r$lengths[k])
}
ep <- rbind(episodes(oni$a1 >= THR, "El Nino"), episodes(oni$a1 <= -THR, "La Nina"))
ep$onset_year <- oni$YR[ep$i0]; ep$onset_season <- oni$SEAS[ep$i0]
ep$end_year <- oni$YR[ep$i1];   ep$end_season <- oni$SEAS[ep$i1]
ep$peak <- mapply(function(a, b, ph) if (ph == "El Nino") max(oni$a1[a:b]) else min(oni$a1[a:b]), ep$i0, ep$i1, ep$phase)
ep$x0 <- oni$t[ep$i0] - 0.5 / 12; ep$x1 <- oni$t[ep$i1] + 0.5 / 12
oni$state <- "neutral"
for (i in seq_len(nrow(ep))) oni$state[ep$i0[i]:ep$i1[i]] <- ep$phase[i]

ep <- ep[order(ep$i0), ]
ep_w <- subset(ep, onset_year >= Y0 & onset_year <= Y1)
tsv_write(ep_w[, c("phase", "onset_year", "onset_season", "end_year", "end_season", "seasons", "peak")],
          file.path(RES, "tables", "enso_episodes_1968_2020.txt"))

# Onsets and share of seasons by period (periods of the models, then decades)
PER <- list("1968-1974" = 1968:1974, "1975-2000" = 1975:2000, "2001-2020" = 2001:2020,
            "1968-1979" = 1968:1979, "1980-1989" = 1980:1989, "1990-1999" = 1990:1999,
            "2000-2009" = 2000:2009, "2010-2020" = 2010:2020)
by_period <- do.call(rbind, lapply(names(PER), function(nm) {
  y <- PER[[nm]]; k <- oni$YR %in% y
  on <- function(ph) sum(ep$phase == ph & ep$onset_year %in% y)
  data.frame(period = nm, years = length(y),
             EN_onsets = on("El Nino"), LN_onsets = on("La Nina"),
             EN_onsets_per_decade = round(10 * on("El Nino") / length(y), 1),
             LN_onsets_per_decade = round(10 * on("La Nina") / length(y), 1),
             pct_seasons_EN = round(100 * mean(oni$state[k] == "El Nino")),
             pct_seasons_LN = round(100 * mean(oni$state[k] == "La Nina")),
             mean_ONI = round(mean(oni$ANOM[k]), 2),
             EN_peak_ge_1.5 = sum(ep$phase == "El Nino" & ep$onset_year %in% y & ep$peak >= 1.5),
             LN_peak_le_m1.5 = sum(ep$phase == "La Nina" & ep$onset_year %in% y & ep$peak <= -1.5))
}))
tsv_write(by_period, file.path(RES, "tables", "enso_episodes_by_period.txt"))

# Exact Poisson comparison of the onset rates of the two modelled regimes (few events, low power)
rt <- function(col) poisson.test(c(by_period[by_period$period == "1975-2000", col], by_period[by_period$period == "2001-2020", col]),
                                 c(26, 20))$p.value
cat(sprintf("Onset rate 1975-2000 vs 2001-2020, exact Poisson test: El Nino p = %.2f, La Nina p = %.2f\n", rt("EN_onsets"), rt("LN_onsets")))
cat(sprintf("Episodes with onset in %d-%d: %d El Nino, %d La Nina\n", Y0, Y1, sum(ep_w$phase == "El Nino"), sum(ep_w$phase == "La Nina")))

# Figure
d <- subset(oni, YR >= Y0 & YR <= Y1)
XL <- c(Y0, Y1 + 1); YL <- c(-2.3, 3.1)
rects <- subset(ep, x1 > XL[1] & x0 < XL[2]); rects$x0 <- pmax(rects$x0, XL[1]); rects$x1 <- pmin(rects$x1, XL[2])
p <- ggplot() +
  geom_rect(data = rects, aes(xmin = x0, xmax = x1, ymin = YL[1], ymax = 2.75, fill = phase), alpha = 0.18) +
  geom_hline(yintercept = c(-THR, THR), linewidth = 0.25, colour = "grey55") +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey35") +
  geom_vline(xintercept = c(1975, 2001), linetype = "dotted", linewidth = 0.4) +
  geom_line(data = d, aes(t, ANOM), linewidth = 0.45, colour = "black") +
  annotate("text", x = c(1971.5, 1988, 2011), y = 2.98, size = 8 / .pt,
           label = c("1968–1974", "1975–2000", "2001–2020")) +
  scale_fill_manual(values = c("El Nino" = EN_COL, "La Nina" = LN_COL), labels = c("El Niño", "La Niña")) +
  scale_x_continuous(breaks = seq(1970, 2020, 5), limits = XL, expand = c(0, 0)) +
  scale_y_continuous(breaks = -2:3, limits = YL, expand = c(0, 0)) +
  labs(x = "Year", y = "Oceanic Niño Index (°C)") +
  theme_pub(9) + theme_axes_fig3 +
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(-6, 0, 0, 0)) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.35)))
save_fig(p, "ENSO_episodes", 18, 7.5)
