# Smooth the climate indices used as predictors (smooth.spline, spar = 0.35).
#
# ENSO 3.4: data/climate/enso34.txt holds monthly Nino 3.4 SST (HadISST1, 1870-2021;
#   columns 1..12 = months, Sum/Aut/Win/Spring = seasonal means, ENSO34 = annual mean).
#   The predictor is the annual mean (column ENSO34), smoothed over the full 1870-2021
#   record and then restricted to 1960-2021.
# ONI:      data/climate/ONI.txt, annual mean of the monthly Oceanic Nino Index, smoothed
#   over 1950-2022 and restricted to 1960 onwards. ONI is a column of the analysis dataset and
#   is not used in the models.
# SST:      produced by 00_extract_SST_HadISST.R (column CHU of SSTchileStocks_2022.txt,
#   HadISST1 cells inside the fishery management unit polygon).
if (!exists("PROJ")) source(file.path(rprojroot::find_root(rprojroot::has_file("R/config.R")), "R", "config.R"))

SPAR <- 0.35

# ENSO 3.4 (HadISST1 Nino 3.4 SST), annual mean, smoothed over the full 1870-2021 record and
# restricted to 1960-2021 afterwards.
enso_raw <- tsv_read(file.path(PROJ, "data", "climate", "enso34.txt"))
names(enso_raw)[2:13] <- month.abb
stopifnot(max(abs(enso_raw$ENSO34 - rowMeans(enso_raw[, month.abb]))) < 1e-6)   # ENSO34 = annual mean
enso <- data.frame(Year = enso_raw$Year, enso34_raw = enso_raw$ENSO34, enso34 = smooth.spline(enso_raw$ENSO34, spar = SPAR)$y)
enso <- subset(enso, Year >= 1960 & Year <= 2021)
tsv_write(enso, file.path(PROJ, "data", "derived", "enso34_smoothed.txt"))

# ONI (annual mean of monthly values)
oni_raw <- tsv_read(file.path(PROJ, "data", "climate", "ONI.txt"))
oni_fit <- smooth.spline(oni_raw$mean, spar = SPAR)
oni <- data.frame(Year = oni_raw$year, oni_raw = oni_raw$mean, oni35 = oni_fit$y)
oni <- subset(oni, Year >= 1960)
tsv_write(oni, file.path(PROJ, "data", "derived", "oni_smoothed.txt"))


# Figure of the raw and smoothed series (Appendix S1: Figure S2)
ragg::agg_png(file.path(RES, "figs", "climate_smoothing.png"), width = 7, height = 7, units = "in", res = 600)
par(mfrow = c(2, 1), mar = c(4, 4.5, 1, 1), bty = "l")
plot(enso$Year, enso$enso34_raw, type = "b", pch = 19, cex = 0.6, col = "grey40", ylim = range(enso[, -1]),
     xlab = "Year", ylab = "Nino 3.4 SST (C)")
lines(enso$Year, enso$enso34, lwd = 2.5, col = "red")
legend("topleft", c("HadISST1 annual mean", "smooth.spline, spar = 0.35"),
       col = c("grey40", "red"), lwd = c(1, 2.5), pch = c(19, NA), bty = "n", cex = 0.9)
sst <- tsv_read(file.path(PROJ, "data", "derived", "SSTchileStocks_raw_2022.txt"))
sst_s <- tsv_read(file.path(PROJ, "data", "derived", "SSTchileStocks_2022.txt"))
plot(sst$Years, sst$CHU, type = "b", pch = 19, cex = 0.6, col = "grey40",
     xlab = "Year", ylab = "SST hake management unit (C)")
lines(sst_s$Years, sst_s$CHU, lwd = 2.5, col = "red")
legend("topleft", c("HadISST1 annual mean, fishery management unit", "smooth.spline, spar = 0.35"),
       col = c("grey40", "red"), lwd = c(1, 2.5), pch = c(19, NA), bty = "n", cex = 0.9)
dev.off()
