# Structural breaks in the hake total-biomass series (strucchange, Zeileis et al. 2002).
# The analysis uses the complete IFOP assessment series 1960-2020
# (data/fisheries/chile_stocks_admb_outs_AES2.txt, Specie == "Chilean hake", column bt).
# The BIC-optimal segmentation gives breaks at 1974, 1983, 1992 and 2001.
if (!exists("PROJ")) source(file.path(rprojroot::find_root(rprojroot::has_file("R/config.R")), "R", "config.R"))
library(strucchange)

dem  <- tsv_read(file.path(PROJ, "data", "fisheries", "chile_stocks_admb_outs_AES2.txt"))
hake <- subset(dem, Specie == "Chilean hake" & !is.na(bt))
tsv_write(hake[, c("Year", "bt", "ssb", "Fi")],
          file.path(PROJ, "data", "derived", "hake_biomass_IFOP_1960_2020.txt"))

bt  <- ts(hake$bt, start = min(hake$Year))
bp  <- breakpoints(bt ~ 1)                 # all segmentations, h = 0.15 (default)
fs  <- Fstats(bt ~ 1)                      # supF test for a single break
sfs <- sctest(fs)

opt <- summary(bp)
tab <- data.frame(n_breaks = 0:(nrow(opt$breakpoints)),
                  RSS = opt$RSS["RSS", ], BIC = opt$RSS["BIC", ])
tab$break_years <- c("", sapply(1:nrow(opt$breakpoints), function(m)
  paste(hake$Year[na.omit(opt$breakpoints[m, ])], collapse = ", ")))
tsv_write(tab, file.path(RES, "tables", "breakpoints_hake_biomass.txt"))

best <- hake$Year[bp$breakpoints]
cat("BIC-optimal breaks:", best, "\n")
cat(sprintf("supF test, single break: F = %.2f, p = %.3g, break = %d\n",
            sfs$statistic, sfs$p.value, hake$Year[breakpoints(fs)$breakpoints]))
stopifnot(all(best == BREAKS))

ragg::agg_png(file.path(RES, "figs", "breakpoints.png"), width = 7, height = 7, units = "in", res = 600)
par(mfrow = c(2, 1), mar = c(4, 4.5, 1, 1), bty = "l")
plot(bp, main = "")
plot(hake$Year, hake$bt / 1000, type = "b", pch = 19, cex = 0.7,
     xlab = "Year", ylab = "Hake total biomass (1000 t)")
abline(v = best, lty = 3, col = "darkgrey")
lines(hake$Year, fitted(bp) / 1000, col = "red", lwd = 2)
dev.off()
