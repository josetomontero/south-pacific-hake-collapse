# Build the analysis dataset: hake and jumbo squid stock-assessment biomass, catches,
# exploitation rate, reproductive rate and the environmental predictors with lags.
#
# Source of data/fisheries/ch_ts_sst.txt (one row per year, 1968-2022):
#   Catchch  hake catch (t), SERNAPESCA, same series used by the IFOP assessment
#   btch     hake total biomass (t), IFOP statistical catch-at-age assessment
#   U        exploitation rate = Catchch / btch
#   Rch      reproductive rate = log(btch[t+1] / btch[t])
#   DGbt     jumbo squid total biomass (thousand t), IFOP assessment (Table 10, 2001-2022)
#   SST0 in ch_ts_sst.txt is not used: SST is taken from data/derived/SSTchileStocks_2022.txt,
#            column CHU (management unit polygon, smoothed); SST1, SST2 = lags
if (!exists("PROJ")) source(file.path(rprojroot::find_root(rprojroot::has_file("R/config.R")), "R", "config.R"))

h <- tsv_read(file.path(PROJ, "data", "fisheries", "ch_ts_sst.txt"))
h <- subset(h, !is.na(year) & year < 2022)

# Re-derive the computed columns so the definitions are explicit (and check them)
stopifnot(max(abs(h$U - h$Catchch / h$btch), na.rm = TRUE) < 1e-8)
R_check <- c(log(h$btch[-1] / h$btch[-nrow(h)]), NA)
stopifnot(max(abs(h$Rch - R_check), na.rm = TRUE) < 1e-8)

sst <- tsv_read(file.path(PROJ, "data", "derived", "SSTchileStocks_2022.txt"))
sst_col <- "CHU"   # HadISST1 cells inside the fishery management unit polygon (00_extract_SST_HadISST.R)

enso <- tsv_read(file.path(PROJ, "data", "derived", "enso34_smoothed.txt"))
oni  <- tsv_read(file.path(PROJ, "data", "derived", "oni_smoothed.txt"))

lagged <- function(x, years, target, k) x[match(target - k, years)]

d <- data.frame(
  year     = h$year,
  catch    = h$Catchch,
  biomass  = h$btch,
  logB     = log(h$btch),
  R        = h$Rch,
  U        = h$U,
  squid_B  = h$DGbt,
  log_squid_B = log(h$DGbt),
  SST0     = lagged(sst[[sst_col]], sst$Years, h$year, 0),
  SST1     = lagged(sst[[sst_col]], sst$Years, h$year, 1),
  SST2     = lagged(sst[[sst_col]], sst$Years, h$year, 2),
  ENSO34   = lagged(enso$enso34, enso$Year, h$year, 0),
  ENSO341  = lagged(enso$enso34, enso$Year, h$year, 1),
  ENSO342  = lagged(enso$enso34, enso$Year, h$year, 2),
  ONI0     = lagged(oni$oni35,   oni$Year,  h$year, 0),
  ONI1     = lagged(oni$oni35,   oni$Year,  h$year, 1)
)

tsv_write(d, file.path(PROJ, "data", "derived", "hake_dataset.txt"))
cat("rows:", nrow(d), " years:", range(d$year), "\n")
