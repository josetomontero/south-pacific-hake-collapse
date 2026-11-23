# Extract yearly HadISST SST (Rayner et al. 2003) inside each stock's fishing extent
# and smooth with smooth.spline(spar = 0.35). Hake uses extent CH (-74,-71 / -39,-29.5).
# Input : HadISST_sst.nc (not in repo, ~485 MB). Download: https://www.metoffice.gov.uk/hadobs/hadisst/
#         Set its location in R/config.R (HADISST_NC) or the env var HADISST_NC.
# Output: data/climate/SSTchileStocks_2022.txt (smoothed series, all stocks; column CH = hake)
#         data/climate/SSTchileStocks_raw_2022.txt (unsmoothed yearly means)
if (!exists("PROJ")) source(file.path(rprojroot::find_root(rprojroot::has_file("R/config.R")), "R", "config.R"))

###############################################################################################################
# Chile Sea Surface SST Area extraction
library(raster)
library(xts)
library(caTools)
library(fields)

startYear <- 1930                                            # start of the period
endYear   <- 2022                                            # end of the period
subp      <- '1951-01-01/1980-12-01'                         # period for the climatology calculatio

# load HadliSST .nc file
sst    <- brick(HADISST_NC)
Date   <- substr(names(sst),2,11)                            # subset the desired dates
Date   <- gsub('\\.', '\\-', Date)                           # give format to the date
Date   <- as.Date(Date)                                      # date formating as.Date
dstart <- paste(startYear,'01','16',sep='-')                 # Create the start date for the period
dstart <- grep(dstart, Date)                                 # Create the start date for the period
dend   <- paste(endYear,'12','16',sep='-')                   # create the end date for the period
dend   <- grep(dend, Date)                                   # create the end date for the period
sst    <- subset(sst, dstart:dend)                           # subset the SST for the target period
Date   <- Date[dstart:dend]                                  # subset the dates for the target period


# Pelagic fisheries distribution
ACPC <- extent(-77,-70, -24, -16) # original -76
ACN  <- extent(-77,-70, -32, -25) # original -76
ACS  <- extent(-77,-70, -41, -33) # original -76
CS   <- extent(-77,-70, -41, -33) # original -76
JM   <- extent(-77,-70, -42, -16) # original -76

# Demersal Fisheries distribution (North and South)
PDA  <- extent(-77,-69, -57, -41)
PDAN <- extent(-77,-69, -47, -41) 
PDAS <- extent(-77,-70, -57, -47)

# Common hake fisheries distribution
CH   <- extent(-74,-71, -39, -29.5)

# Crustacean fisheries distribution
CNN <- extent(-77,-69, -32, -21.5) # original -75
CNS <- extent(-77,-69, -38.5, -32)
LAN <- extent(-77,-69, -30.5, -26.05)
LAS <- extent(-77,-69, -38.8, -30.05)
LCN <- extent(-77,-69, -30.5, -26.05)
LCS <- extent(-77,-69, -38.8, -26.05)



# function to keep the list of extents with names as created
LIST <- function(...) {
  nms <- sapply(as.list(substitute(list(...))), deparse)[-1]
  setNames(list(...), nms)
}

# make a list with all the extent for cropping the SST values
ChiFishExt <- LIST(ACPC, ACN, ACS, CS, JM,       # pelagic 
                   CNN, CNS, LAN, LAS, LCN, LCS, # crustaceans
                   PDA, PDAN, PDAS, CH)          # demersal 

# make a list with all the extent for cropping the SST values
ChiFishNa <- c("ACPC", "ACN", "ACS", "CS", "JM",         # Pelagic 
               "CNN" , "CNS", "LAN", "LAS","LCN", "LCS", # Crustaceans
               "PDA" , "PDAN", "PDAS", "CH")             # Demersal 

#-------------------------------------------------------------------------------
# CROP the SST IN A LOOP for each extension
#-------------------------------------------------------------------------------
# Parallel processing libraries
library(doParallel)
library(foreach)

# register the cores to use
registerDoParallel(cores=max(1, parallel::detectCores() - 1))

# for each loop to crop the Stocks SST layers
sstcrops <- foreach(i = 1:length(ChiFishExt)) %dopar% {
  crop(sst, ChiFishExt[[i]])
}

# stopCluster()

# rename the list of SST extents with the Stock area name
names(sstcrops) <- names(ChiFishExt)

# Hake fishery management unit (Tascheri 2022): HadISST cells whose centre falls inside the
# polygon built by 00_study_area_polygons.R. Series name CHU (canonical SST of the analysis).
unit <- sf::as_Spatial(sf::st_read(file.path(PROJ, "data", "derived", "fishery_unit_polygon.gpkg"), quiet = TRUE))
sstcrops$CHU <- mask(crop(sst, extent(-77, -69, -43, -28)), unit)
cat("HadISST cells inside the management unit:", sum(!is.na(getValues(sstcrops$CHU[[1]]))), "\n")
#-------------------------------------------------------------------------------
# Calculate layer and cell statistics
#-------------------------------------------------------------------------------
secu1 <- seq(1,1116,12)
secu2 <- seq(12,1116,12)
# 
# m1 <- list()
# x  <- list()

# for(i in 1:length(sstcrops)){
#   for(j in 1:length(secu1)){
#     print(i)
#     print(length(sstcrops[[i]][[secu1[j]:secu2[j]]]@data@names))
#     
#     x[[j]] <- mean(sstcrops[[i]][[secu1[j]:secu2[j]]])
#   }
#    m1[[i]] <- x
# } # old nested for loop 

# foreach nested loop to calculate the yearly SST mean for each raster brick
sstCM <- foreach(i=1:length(sstcrops)) %:%
  foreach(j=1:length(secu1), .combine='c') %dopar%{ # add combine function "c" 
    x <- mean(sstcrops[[i]][[secu1[j]:secu2[j]]])
    y <- cellStats(x, "mean")
 return(y)
}

# rename the list of SST extents with the Stock area name
names(sstCM) <- names(sstcrops)

# Make time series data frame normal SST
Years=seq(1930,2022,1)
SSTchileStocks <- data.frame(Years=Years, sstCM)

# (column names already set: Years + one column per stock extent)

# END
write.table(SSTchileStocks, file.path(PROJ, "data", "derived", "SSTchileStocks_raw_2022.txt"), sep="\t", row.names = FALSE)
#-------------------------------------------------------------------------------
# smooth individual time series to reduce noise (smooth splines)
#-------------------------------------------------------------------------------
# loop to estimate the smooth splines with  spar=0.2
sstSmooth <- foreach(i=1:length(sstCM))%do%{
  x <- smooth.spline(sstCM[[i]], spar = 0.35)
  y <- x[["y"]]
  return(y)
}

# rename the list of SST extents with the Stock area name
names(sstSmooth) <- names(sstcrops)

SSTchileStocksSplines <- data.frame(Years=Years, sstSmooth)
plot(SSTchileStocksSplines$Years, SSTchileStocksSplines$CH,t='b')

write.table(SSTchileStocksSplines,
            file.path(PROJ, "data", "derived", "SSTchileStocks_2022.txt"),
            sep="\t", row.names = FALSE)
