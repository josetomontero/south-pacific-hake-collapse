# Study-area polygons used by the SST extraction and by Figure 1.
#   Hake fishery management unit (Unidad de Pesqueria de merluza comun; Tascheri 2022, IFOP):
#   the sea between the northern limit of the Coquimbo Region (29 02' S) and 41 28.6' S, from
#   the coast to 60 nautical miles. Built as a 60 nm buffer of the coastline (UTM 19S), clipped
#   to the latitude band, land removed, largest piece kept.
#   Humboldt Current LME polygon (NOAA LME programme), re-exported as GeoJSON.
# Outputs: data/derived/fishery_unit_polygon.gpkg/.geojson, data/derived/LME_Humboldt.geojson
if (!exists("PROJ")) source(file.path(rprojroot::find_root(rprojroot::has_file("R/config.R")), "R", "config.R"))
suppressPackageStartupMessages(library(sf)); sf_use_s2(FALSE)

world <- st_transform(st_make_valid(st_as_sf(maps::map("world", plot = FALSE, fill = TRUE))), 4326)
hum   <- st_read(file.path(PROJ, "data", "external", "LME_Humboldt.gpkg"), quiet = TRUE)

LAT_N <- -(29 + 2 / 60); LAT_S <- -(41 + 28.6 / 60)
land  <- st_make_valid(st_union(st_geometry(st_crop(world, st_bbox(c(xmin = -78, xmax = -66, ymin = -46, ymax = -25), crs = 4326)))))
chile <- st_transform(land, 32719)
sea60 <- st_transform(st_difference(st_buffer(chile, 60 * 1852), chile), 4326)
band  <- st_as_sfc(st_bbox(c(xmin = -80, xmax = -70, ymin = LAT_S, ymax = LAT_N), crs = 4326))
unit  <- st_intersection(st_make_valid(sea60), band)
unit  <- st_sf(geometry = st_cast(st_union(unit), "POLYGON"))
unit  <- unit[which.max(st_area(unit)), ]
cat(sprintf("fishery management unit: %.0f km2, lat %.2f to %.2f, lon %.2f to %.2f\n", as.numeric(st_area(unit)) / 1e6,
            st_bbox(unit)["ymin"], st_bbox(unit)["ymax"], st_bbox(unit)["xmin"], st_bbox(unit)["xmax"]))
st_write(unit, file.path(PROJ, "data", "derived", "fishery_unit_polygon.gpkg"), quiet = TRUE, delete_dsn = TRUE)
st_write(unit, file.path(PROJ, "data", "derived", "fishery_unit_polygon.geojson"), driver = "GeoJSON", quiet = TRUE, delete_dsn = TRUE)
st_write(st_transform(hum, 4326), file.path(PROJ, "data", "derived", "LME_Humboldt.geojson"), driver = "GeoJSON", quiet = TRUE, delete_dsn = TRUE)
