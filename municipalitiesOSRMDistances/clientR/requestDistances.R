library(sf)
library(osrm)

year <- 2022

shp_file <- read_sf(sprintf("data/STATISTIK_AUSTRIA_GEM_MP_%s0101/", year),
                    options = "ENCODING=WINDOWS-1252")

shp_file <- st_transform(shp_file, 4326)
shp_file <- shp_file[which(shp_file$g_name != "Schwanberg"), ]
coordinates <- st_coordinates(shp_file[, "geometry"]) |>
  as.data.frame()

openstreetmap_response <- osrmTable(src = coordinates,
                                    dst = coordinates,
                                    measure = c("duration", "distance"),
                                    osrm.server = "server:5000/")

durations <- openstreetmap_response[["durations"]]
distances <- openstreetmap_response[["distances"]]

datetime <- format(Sys.time(), "%Y%m%d")
year_of_map <- 2015
write.csv(durations, sprintf("data/outputOSRM/durations_%s_%s.csv", year_of_map, datetime))
write.csv(distances, sprintf("data/outputOSRM/distances_%s_%s.csv", year_of_map, datetime))
