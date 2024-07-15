library(sf)
library(osrm)

year <- 2020

shp_file <- read_sf(sprintf("STATISTIK_AUSTRIA_GEM_MP_%s0101/", year),
                    options = "ENCODING=WINDOWS-1252")

shp_file <- st_transform(shp_file, 4326)
shp_file <- shpFile[which(shp_file$g_name != "Schwanberg"), ]
coordinates <- st_coordinates(shp_file[, "geometry"]) |>
  as.data.frame()

openstreetmap_response <- osrmTable(src = coordinates,
                                    dst = coordinates,
                                    measure = c("duration", "distance"),
                                    osrm.server = "http://localhost:5000/")

durations <- openstreetmapResponse[["durations"]]
distances <- openstreetmapResponse[["distances"]]

datetime <- format(Sys.time(), "%Y%m%d_%H%M")
write.csv(durations, sprintf("durations_%s_%s.csv", year, datetime))
write.csv(distances, sprintf("distances_%s_%s.csv", year, datetime))