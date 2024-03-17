install.packages("devtools")
devtools::install_version("sf", version = "1.0-15")
devtools::install_version("osrm", version = "4.1.1")

library(sf)
library(osrm)

shpFile = read_sf("/Users/johnnymosby/repos/thesis_demography/data/STATISTIK_AUSTRIA_GEM_MP_20220101/", options = "ENCODING=WINDOWS-1252")
shpFile = st_transform(shpFile, 4326)
coordinates = st_coordinates(shpFile[, "geometry"]) |> 
    as.data.frame()

