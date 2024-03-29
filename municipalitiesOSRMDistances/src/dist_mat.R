################################################################################
# Calculate travel times and distances between Austrian municipalities
################################################################################
library(tidyverse)
library(sf)
#install.packages("osrm")
library(osrm)
################################################################################
# Load map of Austrian municipalities' gravity points (pop-weighted) from: 
# https://data.statistik.gv.at/web/meta.jsp?dataset=OGDEXT_GEM_MP_1
aut2022_muni <- read_sf("Shapefiles/STATISTIK_AUSTRIA_GEM_MP_20220101.shp", options = "ENCODING=WINDOWS-1252") %>% 
  rename(id = g_id, name = g_name)

muni_names <- aut2022_muni %>% 
  as_tibble() %>% 
  select(id, name)
################################################################################
# Load potential neighbors based on distances-matrix provided by Ruslan
dist_mat <- read_csv(file = "Municipalities' distances/results/distances.csv", col_names = T)

dist_mat_long <- dist_mat %>%
  mutate(id = names(dist_mat)) %>% 
  left_join(muni_names, by = 'id') %>% 
  gather(key = "Neighbors", value = 'distance', -c(id, name)) %>% 
  filter(!is.na(distance)) %>%
  left_join(muni_names, by = c('Neighbors' = 'id'), suffix = c('id', '_Neighbor'))  
################################################################################
my_src <- dist_mat_long %>%
  left_join(aut2022_muni, by = "id") %>% 
  st_sf()

my_dst <- dist_mat_long %>%
  left_join(aut2022_muni, by = c("Neighbors" = "id")) %>% 
  st_sf()
################################################################################
# Extract travel times and distances with OpenStreetMap
my_trip = NULL
for (i in 1:floor(nrow(my_dst)/50 + 1)){
  
  i1 = 50 * (i - 1) + 1
  i2 = min(i1 + 49, nrow(my_dst))
  message(paste0(i1, ":", i2))
  
  i_src = my_src[i1:i2,]
  i_dst = my_dst[i1:i2,]
  
  tmp1 <- osrmTable(src = i_src, 
                    dst = i_dst, 
                    measure = c('duration', 'distance'))
  
  my_dur <- tmp1$durations %>%
    as_tibble() %>% 
    mutate(trip = paste(i_src$nameid, i_dst$name_Neighbor, sep="--")) %>% 
    gather(src, mins, -trip) %>% 
    separate(trip, into = c("src", "dst"), sep = "--") %>% 
    distinct(src, dst, .keep_all = TRUE)
  
  my_trip <- tmp1$distances %>%
    as_tibble() %>% 
    mutate(trip = paste(i_src$nameid, i_dst$name_Neighbor, sep="--")) %>% 
    gather(src, meters, -trip) %>% 
    separate(trip, into = c("src", "dst"), sep = "--") %>% 
    distinct(src, dst, .keep_all = TRUE) %>% 
    left_join(my_dur) %>% 
    bind_rows(my_trip)
  
}
################################################################################











# test with Viennese districts
vie2022_muni <- aut2022_muni %>% 
  mutate(g_id = as.numeric(g_id)) %>% 
  filter(g_id > 90000)

villages <- c("Sölden", "Innsbruck", "Moosburg", "Krumpendorf am Wörthersee", "Pörtschach am Wörther See", "Weitensfeld im Gurktal")
tmap1 <- aut2022_muni %>% 
  filter(g_name %in% villages) %>% 
  column_to_rownames(var = "g_name") %>% 
  st_sf()
################################################################################
tmp1 <- osrmTable(src = tmap1, 
                  dst = vie2022_muni, 
                  measure = c('duration', 'distance'))

################################################################################
# install.packages("gmapsdistance")
library(gmapsdistance)
# # remotes::install_github("statistikat/STATcubeR")
# library(STATcubeR)
# install.packages('sf')
library(sf)
################################################################################
set.api.key("AIzaSyCZzB6XGv9J-S9PMHTGDDyjKyH6hwPaOI8")
################################################################################
# Some initial tests
gmapsdistance(origin = "New York", destination = "Washington DC")
gmapsdistance(origin = "Velden", destination = "Klagenfurt", mode = "driving")$Time / 60

gmapsdistance(origin = "Pörtschach am Wörthersee", destination = "Klagenfurt", mode = "driving")$Time / 60
gmapsdistance(origin = "Velden", destination = "Klagenfurt", mode = "driving")$Time / 60
gmapsdistance(origin = "Velden", destination = "Pörtschach am Wörthersee", mode = "driving")$Time / 60
gmapsdistance(origin = "Pörtschach am Wörthersee", destination = "Wien", mode = "driving")$Time / 60

gmapsdistance(origin = c("Pörtschach am Wörthersee", "Velden"),
              destination = c("Velden", "Klagenfurt", "Wien"), 
              mode = "driving", shape = "long")
################################################################################
# Load map of Austrian municipalities to get centroid GEO-coordinates
aut2022_muni <- sf::read_sf("STATISTIK_AUSTRIA_GEM_20220101.shp")

# Get names of Austrian municipalities
muni_names <- aut2022_muni$name
sum(duplicated(muni_names)) # some names appear multiple times
muni_names <- muni_names[grep("Wien", muni_names)]

# # Get names of Austrian municipalities from a random StatAT data set 
# # all_datasets <- od_list() # fetch a list of datasets from Statistik Austria
# population <- od_table("OGD_bevstandjbab2002_BevStand_2022")$tabulate()
# muni_names <- population %>% 
#   rename(Gemeinde = `Commune (aggregation by political district)`) %>% 
#   distinct(Gemeinde, `Time section`) %>% 
#   separate_wider_delim(cols = Gemeinde, names = c("Municipality", "GKZ"), delim = "<") %>% 
#   mutate(Municipality = trimws(Municipality)) %>% 
#   as_tibble() %>% 
#   pull(Municipality)
  
# derive distances between every combination of municipalities
muni_dist <- gmapsdistance(origin = muni_names, 
                           destination = muni_names, 
                           mode = "driving", shape = "long",
                           dep_date = "2023-11-16", dep_time = "18:30:00")
################################################################################
