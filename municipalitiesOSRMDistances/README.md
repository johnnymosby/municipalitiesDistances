# Distances between municipalities

In this part of the thesis, the distances between the Austrian municipalities are calculated by taking the population density adjusted centres as points representing the municipalities (publicly accessible from Statistics Austria). The distances are provided by OpenStreetMap Routing Machine, which is installed locally for (potentially (no comparison with requesting an official routing machine was conducted)) faster requesting of the data.

## Tools used
1. **Map data** - OpenStreetMap provided as snapshots by [Geofabrik.de](https://download.geofabrik.de/).\
OpenStreetMap is a geographic database that changes every time someone edits it. The only way to get the data directly from OSM as a single file is to request it for [the whole world](https://wiki.openstreetmap.org/wiki/Planet.osm), which is prohibitively huge. Thankfully, there are extract providers such as Geofabrik that allow access to OSM data in specific countries. In the case of Geofabrik, the data is available from days to years, depending on how far back in time you go.
2. **Program to operate the map data** - [Open Source Routing Machine](https://project-osrm.org/).\
OSRM is a program which constructs a route network (graph) based on the OSM data provided. It optimises the network and returns a route between any two points extremely fast (within a few milliseconds).
3. **Program to request the distances** - [`osrm` library in R](https://cran.r-project.org/web/packages/osrm/index.html).\
For more straightforward interaction with OSRM API, I recommend using `osrm` library in R. Namely, `osrmTable` function returns travel time and distance tables between all provided points and does it within minutes.
4. **Optional: replicating the computing environment** - [Docker](https://www.docker.com/).\
To facilitate the replication of the results by other people, I used Docker containerization. In essence, it allows the computation to be isolated from other computer resources. By doing so, you can be sure that your setup will work on any computer with Docker installed, and you will not have dependency issues or other problems.
