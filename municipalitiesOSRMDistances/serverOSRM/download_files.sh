#!/bin/bash

# You may try to change the tag to a date of format YYMMDD.
# For example, 17 June 2024 becomes 240617.
# Alas, as of 17.06.2024, daywise precision is available only for the last week,
# monthwise (the first day) for the current year, and for the earlier years,
# only the first day of January is available. This concerns only the ready-to-use snapshots.
tag=140101

# Austria and the bordering countries. The bordering countries are limited only to the
# respective ordering regions, if possible, to save time on the data download and calculation.
# THE TOTAL SIZE is around 3.5GB while only Austria is around 700MB
links=(
	"https://download.geofabrik.de/europe/austria-${tag}.osm.pbf"
	# "https://download.geofabrik.de/europe/czech-republic-${tag}.osm.pbf"
	# "https://download.geofabrik.de/europe/germany/bayern-${tag}.osm.pbf"
	# "https://download.geofabrik.de/europe/hungary-${tag}.osm.pbf"
	# "https://download.geofabrik.de/europe/italy/nord-est-${tag}.osm.pbf"
	# "https://download.geofabrik.de/europe/liechtenstein-${tag}.osm.pbf"
	# "https://download.geofabrik.de/europe/slovakia-${tag}.osm.pbf"
	# "https://download.geofabrik.de/europe/slovenia-${tag}.osm.pbf"
	# "https://download.geofabrik.de/europe/switzerland-${tag}.osm.pbf"
)

regions=(
	"austria"
	# "czech-republic"
	# "bayern"				# Bavaria in Germany
	# "hungary"
	# "nord-est"				# Northeastern Italy
	# "liechtenstein"
	# "slovakia"
	# "slovenia"
	# "switzerland"
)

map_files=( "${regions[@]/%/-${tag}.osm.pbf}" )

echo "Downloading files"
for link in "${links[@]}"
do
    wget "$link"
done