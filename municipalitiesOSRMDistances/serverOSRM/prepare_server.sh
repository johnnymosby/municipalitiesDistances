#!/bin/bash

# You may try to change the tag to a date of format YYMMDD.
# For example, 17 June 2024 becomes 240617.
# Alas, as of 17.06.2024, daywise precision is available only for the last week,
# monthwise (the first day) for the current year, and for the earlier years,
# only the first day of January is available. This concerns only the ready-to-use snapshots.
tag=180101

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


osmium merge "${map_files[*]}" -o merged_map.osm
osrm-extract merged_map.osm -p car.lua
osrm-partition merged_map.osrm
osrm-customize merged_map.osrm

# deleting files used to produce merged_map.osrm
rm merged_map.osm
rm -f "${map_files[*]}"
