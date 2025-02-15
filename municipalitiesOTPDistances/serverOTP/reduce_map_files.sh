#!/bin/bash

# You may try to change the tag to a date of format YYMMDD.
# For example, 17 June 2024 becomes 240617.
# Alas, as of 17.06.2024, daywise precision is available only for the last week,
# monthwise (the first day) for the current year, and for the earlier years,
# only the first day of January is available. This concerns only the ready-to-use snapshots.
tag=latest

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

# Osmium filter command
filter_command="osmium tags-filter -o temp.osm.pbf -f pbf,add_metadata=false"
filter_tags="w/highway wa/public_transport=platform wa/railway=platform w/park_ride=yes r/type=restriction r/type=route nwr/timezone"

# Process each file
for file in "${map_files[@]}"; do
    echo "Filtering $file..."
    $filter_command $file $filter_tags
    mv temp.osm.pbf "$file"
done

osmium merge "${map_files[*]}" -o /var/merged_map.osm.pbf
rm -f "${map_files[*]}"
