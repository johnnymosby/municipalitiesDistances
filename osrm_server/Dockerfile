# syntax=docker/dockerfile:1

FROM osrm/osrm-backend:v5.25.0

WORKDIR /opt

# The OSRM image is based on archived version of Debian
RUN echo "deb http://archive.debian.org/debian stretch main" > /etc/apt/sources.list && \
	apt-get update && \
	apt-get install -y wget && \
	wget https://download.geofabrik.de/europe/austria-240315.osm.pbf
RUN osrm-extract austria-240315.osm.pbf -p car.lua && \
	osrm-partition austria-240315.osrm && \
	osrm-customize austria-240315.osrm
RUN echo "osrm-routed austria-240315.osrm --algorithm mld --max-table-size 8000" > ./runserver.sh && \
	chmod +x ./runserver.sh


ENTRYPOINT ./runserver.sh