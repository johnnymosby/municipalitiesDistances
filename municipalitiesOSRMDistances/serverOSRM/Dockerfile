# syntax=docker/dockerfile:1

FROM osrm/osrm-backend:v5.25.0

WORKDIR /opt

COPY ./download_files.sh ./download_files.sh
COPY ./prepare_server.sh ./prepare_server.sh

# The OSRM image is based on archived version of Debian
RUN echo "deb http://archive.debian.org/debian stretch main" > /etc/apt/sources.list && \
    apt-get update && \
    apt-get install -y \
			wget \
			osmium-tool
RUN chmod +x ./download_files.sh && \
	bash ./download_files.sh
RUN chmod +x ./prepare_server.sh && \
	bash ./prepare_server.sh

ENTRYPOINT [ "osrm-routed", "merged_map.osrm" ]
CMD [ "--algorithm mld", "--max-table-size 8000" ]