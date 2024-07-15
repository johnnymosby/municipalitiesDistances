FROM osrm/osrm-backend:v5.25.0

WORKDIR /opt

COPY ./download_files.sh ./download_files.sh
COPY ./prepare_server.sh ./prepare_server.sh

# The OSRM image is based on archived version of Debian
RUN echo "deb http://archive.debian.org/debian stretch main" > /etc/apt/sources.list && \
    apt-get update && \
    apt-get install -y \
			wget \
			osmium-tool \
    && chmod +x ./download_files.sh \
    && ./download_files.sh \
    && chmod +x ./prepare_server.sh \
    && ./prepare_server.sh \
    &&  echo "osrm-routed merged_map.osrm --algorithm mld --max-table-size 8000" > ./runserver.sh \
    && chmod +x ./runserver.sh

ENTRYPOINT ./runserver.sh
