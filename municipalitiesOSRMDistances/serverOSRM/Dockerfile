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

RUN echo "osrm-routed merged_map.osrm --algorithm mld --max-table-size 8000 --port 1234" > ./runserver.sh && \
	chmod +x ./runserver.sh

ENTRYPOINT ./runserver.sh

# CMD ["tail", "-f", "/dev/null"]

# CMD [ "osrm-routed", "merged_map.osrm", "--algorithm mld", "--max-table-size 8000" ] 
# , "--algorithm mld", "--max-table-size 8000"
