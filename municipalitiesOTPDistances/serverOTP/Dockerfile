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
    && ./prepare_server.sh


# Use OpenJDK 11 as the base image
FROM opentripplanner:latest

COPY --from=0 /opt/merged_map.osrm merged_map.osrm 

RUN $ java -Xmx2G -jar *.jar --build --save --serve /var/opentripplanner/ /graph.obj

FROM opentripplanner:latest

COPY --from=1 /graph.obj /graph.obj

CMD [ "java", "-Xmx4G", "-jar", "*.jar", "--load", "/graph.obj" ]
