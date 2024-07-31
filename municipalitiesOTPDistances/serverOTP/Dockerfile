# FROM debian:bullseye AS filtering

# RUN apt-get update && \
#     apt-get install -y \
#       wget \
#       osmium-tool

# COPY ./download_files.sh ./download_files.sh
# RUN chmod +x ./download_files.sh \
#     && ./download_files.sh

# RUN echo 'Reducing map files'
# COPY ./reduce_map_files.sh ./reduce_map_files.sh
# RUN chmod +x ./reduce_map_files.sh  \
#     && ./reduce_map_files.sh 
# ENTRYPOINT [ "tail", "-f", "/dev/null" ]


# Use OpenJDK 11 as the base image
FROM openjdk:21-slim AS builder

RUN apt-get update && apt-get install -y wget
RUN wget -O otp.jar https://repo1.maven.org/maven2/org/opentripplanner/otp/2.5.0/otp-2.5.0-shaded.jar

# COPY --from=filtering /var/merged_map.osm.pbf /var/opentripplanner/merged_map.osm.pbf
# COPY ./GTFS/merged_map.osm.pbf  /var/opentripplanner/merged_map.osm.pbf
COPY ./GTFS /var/opentripplanner

# WORKDIR /var/opentripplanner

# RUN java -Xmx16G -jar /otp.jar --build --save  .

RUN java -Xmx16G -jar /otp.jar --build --save /var/opentripplanner
ENTRYPOINT [ "tail", "-f", "/dev/null" ]
# ENTRYPOINT [ "java", "-Xmx16G", "-jar", "/otp.jar", "--build", "--save", "/var/opentripplanner" ]

# RUN $ java -Xmx2G -jar *.jar --build --save --serve /var/opentripplanner/ /graph.obj



FROM opentripplanner/opentripplanner:latest

COPY --from=builder /var/opentripplanner/graph.obj /var/opentripplanner/graph.obj
RUN wget -O otp.jar https://repo1.maven.org/maven2/org/opentripplanner/otp/2.5.0/otp-2.5.0-shaded.jar
# RUN java -Xmx2G -jar otp.jar --serve /graph.obj

ENTRYPOINT [ "java", "-Xmx4G", "-jar", "otp.jar", "--load", "/var/opentripplanner/" ]
