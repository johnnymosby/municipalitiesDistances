FROM r-base:4.4.0

WORKDIR /script

COPY download_libraries.R .
COPY requestDistances.R .

RUN apt-get update \
    && apt-get install -y \
      r-cran-devtools \
      r-cran-sf \
      iputils-ping \
    && Rscript download_libraries.R \
    && chmod 777 requestDistances.R

WORKDIR /var/distances

ENTRYPOINT [ "Rscript", "/script/requestDistances.R" ]
