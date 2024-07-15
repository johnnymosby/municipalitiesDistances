FROM r-base:4.4.0

WORKDIR /script

RUN apt-get update \
    && apt-get install -y \
      r-cran-devtools \
      r-cran-sf \
      iputils-ping

COPY download_libraries.R .
RUN Rscript download_libraries.R
COPY requestDistances.R .
RUN chmod 777 requestDistances.R


WORKDIR /var/distances

ENTRYPOINT [ "Rscript", "/script/requestDistances.R" ]

# CMD ["tail", "-f", "/dev/null"]