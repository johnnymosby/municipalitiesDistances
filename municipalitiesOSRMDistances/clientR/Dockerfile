FROM r-base:4.4

WORKDIR /script

COPY download_libraries.R .
COPY requestDistances.R .

RUN Rscript download_libraries.R \
	&& chmod 777 requestDistances.R

WORKDIR /var/distances

ENTRYPOINT [ "Rscript", "/script/requestDistances.R" ]
