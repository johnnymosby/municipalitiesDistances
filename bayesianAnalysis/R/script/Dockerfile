FROM r-base:4.3

WORKDIR /${RPROJECTNAME}

COPY . .

RUN R -e "install.packages('renv', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv')"

