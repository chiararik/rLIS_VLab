# Base image https://hub.docker.com/u/rocker/
FROM rocker/verse

MAINTAINER "Carl Boettiger and Dirk Eddelbuettel" rocker-maintainers@eddelbuettel.com

RUN apt-get update  \
    && apt-get install -yq --no-install-recommends groff \
    && rm -rf /var/lib/apt/lists/* 


RUN apt-get update
RUN apt-get -y install libgdal-dev
RUN apt-get -y install libudunits2-dev

RUN mkdir -p /02_code

## copy files
COPY install_packages.R /02_code/install_packages.R


## install R-packages
RUN Rscript /02_code/install_packages.R

