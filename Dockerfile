# Base image https://hub.docker.com/u/rocker/
FROM rocker/verse

MAINTAINER "Carl Boettiger and Dirk Eddelbuettel" rocker-maintainers@eddelbuettel.com

RUN apt-get update  \
    && apt-get install -yq --no-install-recommends groff \
    && rm -rf /var/lib/apt/lists/* 
	
## create directories
RUN mkdir -p /01_data
RUN mkdir -p /02_code

## copy files
COPY /02_code/install_packages.R /02_code/install_packages.R
COPY /02_code/snowCover.R /02_code/snowCover.R

## install R-packages
RUN Rscript /02_code/install_packages.R
<<<<<<< HEAD

=======
>>>>>>> d554cf970d6379eef7b8fcbb7c1b21dda4073f2d
