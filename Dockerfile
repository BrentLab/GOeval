# get the base image, the rocker/verse has R, RStudio and pandoc
FROM r-base
# Get and install system dependencies

ENV RENV_VERSION 0.16.0
RUN R -e "install.packages(c('remotes'), repos = c(CRAN = 'https://cran.wustl.edu'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

# it would be a good idea to figure out what
# each of these do. My method of dealing with this
# is to reduce it, try installing, and then add
# packages to deal with errors. I don't adjust this
# as frequently as I should.
RUN  apt-get update && \
     apt-get install -y --no-install-recommends \
      software-properties-common \
      dirmngr \
      wget \
        build-essential \
        libssl-dev \
        libxml2-dev \
        libglpk-dev \
      libcurl4-openssl-dev \
      libfontconfig1-dev \
      libharfbuzz-dev \
      libfribidi-dev \
      libtiff-dev


# Clean up
RUN apt-get autoremove -y

WORKDIR /project
COPY renv.lock renv.lock
COPY renv/ renv/

COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R
# COPY renv/settings.dcf renv/settings.dcf

# note: update this path as necessary based on the r-base r version
# and what you make your WORKDIR
ENV R_LIBS /project/renv/library/R-4.2/x86_64-pc-linux-gnu

# Using renv takes a really, really long time. I don't know
# renv terribly well, so it is possible that there are settings
# that could reduce the time it takes to install.
RUN R -e "renv::restore()"

WORKDIR /project/src
COPY . .
WORKDIR /project
RUN R -e "renv::activate();renv::install('./src', dependencies = TRUE)"
RUN rm -rf src

##### METHOD 2 install from github releases using remotes::install_github #####
#
# You can skip the renv portion and just use remotes::install_github(..., dependencies=TRUE)
