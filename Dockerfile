FROM rocker/shiny-verse:4.5.1
ARG DEBIAN_FRONTEND=noninteractive
ENV SHINY_LOG_STDERR=1
ENV TZ=UTC

RUN apt-get update \     
    && apt-get clean \
    && apt-get autoremove -y \
    && apt-get autoclean -y \
    && apt-get install -y tzdata wget libxml2-dev libglpk-dev libnetcdf-dev perl-base gcc-14 g++-14 \
    libstdc++-14-dev libc++-dev libc++abi-dev gfortran-14 libgfortran5 libcurl4-openssl-dev libssl-dev

# Set versions for gcc, g++, and gfortran
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-14 100 \
    && update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-14 100 \
    && update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/gfortran-14 100 \
    && update-alternatives --install /usr/bin/cc cc /usr/bin/gcc 30 \
    && update-alternatives --set cc /usr/bin/gcc \
    && update-alternatives --set gcc /usr/bin/gcc-14 \
    && update-alternatives --set g++ /usr/bin/g++-14 \        
    && update-alternatives --set gfortran /usr/bin/gfortran-14

# Application specific packages
RUN mkdir -p /packages/R \
    && mkdir -p /srv/data
COPY renv.lock /packages/R

ENV R_REMOTES_UPGRADE="never"
RUN R -e "install.packages('renv', repos = 'https://cloud.r-project.org/')"
RUN R -e "renv::restore('/packages/R')"

COPY R/BRIDGE /packages/R/BRIDGE
RUN R -e "devtools::install('/packages/R/BRIDGE')" 

RUN chown shiny:shiny /var/lib/shiny-server \ 
    && rm -rf /srv/shiny-server/* \
    && chown shiny:shiny /srv/shiny-server

ADD --chown=shiny:shiny ./docker/entrypoint.sh /entrypoint.sh
ADD --chown=shiny:shiny ./app.R /srv/shiny-server/server.R
#ADD --chown=shiny:shiny ./R/modules /srv/shiny-server/R/modules
#ADD --chown=shiny:shiny ./R/imports.R /srv/shiny-server/R/imports.R
#ADD --chown=shiny:shiny ./R/global.R /srv/shiny-server/R/global.R

RUN chmod +x /entrypoint.sh

EXPOSE 3838
USER shiny

# CMD R -e 'shiny::runApp("/srv/shiny-server/server.R", port = 3838, host = "0.0.0.0")'
ENTRYPOINT ["/entrypoint.sh"]
