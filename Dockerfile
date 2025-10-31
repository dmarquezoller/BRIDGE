FROM rocker/shiny-verse:4.5.1
ARG DEBIAN_FRONTEND=noninteractive
ENV SHINY_LOG_STDERR=1
ENV TZ=UTC

RUN apt update -y\     
    && apt clean -y\
    && apt autoremove -y \
    && apt autoclean -y \
    && apt update -y \
    && apt upgrade -y
RUN rm -rf /var/lib/apt/lists/* \
    && apt update -y && apt install -y build-essential software-properties-common
RUN apt install -y tzdata wget libxml2-dev libglpk-dev libnetcdf-dev perl-base gcc-14 g++-14 libstdc++-14-dev libc++-dev libc++abi-dev gfortran-14 libgfortran5 libcurl4-openssl-dev libssl-dev liblzma-dev

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
RUN R -e "install.packages('renv', repos = 'https://rstudio.r-universe.dev')"
RUN R -e "renv::restore('/packages/R')"

COPY R/BRIDGE /packages/R/BRIDGE
RUN R -e "devtools::install('/packages/R/BRIDGE')" 

RUN chown shiny:shiny /var/lib/shiny-server \ 
    && rm -rf /srv/shiny-server/* \
    && chown shiny:shiny /srv/shiny-server

ENV BRIDGE_PORT=3838
ENV BRIDGE_DB_PATH=/srv/data/database.db
RUN touch /srv/data/database.db && chown -R shiny:shiny /srv/data && chmod -R 775 /srv/data/database.db
# ADD --chown=shiny:shiny ./docker/entrypoint.sh /entrypoint.sh
# RUN chmod +x /entrypoint.sh
ADD --chown=shiny:shiny ./app.R /srv/shiny-server/app.R

EXPOSE 3838
USER shiny

CMD ["R", "-e", "options(shiny.host='0.0.0.0', shiny.port=3838); shiny::runApp('/srv/shiny-server/app.R')"]
# CMD ["R", "-e", "options(shiny.host='0.0.0.0', shiny.port=3838, shiny.fullstacktrace=TRUE, shiny.trace=TRUE); shiny::runApp('/srv/shiny-server/app.R')"]  # Debug with full stack trace
# ENTRYPOINT ["/entrypoint.sh"]
# Build with docker build --network=host -t bridge:latest .