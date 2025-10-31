#!/bin/sh

# start shiny server
exec Rscript /srv/shiny-server/server.R /srv/data/database.db 2>&1
#exec shiny-server 2>&1