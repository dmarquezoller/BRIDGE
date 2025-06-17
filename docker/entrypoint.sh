#!/bin/sh

# start shiny server
exec Rscript /srv/shiny-server/server.R /srv/data/database.sql 2>&1
