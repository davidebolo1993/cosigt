#!/bin/bash

#server
Rscript -e "shiny::runApp('app.r', host = '0.0.0.0', port = 3838)" $1
#local
#ssh -L 3838:localhost:3838 <user>@<ip>
