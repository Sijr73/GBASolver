# Base image https://hub.docker.com/u/rocker/
FROM --platform=linux/amd64 rocker/shiny:latest
ENV DEBIAN_FRONTEND=noninteractive

WORKDIR /app


RUN apt-get update -qq && apt-get -y dist-upgrade && apt-get clean && apt-get -y autoremove && rm -rf /var/lib/apt/lists/*

    

RUN R -e "install.packages(c('rjson','dplyr', 'rstudioapi', 'nloptr', 'readODS', 'Matrix', 'MASS', 'lpSolve', 'shiny', 'shinyMatrix', 'reactable', 'shinydisconnect', 'htmltools', 'shinyjs', 'shinybusy', 'apexcharter', 'reshape2', 'tidyr', 'ggplot2', 'shinyalert', 'rintrojs', 'later', 'ipoptr'))" && rm -rf /tmp/*
COPY GBApp/ .

CMD ["R", "-e", "shiny::runApp('/app/app.R', host='0.0.0.0', port=3838)"]