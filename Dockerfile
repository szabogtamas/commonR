FROM rocker/tidyverse

# Not recommended with Rocker, but this is the only workaround that fixes DESeq2
#RUN sudo apt update
#RUN sudo apt-get install -y r-cran-lattice

RUN sudo apt-get update -y
RUN sudo apt-get install -y libxt-dev

RUN install2.r --error \
    --deps TRUE \
    devtools \
    rlang \
    proj4

RUN R -e "BiocManager::install('EnhancedVolcano')"

ADD ./ /commonR