FROM rocker/tidyverse

# Not recommended with Rocker, but this is the only workaround that fixes DESeq2
#RUN sudo apt update
#RUN sudo apt-get install -y r-cran-lattice

RUN sudo apt-get update -y
#RUN sudo apt-get remove -y libxt-dev
RUN sudo apt-get install -y libxt-dev
#RUN sudo /sbin/ldconfig -v

RUN install2.r --error \
    --deps TRUE \
    devtools \
    rlang \
    proj4

RUN R -e "remove.packages('proj4')"
RUN R -e "BiocManager::install('EnhancedVolcano')"

RUN ls /usr/local/lib/R/site-library/proj4/libs/

ADD ./ /commonR