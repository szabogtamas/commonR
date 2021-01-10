FROM rocker/tidyverse:4.0.3

RUN sudo apt-get update -y
RUN sudo apt-get install -y libxt-dev

RUN R -e "BiocManager::install('DESeq2')"
RUN R -e "BiocManager::install('progeny')"
RUN R -e "BiocManager::install('cowplot')"

RUN install2.r --error \
    --deps TRUE \
    pheatmap

ADD ./ /commonR