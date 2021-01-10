FROM rocker/tidyverse:4.0.3

RUN sudo apt-get update -y
RUN sudo apt-get install -y libxt-dev

RUN install2.r --error \
    --deps TRUE \
    pheatmap

RUN R -e "BiocManager::install('DESeq2')"
RUN R -e "BiocManager::install('progeny')"
RUN R -e "BiocManager::install('EnhancedVolcano')"

ADD ./ /commonR