FROM rocker/tidyverse:4.0.1

RUN install2.r --error \
    --deps TRUE \
    pheatmap

RUN R -e "BiocManager::install('DESeq2')"
RUN R -e "BiocManager::install('progeny')"

ADD ./ /commonR