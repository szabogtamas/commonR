FROM rocker/tidyverse:4.0.3

RUN R -e "BiocManager::install('DESeq2')"
RUN R -e "BiocManager::install('progeny')"
RUN R -e "BiocManager::install('EnhancedVolcano')"

ADD ./ /commonR