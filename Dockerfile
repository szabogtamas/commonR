FROM rocker/tidyverse:3.6.3


RUN install2.r --error \
    --deps TRUE \
    devtools \
    rlang \
    optparse \
    pheatmap \
    ggplotify \
    cowplot


RUN R -e "BiocManager::install('clusterProfiler')"
RUN R -e "BiocManager::install('EnhancedVolcano')"

ADD ./ /commonR