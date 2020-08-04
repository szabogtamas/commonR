FROM rocker/tidyverse:3.6.3


RUN install2.r --error \
    --deps TRUE \
    survival \
    survminer \
    msigdbr \
    ggpubr \
    optigrab


RUN R -e "BiocManager::install('clusterProfiler')"
