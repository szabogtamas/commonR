FROM rocker/tidyverse:3.6.3


RUN install2.r --error \
    --deps TRUE \
    optparse \
    docstring \
    survival \
    survminer \
    msigdbr \
    pheatmap \
    ggpubr \
    ggplotify \
    cowplot


RUN R -e "BiocManager::install('clusterProfiler')"
RUN R -e "BiocManager::install('EnhancedVolcano')"
