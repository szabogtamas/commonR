FROM rocker/tidyverse:3.6.3


RUN install2.r --error \
    --deps TRUE \
    optparse \
    docstring \
    survival \
    survminer \
    msigdbr \
    pheatmap \
    EnhancedVolcano \
    ggpubr \
    ggplotify \
    cowplot


RUN R -e "BiocManager::install('clusterProfiler')"
