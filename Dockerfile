FROM rocker/tidyverse:3.6.3

RUN install2.r --error \
    --deps TRUE \
    devtools \
    rlang \
    optparse \
    docstring \
    openxlsx \
    survival \
    survminer \
    msigdbr \
    pheatmap


RUN R -e "BiocManager::install('HPAanalyze')"
RUN R -e "BiocManager::install('pathview')"
RUN R -e "BiocManager::install('clusterProfiler')"
RUN R -e "BiocManager::install('edgeR')"
RUN R -e "BiocManager::install('EnhancedVolcano')"
RUN R -e "BiocManager::install('org.Hs.eg.db')"
RUN R -e "BiocManager::install('cowplot')"
RUN R -e "devtools::install_github('GuangchuangYu/scatterpie')"
RUN R -e "devtools::install_github('GuangchuangYu/ggplotify')"

ADD ./ /commonR