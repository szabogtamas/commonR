FROM rocker/tidyverse:4.0.1

RUN sudo apt-get update -y
RUN sudo apt-get install -y libxt-dev

RUN install2.r --error \
    --deps TRUE \
    devtools \
    rlang \
    optparse \
    docstring \
    openxlsx \
    #survival \
    #survminer \
    msigdbr \
    pheatmap

RUN R -e "BiocManager::install('DESeq2')"
RUN R -e "BiocManager::install('edgeR')"
RUN R -e "BiocManager::install('DEFormats')"
RUN R -e "BiocManager::install('HTqPCR')"
RUN R -e "BiocManager::install('org.Hs.eg.db')"
RUN R -e "BiocManager::install('clusterProfiler')"
RUN R -e "BiocManager::install('pathview')"
RUN R -e "BiocManager::install('progeny')"
RUN R -e "BiocManager::install('HPAanalyze')"
RUN R -e "BiocManager::install('karyoploteR')"
RUN R -e "BiocManager::install('ggbio')"
RUN R -e "BiocManager::install('maftools')"
RUN R -e "BiocManager::install('cowplot')"
RUN R -e "devtools::install_github('GuangchuangYu/scatterpie')"
RUN R -e "devtools::install_github('GuangchuangYu/ggplotify')"

ADD ./ /commonR