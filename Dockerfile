FROM rocker/tidyverse

# Not recommended with Rocker, but this is the only workaround that fixes DESeq2
#RUN sudo apt update
#RUN sudo apt-get install -y r-cran-lattice

RUN sudo apt-get update -y
RUN sudo apt-get install -y libxt-dev

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
    pheatmap \
    proj4


RUN R -e "BiocManager::install('HPAanalyze')"
RUN R -e "BiocManager::install('karyoploteR')"
RUN R -e "BiocManager::install('ggbio')"
RUN R -e "BiocManager::install('maftools')"
RUN R -e "BiocManager::install('progeny', version='3.12')"
RUN R -e "BiocManager::install('pathview')"
RUN R -e "BiocManager::install('clusterProfiler')"
RUN R -e "BiocManager::install('edgeR')"
RUN R -e "BiocManager::install('DESeq2')"
RUN R -e "BiocManager::install('DEFormats')"
RUN R -e "BiocManager::install('HTqPCR')"
RUN R -e "BiocManager::install('EnhancedVolcano', version='3.9')"
RUN R -e "BiocManager::install('org.Hs.eg.db')"
RUN R -e "BiocManager::install('cowplot')"
RUN R -e "devtools::install_github('GuangchuangYu/scatterpie')"
RUN R -e "devtools::install_github('GuangchuangYu/ggplotify')"

#RUN sudo R -e "devtools::install_github('kevinblighe/EnhancedVolcano')"

ADD ./ /commonR