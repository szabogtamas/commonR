FROM rocker/tidyverse:4.0.1

RUN sudo apt-get update -y
RUN sudo apt-get install -y libxt-dev

RUN install2.r --error \
    --deps TRUE \
    devtools \
    rlang \
    optparse \
    docstring \
    openxlsx

RUN R -e "BiocManager::install('DESeq2')"
RUN R -e "BiocManager::install('progeny')"

ADD ./ /commonR