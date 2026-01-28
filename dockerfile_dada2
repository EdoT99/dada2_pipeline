FROM bioconductor/bioconductor_docker:RELEASE_3_18

# FROM rocker/r-ver:4.3.2
# # 1.  run this to ensure the OS has what the R packages need
# RUN apt-get update && apt-get install -y \
#     libxml2-dev \
#     libssl-dev \
#     libcurl4-openssl-dev \
#     zlib1g-dev \
#     && rm -rf /var/lib/apt/lists/*

RUN R -e "BiocManager::install(c('dada2', 'jsonlite', 'ShortRead'), ask=FALSE, Ncpus=4)"

WORKDIR /scripts
