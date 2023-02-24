#!/bin/bash

# wd
wd=/camp/project/tracerX/working/VHL_GERMLINE/tidda/vhl/ ; cd $wd

# create conda env (run on login node)
source ~/.bashrc
conda create -n vhl -c conda-forge -c bioconda -c r \
  r-environment \
  r-essentials \
  r-base \
  r-devtools \
  bioconductor-celldex \
  r-clustree \
  bioconductor-dittoseq \
  r-ggforce \
  r-gridextra \
  r-intrinsicdimension \
  r-knitr \
  r-lemon \
  r-magrittr \
  r-matrix \
  r-patchwork \
  r-pheatmap \
  r-purrr \
  bioconductor-scdblfinder \
  r-seurat \
  bioconductor-singler


