#!/bin/bash

# sbatch --mem=20GB --time=1-00:00:0 src/setup_vhl.sh

# wd
wd=/camp/project/tracerX/working/VHL_GERMLINE/tidda/vhl/ ; cd $wd

# create conda env (run on login node)
source ~/.bashrc
conda create -n vhl -c conda-forge -c bioconda -c r \
  r-essentials \
  r-base \
  r-rmarkdown \
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
  r-patchwork \
  r-purrr \
  bioconductor-scdblfinder \
  r-seurat \
  bioconductor-singler \
  r-janitor \
  r-bookdown \
  r-monocle3
