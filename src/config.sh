#!/bin/bash

# directories
base_dir=/camp/project/tracerX/working/VHL_GERMLINE/tidda/
wd=$base_dir/vhl/
cd $wd

. ~/.bashrc

# load conda env
ml purge
conda activate vhl2
