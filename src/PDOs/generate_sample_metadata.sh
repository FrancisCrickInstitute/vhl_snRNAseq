#!/bin/bash
cd /camp/project/tracerX/working/VHL_GERMLINE/tidda/vhl/
data_dir=/nemo/lab/turajlics/home/users/dengd/scRNA_seq/proj_PDOT_Parse/output/Parse_lib_1/

echo -e 'sample\tnih_pid' > out/PDOs/sample_metadata.tsv
for file in $data_dir/K*html ; do
  basename ${file/_analysis_summary.html/} |
  awk -F'_' '{print $0"\t"$1}' \
  >> out/PDOs/sample_metadata.tsv
done
