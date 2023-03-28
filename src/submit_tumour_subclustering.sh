#!/bin/bash
# generate QC report for all tumour / sample subclusters
run_info=tumour_subclustering
script=src/sbatch/sbatch_${run_info}.sh
echo $run_info

cat << EOF > $script
#!/bin/bash
#SBATCH --job-name=parse_${run_info}
#SBATCH --time=1-00:00:0
#SBATCH --mem=120GB
#SBATCH -o /camp/project/tracerX/working/VHL_GERMLINE/tidda/vhl/log/${run_info}.out
#SBATCH -e /camp/project/tracerX/working/VHL_GERMLINE/tidda/vhl/log/${run_info}.err

# get directories
cd $(pwd) ; . src/config.sh

# run generate_QC_report.R
Rscript src/run_tumour_subclustering.R
EOF

# submit the script
sbatch $script
