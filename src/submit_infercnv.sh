#!/bin/bash
# run inferCNV
run_info=run_infercnv_${nih_pid}_${analysis_mode}
script=src/sbatch/sbatch_${run_info}.sh
echo $run_info

cat << EOF > $script
#!/bin/bash
#SBATCH --job-name=parse_${run_info}
#SBATCH --time=3-00:00:0
#SBATCH --mem=120GB
#SBATCH -o /camp/project/tracerX/working/VHL_GERMLINE/tidda/vhl/log/${run_info}.out
#SBATCH -e /camp/project/tracerX/working/VHL_GERMLINE/tidda/vhl/log/${run_info}.err

# codify params
nih_pid=$nih_pid
analysis_mode=$analysis_mode

# get directories
cd $(pwd) ; . src/config.sh

# run generate_QC_report.R
Rscript src/run_infercnv.R $nih_pid $analysis_mode
EOF

# submit the script
sbatch $script
