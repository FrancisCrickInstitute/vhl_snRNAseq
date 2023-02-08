#!/bin/bash
# create and submit the split-pipe comb sbatch script for each experiment / genome
run_info=qc_${experiment}_${genome}_${sublibrary}_${parse_analysis_subdir/\//_}_int${do_integration}
script=src/sbatch/sbatch_vhl_${run_info}.sh
echo $run_info

cat << EOF > $script
#!/bin/bash
#SBATCH --job-name=parse_${run_info}
#SBATCH --time=1-00:00:0
#SBATCH --mem=50GB
#SBATCH -o /camp/project/tracerX/working/VHL_GERMLINE/tidda/vhl/log/parse_${run_info}_%j.out
#SBATCH -e /camp/project/tracerX/working/VHL_GERMLINE/tidda/vhl/log/parse_${run_info}_%j.err

echo "------------------------------"
echo "Experiment: $experiment"
echo "Genome: $genome"
echo "Sublibrary: $sublibrary"
echo "Parse analysis subdirectory: $parse_analysis_subdir"
echo "Script: $script"
echo "------------------------------"

# get directories
cd $(pwd) ; . src/config.sh

# run generate_QC_report.R
Rscript src/run_generate_qc_report.R \
  $experiment \
  $genome \
  $sublibrary \
  $parse_analysis_subdir \
  $do_integration
EOF

# submit the script
sbatch $script
