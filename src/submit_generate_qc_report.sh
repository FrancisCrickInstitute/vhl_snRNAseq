#!/bin/bash
# generate QC report for snRNAseq experiment
run_info=generate_qc_report_${experiment/\,/_x_}_${genome}_${sublibrary/\,/_x_}_${parse_analysis_subdir/\//_}_int${do_integration}
script=src/sbatch/sbatch_${run_info}.sh
echo $run_info

cat << EOF > $script
#!/bin/bash
#SBATCH --job-name=parse_${run_info}
#SBATCH --time=1-00:00:0
#SBATCH --mem=150GB
#SBATCH -o /camp/project/tracerX/working/VHL_GERMLINE/tidda/vhl/log/${run_info}.out
#SBATCH -e /camp/project/tracerX/working/VHL_GERMLINE/tidda/vhl/log/${run_info}.err

echo "------------------------------"
echo "Experiment: $experiment"
echo "Genome: $genome"
echo "Sublibrary: $sublibrary"
echo "Parse analysis subdirectory: $parse_analysis_subdir"
echo "Integration: $do_integration"
echo "Sample subset: $sample_subset"
echo "Script: $script"
echo "------------------------------"

# codify params
experiment=$experiment
genome=$genome
sublibrary=$sublibrary
parse_analysis_subdir=$parse_analysis_subdir
do_integration=$do_integration
sample_subset=$sample_subset

# get directories
cd $(pwd) ; . src/config.sh

# run generate_QC_report.R
Rscript src/run_generate_qc_report.R \
  $experiment \
  $genome \
  $sublibrary \
  $parse_analysis_subdir \
  $do_integration \
  $sample_subset
EOF

# submit the script
sbatch $script