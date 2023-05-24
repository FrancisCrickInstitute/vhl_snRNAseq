#!/bin/bash
# generate QC report for snRNAseq experiment
run_info=generate_qc_report_PDOs
script=src/sbatch/sbatch_${run_info}.sh
mem=62
echo $run_info

cat << EOF > $script
#!/bin/bash
#SBATCH --job-name=parse_${run_info}
#SBATCH --time=1-00:00:0
#SBATCH --mem=${mem}GB
#SBATCH -o /camp/project/tracerX/working/VHL_GERMLINE/tidda/vhl/log/${run_info}.out
#SBATCH -e /camp/project/tracerX/working/VHL_GERMLINE/tidda/vhl/log/${run_info}.err

echo "------------------------------"
echo "Experiment: Daqi's PDOs"
echo "Genome: hg38"
echo "Output directory: out/PDOs/"
echo "Script: $script"
echo "Memory allocation: $mem GB"
echo "------------------------------"

# get directories
cd $(pwd) ; . src/config.sh

# run generate_QC_report.R
Rscript src/PDOs/run_generate_qc_report.R
EOF

# submit the script
sbatch $script
