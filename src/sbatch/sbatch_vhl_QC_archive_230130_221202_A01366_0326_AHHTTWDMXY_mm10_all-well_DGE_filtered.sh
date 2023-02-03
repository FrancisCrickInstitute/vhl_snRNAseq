#!/bin/bash
#SBATCH --job-name=parse_QC_archive_230130_221202_A01366_0326_AHHTTWDMXY_mm10_all-well_DGE_filtered
#SBATCH --time=1-00:00:0
#SBATCH --mem=50GB
#SBATCH -o /camp/project/tracerX/working/VHL_GERMLINE/tidda/vhl/log/parse_QC_archive_230130_221202_A01366_0326_AHHTTWDMXY_mm10_all-well_DGE_filtered_%j.out
#SBATCH -e /camp/project/tracerX/working/VHL_GERMLINE/tidda/vhl/log/parse_QC_archive_230130_221202_A01366_0326_AHHTTWDMXY_mm10_all-well_DGE_filtered_%j.err

echo "------------------------------"
echo "Experiment: archive_230130"
echo "Genome: 221202_A01366_0326_AHHTTWDMXY"
echo "Sublibrary: mm10"
echo "Parse analysis subdirectory: all-well/DGE_filtered"
echo "Script: src/sbatch/sbatch_vhl_QC_archive_230130_221202_A01366_0326_AHHTTWDMXY_mm10_all-well_DGE_filtered.sh"
echo "------------------------------"

# get directories
cd /camp/project/tracerX/working/VHL_GERMLINE/tidda/vhl
ml R/4.2.0-foss-2021b

# run generate_QC_report.R
Rscript src/run_generate_QC_report.R   archive_230130   221202_A01366_0326_AHHTTWDMXY   mm10   all-well/DGE_filtered
