#!/bin/bash
#SBATCH --job-name=parse_QC_230127_A01366_0343_AHGNCVDMXY_hg38_SHE5052A12_S165_all-well_DGE_filtered
#SBATCH --time=1-00:00:0
#SBATCH --mem=50GB
#SBATCH -o /camp/project/tracerX/working/VHL_GERMLINE/tidda/vhl/log/parse_QC_230127_A01366_0343_AHGNCVDMXY_hg38_SHE5052A12_S165_all-well_DGE_filtered_%j.out
#SBATCH -e /camp/project/tracerX/working/VHL_GERMLINE/tidda/vhl/log/parse_QC_230127_A01366_0343_AHGNCVDMXY_hg38_SHE5052A12_S165_all-well_DGE_filtered_%j.err

echo "------------------------------"
echo "Experiment: 230127_A01366_0343_AHGNCVDMXY"
echo "Genome: hg38"
echo "Sublibrary: SHE5052A12_S165"
echo "Parse analysis subdirectory: all-well/DGE_filtered"
echo "Script: src/sbatch/sbatch_vhl_QC_230127_A01366_0343_AHGNCVDMXY_hg38_SHE5052A12_S165_all-well_DGE_filtered.sh"
echo "------------------------------"

# get directories
cd /camp/project/tracerX/working/VHL_GERMLINE/tidda/vhl
ml R/4.2.0-foss-2021b

# run generate_QC_report.R
Rscript src/run_generate_QC_report.R   230127_A01366_0343_AHGNCVDMXY   hg38   SHE5052A12_S165   all-well/DGE_filtered
