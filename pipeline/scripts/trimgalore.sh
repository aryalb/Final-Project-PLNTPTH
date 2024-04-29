#!/bin/bash
#SBATCH --account=PAS2700
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-trimgalore-%j.out
set -euo pipefail

# Load the software
module load miniconda3/24.1.2-py310
conda activate /fs/ess/PAS0471/jelmer/conda/trimgalore

# Copy the placeholder variables
if [[ ! "$#" -eq 2 ]]; then
    echo "Error: You provided $# arguments, while 2 are required."
    echo "Usage: trimgalore.sh <R1-FASTQ> <outdir>"
    echo "Example: sbatch trimgalore.sh data/fastq/A_R1.fastq.gz results/trimgalore"
    echo "Your arguments: $*"
    exit 1
fi
R1_in=$1
outdir=$2

# Infer the R2 FASTQ file name
R2_in=${R1_in/_R1/_R2}

# Report start of script and variables
echo "# Starting script trimgalore.sh"
date
echo "# Input R1 FASTQ file:      $R1_in"
echo "# Input R2 FASTQ file:      $R2_in"
echo "# Output dir:               $outdir"
echo

# Create the output dir (with a subdir for Slurm logs)
mkdir -p "$outdir"/logs

# Run TrimGalore
trim_galore \
    --paired \
    --fastqc \
    --output_dir "$outdir" \
    --cores 8 \
    "$R1_in" \
    "$R2_in"

# Rename the output files - TrimGalore's output file names are unwieldy
# with two separate read direction indicators (_R1 and _1).
echo
echo "# Renaming the output files:"
sample_id=$(basename "$R1_in" _R1.fastq.gz)
R1_out_init="$outdir"/"$sample_id"_R1_val_1.fq.gz # This will be the initial R1 output file 
R2_out_init="$outdir"/"$sample_id"_R2_val_2.fq.gz # This will be the initial R2 output file 
R1_out_final="$outdir"/"$sample_id"_R1.fastq.gz
R2_out_final="$outdir"/"$sample_id"_R2.fastq.gz
mv -v "$R1_out_init" "$R1_out_final"
mv -v "$R2_out_init" "$R2_out_final"

# Report software version, end of script, and date
echo
echo "# Used TrimGalore version:"
trim_galore --version | grep version
echo "# Done with script trimgalore.sh"
date
