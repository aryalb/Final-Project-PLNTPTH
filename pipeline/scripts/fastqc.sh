#!/bin/bash
#SBATCH --account=PAS2700
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-fastqc-%j.out
set -euo pipefail

# Load the software
module load miniconda3/24.1.2-py310
conda activate /fs/ess/PAS0471/jelmer/conda/fastqc

# Process the command-line arguments
if [[ ! "$#" -eq 2 ]]; then
    echo "Error: You provided $# arguments, while 4 are required."
    echo "Usage: fastqc.sh <FASTQ-file> <output-dir>"
    echo "Example: sbatch fastqc.sh data/fastq/A_R1.fastq.gz results/fastqc"
    echo "Your arguments: $*"
    exit 1
fi
fastq_file=$1
outdir=$2

# Report start of script and variables
echo "# Starting script fastqc.ch"
date
echo "# Input FASTQ file:   $fastq_file"
echo "# Output dir:         $outdir"
echo

# Create the output dir (with a subdir for Slurm logs)
mkdir -p "$outdir"/logs

# Run FastQC
fastqc --outdir "$outdir" "$fastq_file"

# Report software version, end of script, and date
echo
echo "# Used FastQC version:"
fastqc --version
echo "# Done with script fastqc.sh"
date
echo
