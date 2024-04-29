#!/bin/bash
#SBATCH --account=PAS2700
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-multiqc-%j.out
set -euo pipefail

# Load the software
module load miniconda3/24.1.2-py310
conda activate /fs/ess/PAS0471/jelmer/conda/multiqc

# Process the command-line arguments
if [[ ! "$#" -eq 2 ]]; then
    echo "Error: You provided $# arguments, while 2 are required."
    echo "Usage: multiqc.sh <input-dir> <output-dir>"
    echo "Example: sbatch multiqc.sh results/ results/multiqc"
    echo "Your arguments: $*"
    exit 1
fi
indir=$1
outdir=$2

# Report start of script and variables
echo "# Starting script multiqc.sh"
date
echo "# Input dir:                      $indir"
echo "# Output file:                    $outdir"
echo

# Create the output dir (with a subdir for Slurm logs)
mkdir -p "$outdir"/logs

# Run MultiQC
multiqc \
    --outdir "$outdir" \
    --interactive \
    --force \
    "$indir"

# Explanation of some of the options used:
# --interactive     => Always use interactive plots
# --force           => Overwrite existing MultiQC reports in output dir

# Report software version, end of script, and date
echo
echo "# Used MultiQC version:"
multiqc --version
echo "# Done with script multiqc.sh"
date
