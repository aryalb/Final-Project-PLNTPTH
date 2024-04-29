#!/bin/bash
#SBATCH --account=PAS2700
#SBATCH --time=6:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm-nfc_rnaseq-%j.out
set -euo pipefail

# Settings and constants
WORKFLOW_DIR=software/nfc-rnaseq/3_14_0

# Load the Nextflow Conda environment
module load miniconda3/23.3.1-py310
conda activate /fs/ess/PAS0471/jelmer/conda/nextflow
export NXF_SINGULARITY_CACHEDIR=/fs/ess/PAS0471/containers

# Process command-line arguments
if [[ ! "$#" -eq 5 ]]; then
    echo "Error: wrong number of arguments"
    echo "You provided $# arguments, while 5 are required."
    echo "Usage: nfc-rnaseq.sh <samplesheet> <FASTA> <GTF> <outdir> <workdir>"
    exit 1
fi
samplesheet=$1
fasta=$2
gtf=$3
outdir=$4
workdir=$5

# Report
echo "Starting script nfc-rnaseq.sh"
date
echo "Samplesheet:          $samplesheet"
echo "Reference FASTA:      $fasta"
echo "Reference GTF:        $gtf"
echo "Output dir:           $outdir"
echo "Nextflow work dir:    $workdir"
echo

# Create the output dirs
mkdir -p "$outdir" "$workdir"

# Run the workflow
nextflow run "$WORKFLOW_DIR" \
    --input "$samplesheet" \
    --fasta "$fasta" \
    --gtf "$gtf" \
    --outdir "$outdir" \
    --remove_ribo_rna \
    -work-dir "$workdir" \
    -profile singularity \
    -ansi-log false \
    -resume

# Report
echo "Done with script nfc-rnaseq.sh"
date

