#!/bin/bash
#SBATCH --account=PAS2700
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-star_index-%j.out
set -euo pipefail

# Load the software
module load miniconda3/24.1.2-py310
conda activate /fs/ess/PAS0471/jelmer/conda/star

# Process the command-line arguments
if [[ ! "$#" -eq 3 ]]; then
    echo "Error: You provided $# arguments, while 3 are required."
    echo "Usage: star_index.sh <assembly-fasta> <annotation-gtf> <outdir>"
    echo "Example: sbatch star_index.sh data/ref/my.fna data/ref/my.gtf results/star/index"
    echo "Your arguments: $*"
    exit 1
fi
fasta=$1
gtf=$2
outdir=$3

# Make sure the nr of threads/cores/cpus is however many the job has
threads=$SLURM_CPUS_PER_TASK

# Report start of script and variables
echo "# Starting script star_index.sh"
date
echo "# Input assembly FASTA:           $fasta"
echo "# Input annotation GTF:           $gtf"
echo "# Output dir:                     $outdir"
echo "# Number of threads:              $threads"
echo

# Create the output dir (with a subdir for Slurm logs)
mkdir -p "$outdir"/logs

# Run STAR
# (Note: keeping this command as simple as possible - if using STAR in your research,
#  then look into the --sjdbOverhang and --genomeSAindexNbases options as well.)
STAR \
    --runMode genomeGenerate \
    --genomeFastaFiles "$fasta" \
    --sjdbGTFfile "$gtf" \
    --genomeDir "$outdir" \
    --runThreadN "$threads"

# Explanation of some of the options used:
# --runMode genomeGenerate  => Instead of aligning, the "run mode" is to generate a genome index

# Report software version, end of script, and date
echo
echo "# Used STAR version:"
STAR --version
echo "# Done with script star_index.sh"
date

