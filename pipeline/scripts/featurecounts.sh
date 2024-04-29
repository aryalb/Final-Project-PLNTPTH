#!/bin/bash
#SBATCH --account=PAS2700
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-featurecounts-%j.out
set -euo pipefail

# Load the software
module load miniconda3/24.1.2-py310
conda activate /fs/ess/PAS0471/jelmer/conda/subread

# Process the command-line arguments
if [[ ! "$#" -eq 4 ]]; then
    echo "Error: You provided $# arguments, while 4 are required."
    echo "Usage: featurecounts.sh <input-dir> <annotation-gtf> <strandedness> <outfile>"
    echo "Example: sbatch featurecounts.sh results/star data/ref/my.gtf 2 results/featurecounts/counts.tsv"
    echo "Your arguments: $*"
    exit 1
fi
indir=$1
gtf=$2
strandedness=$3
outfile=$4

# Make sure the nr of threads/cores/cpus is however many the job has
threads=$SLURM_CPUS_PER_TASK

# Report start of script and variables
echo "# Starting script featurecounts.sh"
date
echo "# Input dir:                      $indir"
echo "# Input annotation GTF:           $gtf"
echo "# Strandedness:                   $strandedness"
echo "# Output file:                    $outfile"
echo "# Number of threads:              $threads"
echo

# Create the output dir (with a subdir for Slurm logs)
outdir=$(dirname "$outfile")
mkdir -p "$outdir"/logs

# Run featureCounts
featureCounts \
    -a "$gtf" \
    -o "$outfile" \
    -s "$strandedness" \
    -p \
    --countReadPairs \
    -C \
    -M \
    -T "$threads" \
    "$indir"/*bam

# Explanation of some of the options used:
# -s 2              => Reverse-stranded library like TruSeq
# -p                => paired-end reads
# --countReadPairs  => Count read pairs, not reads
# -C                => Don't count pairs with discordant mates
# -M                => Count multi-mapping reads

# Report software version, end of script, and date
echo
echo "# Used featureCounts/subread version:"
featureCounts -v
echo "# Done with script featurecounts.sh"
date
