#!/bin/bash
#SBATCH --account=PAS2700
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=FAIL
#SBATCH --output=slurm-star_align-%j.out
set -euo pipefail

# Load the software
module load miniconda3/24.1.2-py310
conda activate /fs/ess/PAS0471/jelmer/conda/star

# Process the command-line arguments
if [[ ! "$#" -eq 4 ]]; then
    echo "Error: You provided $# arguments, while 4 are required."
    echo "Usage: star_align.sh <R1-FASTQ> <genome-index-dir> <annotation-gtf> <outdir>"
    echo "Example: sbatch star_align.sh data/fastq/A_R1.fastq.gz results/star/index data/ref/my.gtf results/star"
    echo "Your arguments: $*"
    exit 1
fi
R1_in=$1
index_dir=$2
gtf=$3
outdir=$4

# Make sure the nr of threads/cores/cpus is however many the job has
threads=$SLURM_CPUS_PER_TASK

# Infer the R2 FASTQ file name
R2_in=${R1_in/_R1/_R2}

# Infer the "sample ID" - we need this for the output file specification
sample_id=$(basename "$R1_in" _R1.fastq.gz)

# Report start of script and variables
echo "# Starting script star_align.sh"
date
echo "# Input R1 FASTQ file:            $R1_in"
echo "# Input R2 FASTQ file:            $R2_in"
echo "# Input genome index:             $index_dir"
echo "# Input annotation GTF:           $gtf"
echo "# Output dir:                     $outdir"
echo "# Number of threads:              $threads"
echo

# Create the output dir (with a subdir for Slurm logs)
mkdir -p "$outdir"/logs

# Run STAR
# (Note: keeping this command as simple as possible - if using STAR in your research,
#  then look into, e.g., the --alignIntronMin, --alignIntronMax and --outFilterMultimapNmax
#  options as well.)
STAR \
    --readFilesIn "$R1_in" "$R2_in" \
    --genomeDir "$index_dir" \
    --sjdbGTFfile "$gtf" \
    --runThreadN "$threads" \
    --readFilesCommand zcat \
    --outFileNamePrefix "$outdir"/"$sample_id"_ \
    --outSAMtype BAM SortedByCoordinate

# Explanation of some of the options used:
# --readFilesCommand zcat                     => Tell STAR that the FASTQ files are gzipped
# --outSAMtype BAM SortedByCoordinate         => Request a sorted BAM file as output (instead of unsorted SAM)
# --outFileNamePrefix "$outdir"/"$sample_id"_ => Specify not just an outdir but a "sample ID" prefix, otherwise
#                                                the BAM file would have a generic file name that does not ID the file/sample 

# Report software version, end of script, and date
echo
echo "# Used STAR version:"
STAR --version
echo "# Done with script star_align.sh"
date
