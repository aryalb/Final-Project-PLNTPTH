#!/bin/bash
#SBATCH --account=PAS2700
#SBATCH --output=fastq-%j.out
#SBATCH --cpus-per-task=15
#SBATCH --time=5:00:00
#SBATCH --mail-type=END,FAILls
set -euo pipefail

module load sratoolkit/2.11.2

# Create output directories 
mkdir -p "results/data/Samples_FASTQ"
#mkdir -p "results/data/RefSeq_FASTA"
#mkdir -p "results/data/RefSeq_GTF"

# Accept the accession number as an argument
SRR_ID="$1"

# Download SRA file
    if ! prefetch --output-directory "results/data/Samples_FASTQ" $SRR_ID; then
    echo "Error downloading $SRR_ID..."
    fi

# Convert SRA file to FASTQ format
  if ! fasterq-dump --split-files --outdir "results/data/Samples_FASTQ" --skip-technical $SRR_ID; then
    echo "Error converting $SRR_ID to FASTQ..."
    continue
  fi

# Zip and rename FASTQ files
gzip "results/data/Samples_FASTQ/${SRR_ID}_1.fastq"
gzip "results/data/Samples_FASTQ/${SRR_ID}_2.fastq"
mv "results/data/Samples_FASTQ/${SRR_ID}_1.fastq.gz" "results/data/Samples_FASTQ/$(basename ${SRR_ID})_R1.fastq.gz"
mv "results/data/Samples_FASTQ/${SRR_ID}_2.fastq.gz" "results/data/Samples_FASTQ/$(basename ${SRR_ID})_R2.fastq.gz"