# 2024-04-15, Bikash Aryal
# - A runner script to run a simple RNA-seq analysis, creating a gene count table
#   from RNA-seq reads and a reference genome assembly and annotation.
# - The steps are read QC and trimming, alignment, and quantification.

## Step 0: Download samples in fastq format with their SRR IDs.

# Define SRR IDs
SRR_IDs=("SRR923920" "SRR923921" "SRR923922" "SRR923923" "SRR923924" "SRR923925" "SRR923926" "SRR923927" "SRR923928" "SRR923929" "SRR923930" "SRR923931" "SRR923932" "SRR923933")

# Looping over each accession and submmiting a separate job
for SRR_ID in "${SRR_IDs[@]}"; do
  sbatch scripts/data.sh "$SRR_ID"
done

## Step 1: Download Reference genome 
# Download reference genome in FASTA format directly to RefSeq_FASTA directory
wget -q -P results/data/RefSeq_FASTA ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip results/data/RefSeq_FASTA/*.gz

# Download reference genome in GTF format directly to RefSeq_GTF directory
wget -q -P results/data/RefSeq_GTF ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz 
gunzip results/data/RefSeq_GTF/*.gz

## Step 02: Setting up for analysis
# A) Define the inputs
fastq_dir=results/data/Samples_FASTQ
fasta=results/data/RefSeq_FASTA/Homo_sapiens.GRCh38.dna.primary_assembly.fa
gtf=results/data/RefSeq_GTF/Homo_sapiens.GRCh38.99.gtf

# B) Define settings
strandedness=1       # 1 means forward-stranded (used with featureCounts)
# C) Define the outputs
count_table=results/featurecounts/counts.tsv

# Step 03: QC the reads with FastQC
for fastq_file in "$fastq_dir"/*fastq.gz; do
    sbatch scripts/fastqc.sh "$fastq_file" results/fastqc
done

# Step 04: Trim the reads with TrimGalore
for R1 in "$fastq_dir"/*_R1.fastq.gz; do
    sbatch scripts/trimgalore.sh "$R1" results/trimgalore
done

# Step 05: Align the reads with STAR
# A) Index the reference genomes
sbatch scripts/star_index.sh "$fasta" "$gtf" results/star/index

# B) Align the reads to the index
for R1 in results/trimgalore/*R1.fastq.gz; do
    sbatch scripts/star_align.sh "$R1" results/star/index "$gtf" results/star
done

# Step 06: Create a gene count table with featureCounts
sbatch scripts/featurecounts.sh results/star "$gtf" "$strandedness" results/featurecounts/counts.tsv

# Step 07: QC summaries with MultiQC
sbatch scripts/multiqc.sh results results/multiqc


