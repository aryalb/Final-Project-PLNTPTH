# Activate the conda environment
module load miniconda3/23.3.1-py310
conda activate /fs/ess/PAS0471/jelmer/conda/nextflow

# Create an environment variable for the container dir
export NXF_SINGULARITY_CACHEDIR=/fs/ess/PAS2700/containers

# Download the nf-core rnaseq pipeline files
nf-core download rnaseq \
    --revision 3.14.0 \
    --outdir software/nfc-rnaseq \
    --compress none \
    --container-system singularity \
    --container-cache-utilisation amend \
    --download-configuration

#
# Defining the pipeline outputs
outdir=results/nfc-rnaseq
workdir=/fs/scratch/PAS2700/week6/nfc-rnaseq/$USER

# Defining the pipeline inputs
samplesheet="$outdir"/nfc_samplesheet.csv
fasta=../pipeline/results/data/RefSeq_FASTA/Homo_sapiens.GRCh38.dna.primary_assembly.fa
gtf=../pipeline/results/data/RefSeq_GTF/Homo_sapiens.GRCh38.99.gtf

# Create the dir that will contain the sample sheet
mkdir -p "$outdir"

# Create the sample sheet for the nf-core pipeline
python3 software/nfc-rnaseq/3_14_0/bin/fastq_dir_to_samplesheet.py \
    --strandedness forward \
    --read1_extension "_R1.fastq.gz" \
    --read2_extension "_R2.fastq.gz" \
    ../pipeline/results/data/Samples_FASTQ \
    "$samplesheet"

# Create a config file for batch job submissions
echo "
process.executor='slurm'
process.clusterOptions='--account=PAS2700'
" > nextflow.config

# Submit the script to run the pipeline as a batch job 
sbatch scripts/nfc-rnaseq.sh "$samplesheet" "$fasta" "$gtf" "$outdir" "$workdir"

