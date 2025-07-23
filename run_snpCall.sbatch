#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=QC_Handler
#SBATCH --time=4-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5

module load miniconda3
source activate nextflow

nextflow run gcl_snp_calling/main.nf \
    -profile slurm \
    -resume \
    --bams "./data/bam/*.bam" \
    --reference "./genome/genome.fa" \
    --num_chunks 50 \
    --output_dir "data/snps" \
    --output_vcf "all_samples.vcf.gz"
