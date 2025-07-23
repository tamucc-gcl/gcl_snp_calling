#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=SNP_Caller
#SBATCH --time=4-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5

#module load miniconda3
#source activate nextflow
module load nextflow/23.10.1
module load miniconda3

freebayes_parameters="${1:-gcl_snp_calling/freebayes_parameters.json}"
nChunk="{$2:-50}"

echo ${freebayes_parameters}
cat ${freebayes_parameters}
echo ${nChunk}


nextflow run gcl_snp_calling/main.nf \
    -profile slurm \
    -resume \
    --bams "./data/bam/*.bam" \
    --reference "./genome/genome.fa" \
    --freebayes_config ${freebayes_parameters} \
    --num_chunks ${nChunk} \
    --output_dir "data/snps" \
    --output_vcf "all_samples.vcf.gz"
