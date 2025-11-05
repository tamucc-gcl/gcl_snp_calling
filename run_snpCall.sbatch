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

# Parse command line arguments
freebayes_parameters="${1:-gcl_snp_calling/freebayes_parameters.json}"
bam_filter_parameters="${2:-gcl_snp_calling/bam_filter_parameters.json}"
nChunk="${3:-50}"
ploidy_map="${4:-}"  # NEW: Optional 4th parameter for ploidy map

echo "Freebayes parameters: ${freebayes_parameters}"
cat ${freebayes_parameters}
echo ""

echo "BAM filter parameters: ${bam_filter_parameters}"
cat ${bam_filter_parameters}
echo ""

echo "Number of chunks: ${nChunk}"

# Check if ploidy map was provided
if [ -n "${ploidy_map}" ] && [ -f "${ploidy_map}" ]; then
    echo "Ploidy map: ${ploidy_map}"
    echo "Ploidy map contents:"
    cat ${ploidy_map}
    echo ""
    PLOIDY_MAP_PARAM="--ploidy_map ${ploidy_map}"
else
    echo "No ploidy map provided - using global ploidy from config or default"
    PLOIDY_MAP_PARAM=""
fi

nextflow run gcl_snp_calling/main.nf \
    -profile slurm \
    -resume \
    --run_qc \
    --bams "./data/bam/*.bam" \
    --reference "./genome/genome.fa" \
    --freebayes_config ${freebayes_parameters} \
    --bam_filter_config ${bam_filter_parameters} \
    --num_chunks ${nChunk} \
    ${PLOIDY_MAP_PARAM} \
    --output_dir "results" \
    --output_vcf "raw_variants.vcf.gz"
