#!/bin/bash
#SBATCH --partition=normal
#SBATCH --job-name=SNP_Caller_QC
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
ploidy_map="${4:-}"  # Optional 4th parameter for ploidy map
generate_report="true"

echo "=========================================="
echo "SNP Calling Pipeline with QC Reporting"
echo "=========================================="
echo ""

echo "Freebayes parameters: ${freebayes_parameters}"
if [ -f "${freebayes_parameters}" ]; then
    cat ${freebayes_parameters}
fi
echo ""

echo "BAM filter parameters: ${bam_filter_parameters}"
if [ -f "${bam_filter_parameters}" ]; then
    cat ${bam_filter_parameters}
fi
echo ""

echo "Number of chunks: ${nChunk}"
echo "Generate QC report: ${generate_report}"

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

echo ""
echo "Starting pipeline..."
echo "===================="

# Run the nextflow pipeline with QC reporting
nextflow run gcl_snp_calling/main.nf \
    -profile slurm \
    -resume \
    --bams "./data/bam/*.bam" \
    --reference "./genome/genome.fa" \
    --freebayes_config ${freebayes_parameters} \
    --bam_filter_config ${bam_filter_parameters} \
    --num_chunks ${nChunk} \
    ${PLOIDY_MAP_PARAM} \
    --generate_report ${generate_report} \
    --output_dir "data/snps" \
    --output_vcf "all_samples.vcf.gz"

# Check exit status
EXIT_STATUS=$?

if [ $EXIT_STATUS -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "Pipeline completed successfully!"
    echo "=========================================="
    
    if [ "${generate_report}" = "true" ]; then
        echo ""
        echo "📊 QC Report Summary:"
        echo "--------------------"
        
        # Display key metrics if report was generated
        if [ -f "data/snps/qc/filter_thresholds.json" ]; then
            echo ""
            echo "Filtering Recommendations:"
            python3 -c "
import json
data = json.load(open('data/snps/qc/filter_thresholds.json'))
thresh = data['recommended_thresholds']
impact = data['filtering_impact']
print(f'  • Min Depth:        {thresh[\"min_depth\"]}X')
print(f'  • Min Quality:      {thresh[\"min_quality\"]}')
print(f'  • Max Missingness:  {thresh[\"max_missingness\"]*100:.0f}%')
print(f'  • Expected retention: {impact[\"percent_passing\"]:.1f}% of variants')
"
        fi
        
        echo ""
        echo "📁 Output Files:"
        echo "  • VCF: data/snps/all_samples.vcf.gz"
        echo "  • Report: data/snps/snp_calling_report.html"
        echo "  • Filter script: data/snps/qc/filter_commands.sh"
        echo "  • QC plots: data/snps/qc/plots/"
        
        echo ""
        echo "🎯 Next Steps:"
        echo "  1. Open data/snps/snp_calling_report.html in a web browser"
        echo "  2. Review QC metrics and sample outliers"
        echo "  3. Run data/snps/qc/filter_commands.sh to apply filters"
        echo "  4. Annotate filtered variants with SnpEff/VEP"
    fi
else
    echo ""
    echo "=========================================="
    echo "Pipeline failed with exit status: $EXIT_STATUS"
    echo "=========================================="
    echo "Check the Nextflow log for details"
fi

echo ""
echo "Job completed at: $(date)"