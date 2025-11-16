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

# Script usage function
usage() {
    cat << EOF
================================================================================
USAGE: sbatch run_snpCall.sbatch [OPTIONS]

DESCRIPTION:
    Submit SNP calling pipeline using either FreeBayes or ANGSD

OPTIONS:
    -g, --genotyper <STRING>      Genotyper to use: 'freebayes' or 'angsd' (default: freebayes)
    -c, --config <FILE>           Configuration file (freebayes_parameters.json or angsd_parameters.json)
    -b, --bam-filter <FILE>       BAM filter parameters (default: gcl_snp_calling/bam_filter_parameters.json)
    -n, --chunks <INT>            Number of chunks to split genome (default: 50)
    -p, --ploidy-map <FILE>       Ploidy map file (FreeBayes only, optional)
    -s, --sites <FILE>            Sites file for targeted analysis (ANGSD only, optional)
    -o, --output <STRING>         Output VCF name (used as base name for ANGSD outputs)
    -h, --help                    Show this help message

EXAMPLES:
    # FreeBayes with default settings
    sbatch run_snpCall.sbatch

    # FreeBayes with custom config and ploidy map
    sbatch run_snpCall.sbatch -g freebayes -c freebayes_params.json -p ploidy_map.txt

    # ANGSD for individual samples (uses default output name)
    sbatch run_snpCall.sbatch -g angsd -c angsd_parameters.json

    # ANGSD for pooled samples with custom output name
    sbatch run_snpCall.sbatch -g angsd -c angsd_pooled_parameters.json -n 100 -o pooled.vcf.gz

    # ANGSD with sites file for targeted analysis
    sbatch run_snpCall.sbatch -g angsd -c angsd_parameters.json -s high_quality_sites.txt

    # Both genotypers with same output base name
    sbatch run_snpCall.sbatch -g freebayes -o my_project.vcf.gz  # Creates: my_project.vcf.gz
    sbatch run_snpCall.sbatch -g angsd -o my_project.vcf.gz      # Creates: my_project.beagle.gz, etc.

LEGACY USAGE (still supported):
    sbatch run_snpCall.sbatch <freebayes_params> <bam_filter_params> <nChunks> [ploidy_map]
================================================================================
EOF
    exit 0
}

# Default values
GENOTYPER="freebayes"
CONFIG_FILE=""
BAM_FILTER="gcl_snp_calling/bam_filter_parameters.json"
N_CHUNKS="50"
PLOIDY_MAP=""
SITES_FILE=""
OUTPUT_VCF="raw_variants.vcf.gz"  # Default output name for both genotypers

# Check if using legacy format (positional arguments)
if [[ $# -gt 0 && "$1" != "-"* ]]; then
    echo "Detected legacy command format - converting to new format"

    # Legacy format: <freebayes_params> <bam_filter_params> <nChunks> [ploidy_map]
    CONFIG_FILE="${1:-gcl_snp_calling/freebayes_parameters.json}"
    BAM_FILTER="${2:-gcl_snp_calling/bam_filter_parameters.json}"
    N_CHUNKS="${3:-50}"
    PLOIDY_MAP="${4:-}"
    GENOTYPER="freebayes"
    OUTPUT_VCF="raw_variants.vcf.gz"  # Use default for legacy format

    echo "Legacy parameters interpreted as:"
    echo "  FreeBayes config: ${CONFIG_FILE}"
    echo "  BAM filter: ${BAM_FILTER}"
    echo "  Chunks: ${N_CHUNKS}"
    if [ -n "${PLOIDY_MAP}" ]; then
        echo "  Ploidy map: ${PLOIDY_MAP}"
    fi
else
    # Parse command line arguments using getopt
    TEMP=$(getopt -o g:c:b:n:p:s:o:h \
                  --long genotyper:,config:,bam-filter:,chunks:,ploidy-map:,sites:,output:,help \
                  -n 'run_snpCall.sbatch' -- "$@")

    if [ $? != 0 ]; then
        echo "Error parsing arguments. Use -h for help." >&2
        exit 1
    fi

    eval set -- "$TEMP"

    while true; do
        case "$1" in
            -g|--genotyper)
                GENOTYPER="$2"
                shift 2
                ;;
            -c|--config)
                CONFIG_FILE="$2"
                shift 2
                ;;
            -b|--bam-filter)
                BAM_FILTER="$2"
                shift 2
                ;;
            -n|--chunks)
                N_CHUNKS="$2"
                shift 2
                ;;
            -p|--ploidy-map)
                PLOIDY_MAP="$2"
                shift 2
                ;;
            -s|--sites)
                SITES_FILE="$2"
                shift 2
                ;;
            -o|--output)
                OUTPUT_VCF="$2"
                shift 2
                ;;
            -h|--help)
                usage
                ;;
            --)
                shift
                break
                ;;
            *)
                echo "Internal error!"
                exit 1
                ;;
        esac
    done
fi

# Validate genotyper choice
if [[ "${GENOTYPER}" != "freebayes" && "${GENOTYPER}" != "angsd" ]]; then
    echo "ERROR: Invalid genotyper '${GENOTYPER}'. Must be 'freebayes' or 'angsd'"
    exit 1
fi

# Set default config files if not specified
if [ -z "${CONFIG_FILE}" ]; then
    if [ "${GENOTYPER}" = "freebayes" ]; then
        CONFIG_FILE="gcl_snp_calling/freebayes_parameters.json"
    else
        CONFIG_FILE="gcl_snp_calling/angsd_parameters.json"
    fi
fi

# Derive output prefix for ANGSD from OUTPUT_VCF
OUTPUT_PREFIX=$(echo "${OUTPUT_VCF}" | sed 's/\.\(vcf\|VCF\)\(\.gz\)\?$//')

# Print configuration summary
echo "================================================================"
echo "SNP Calling Pipeline Configuration"
echo "================================================================"
echo "Genotyper:         ${GENOTYPER}"
echo "Config file:       ${CONFIG_FILE}"
echo "BAM filter config: ${BAM_FILTER}"
echo "Number of chunks:  ${N_CHUNKS}"
echo "Output VCF name:   ${OUTPUT_VCF}"
if [ "${GENOTYPER}" = "angsd" ]; then
    echo "Output prefix:     ${OUTPUT_PREFIX} (derived from VCF name)"
fi

if [ -f "${CONFIG_FILE}" ]; then
    echo ""
    echo "Configuration file contents:"
    echo "----------------------------"
    cat "${CONFIG_FILE}"
    echo "----------------------------"
fi

if [ -f "${BAM_FILTER}" ]; then
    echo ""
    echo "BAM filter configuration:"
    echo "------------------------"
    cat "${BAM_FILTER}"
    echo "------------------------"
fi

# Build Nextflow command
CMD="nextflow run gcl_snp_calling/main.nf"
CMD="${CMD} -profile slurm"
CMD="${CMD} -resume"

BAM_PATTERN="./data/bam/*.bam"
CMD="${CMD} --bams '${BAM_PATTERN}'"
#CMD="${CMD} --bams ./data/bam/*.bam"
CMD="${CMD} --reference ./genome/genome.fa"
CMD="${CMD} --genotyper ${GENOTYPER}"
CMD="${CMD} --num_chunks ${N_CHUNKS}"
CMD="${CMD} --output_dir results"
CMD="${CMD} --output_vcf ${OUTPUT_VCF}"  # Used for both genotypers

# Add BAM filter config if it exists
if [ -f "${BAM_FILTER}" ]; then
    CMD="${CMD} --bam_filter_config ${BAM_FILTER}"
fi

# Add genotyper-specific parameters
if [ "${GENOTYPER}" = "freebayes" ]; then
    echo ""
    echo "FreeBayes-specific settings:"
    echo "----------------------------"

    # Add FreeBayes config
    if [ -f "${CONFIG_FILE}" ]; then
        CMD="${CMD} --freebayes_config ${CONFIG_FILE}"
    fi

    echo "Output VCF: ${OUTPUT_VCF}"

    # Add ploidy map if provided
    if [ -n "${PLOIDY_MAP}" ] && [ -f "${PLOIDY_MAP}" ]; then
        echo "Ploidy map: ${PLOIDY_MAP}"
        echo ""
        echo "Ploidy map contents:"
        cat "${PLOIDY_MAP}"
        CMD="${CMD} --ploidy_map ${PLOIDY_MAP}"
    else
        echo "Ploidy: Using global setting from config or FreeBayes default (diploid)"
    fi

elif [ "${GENOTYPER}" = "angsd" ]; then
    echo ""
    echo "ANGSD-specific settings:"
    echo "------------------------"

    # Add ANGSD config
    if [ -f "${CONFIG_FILE}" ]; then
        CMD="${CMD} --angsd_config ${CONFIG_FILE}"

        # Check if config specifies pooled mode
        IS_POOLED=$(grep -o '"is_pooled"[[:space:]]*:[[:space:]]*true' "${CONFIG_FILE}" || echo "")
        if [ -n "${IS_POOLED}" ]; then
            echo "Mode: Pooled sequencing"
        else
            echo "Mode: Individual samples"
        fi
    fi

    echo "Output base name: ${OUTPUT_PREFIX}"

    # Add sites file if provided
    if [ -n "${SITES_FILE}" ] && [ -f "${SITES_FILE}" ]; then
        echo "Sites file: ${SITES_FILE}"
        echo "Number of sites: $(wc -l < ${SITES_FILE})"
        CMD="${CMD} --sites_file ${SITES_FILE}"
    else
        echo "Sites: Will discover from data"
    fi

    # Show expected outputs
    echo ""
    echo "Expected output files:"
    echo "  - ${OUTPUT_PREFIX}.beagle.gz (genotype likelihoods)"
    echo "  - ${OUTPUT_PREFIX}.mafs.gz (allele frequencies)"
    echo "  - ${OUTPUT_PREFIX}.vcf.gz (VCF format)"
    echo "  - ${OUTPUT_PREFIX}.sites (positions file)"
    echo "  - ${OUTPUT_PREFIX}_summary.txt (run summary)" "Sites: Will discover from data"
#    fi

    # Show expected outputs
    echo ""
    echo "Expected output files:"
    if [[ "${OUTPUT_FORMAT}" == "all" || "${OUTPUT_FORMAT}" == "beagle" ]]; then
        echo "  - ${OUTPUT_PREFIX}.beagle.gz (genotype likelihoods)"
        echo "  - ${OUTPUT_PREFIX}.mafs.gz (allele frequencies)"
    fi
    if [[ "${OUTPUT_FORMAT}" == "all" || "${OUTPUT_FORMAT}" == "vcf" ]]; then
        echo "  - ${OUTPUT_PREFIX}.vcf.gz (VCF format)"
    fi
    echo "  - ${OUTPUT_PREFIX}.sites (positions file)"
    echo "  - ${OUTPUT_PREFIX}_summary.txt (run summary)"
fi

echo ""
echo "================================================================"
echo "Executing Nextflow Pipeline"
echo "================================================================"
echo "Command:"
echo "${CMD}"
echo "================================================================"
echo ""

# Execute the pipeline
${CMD}

# Check exit status
EXIT_CODE=$?

if [ ${EXIT_CODE} -eq 0 ]; then
    echo ""
    echo "================================================================"
    echo "Pipeline completed successfully!"
    echo "================================================================"
    echo "Results location: ./results/"

    if [ "${GENOTYPER}" = "freebayes" ]; then
        echo "Main output: ./results/${OUTPUT_VCF}"
    else
        echo "Main outputs:"
        echo "  - ./results/${OUTPUT_PREFIX}.beagle.gz"
        echo "  - ./results/${OUTPUT_PREFIX}.mafs.gz"
        echo "  - ./results/${OUTPUT_PREFIX}.vcf.gz"
    fi

    echo ""
    echo "QC Reports:"
    echo "  - ./results/multiqc_reports/multiqc_raw_bams.html"
    if [ -f "${BAM_FILTER}" ]; then
        echo "  - ./results/multiqc_reports/multiqc_filtered_bams.html"
    fi
    echo "================================================================"
else
    echo ""
    echo "================================================================"
    echo "Pipeline failed with exit code: ${EXIT_CODE}"
    echo "Check the Nextflow logs for details"
    echo "================================================================"
    exit ${EXIT_CODE}
fi
