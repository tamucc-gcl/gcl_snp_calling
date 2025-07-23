// modules/freebayes_chunk.nf - Fixed for multiple BAM files

process FREEBAYES_CHUNK {
    tag "chunk_${chunk_id}"
    
    input:
    tuple val(chunk_id), val(regions_string), path(reference), path(bams), path(bam_indices), path(config_file)
    
    output:
    tuple val(chunk_id), path("chunk_${chunk_id}.vcf.gz")
    
    script:
    // Create the --bam arguments as a Groovy string that will be properly interpolated
    def bam_args = bams.collect{ "--bam ${it.name}" }.join(' ')
    def config_name = config_file.name
    """
    # Debug: Show all BAM files
    echo "Processing chunk ${chunk_id}"
    echo "Number of BAM files: ${bams.size()}"
    echo "BAM files present in work directory:"
    ls -la *.bam 2>/dev/null || echo "No BAM files found!"
    
    # Ensure BAM indices exist
    for bam in *.bam; do
        if [ -f "\$bam" ]; then
            if [ ! -f "\${bam}.bai" ] && [ ! -f "\${bam%.bam}.bai" ]; then
                echo "Creating index for \$bam"
                samtools index "\$bam"
            fi
        fi
    done
    
    echo "Regions: ${regions_string}"
    
    # Parse regions string and create BED file for this chunk
    echo "${regions_string}" | tr ',' '\\n' > regions_list.txt
    
    # Convert regions to BED format for freebayes --targets
    # Create BED file WITHOUT comment lines (freebayes doesn't like them)
    touch chunk_targets.bed
    
    while read region; do
        if [ ! -z "\$region" ]; then
            # Parse region format: chr:start-end
            chrom=\$(echo \$region | cut -d':' -f1)
            positions=\$(echo \$region | cut -d':' -f2)
            start=\$(echo \$positions | cut -d'-' -f1)
            end=\$(echo \$positions | cut -d'-' -f2)
            
            # Convert to 0-based BED format
            bed_start=\$((start - 1))
            
            echo -e "\$chrom\\t\$bed_start\\t\$end" >> chunk_targets.bed
            echo "  Added region: \$chrom:\$start-\$end (BED: \$chrom:\$bed_start-\$end)"
        fi
    done < regions_list.txt
    
    echo "Created BED file for chunk ${chunk_id}:"
    cat chunk_targets.bed
    
    # Verify BED file is not empty
    if [ ! -s chunk_targets.bed ]; then
        echo "ERROR: BED file is empty for chunk ${chunk_id}"
        exit 1
    fi
    
    # Parse freebayes config if provided
    FREEBAYES_OPTS=""
    
    # Check if config file exists and is not the placeholder
    if [ "${config_name}" != "NO_CONFIG" ] && [ -f "${config_file}" ] && [[ "${config_file}" == *.json ]]; then
        echo "Loading freebayes configuration from ${config_file}"
        
        # Parse JSON config using Python (available in most systems)
        FREEBAYES_OPTS=\$(python3 << 'PYTHON_SCRIPT'
import json
import sys

# Read the config file
with open('${config_file}', 'r') as f:
    config = json.load(f)

# Build freebayes options
opts = []

# Algorithm parameters
if 'algorithm_parameters' in config:
    params = config['algorithm_parameters']
    
    if params.get('min_mapping_quality') is not None:
        opts.append(f"--min-mapping-quality {params['min_mapping_quality']}")
    if params.get('min_base_quality') is not None:
        opts.append(f"--min-base-quality {params['min_base_quality']}")
    if params.get('min_supporting_allele_qsum') is not None:
        opts.append(f"--min-supporting-allele-qsum {params['min_supporting_allele_qsum']}")
    if params.get('min_supporting_mapping_qsum') is not None:
        opts.append(f"--min-supporting-mapping-qsum {params['min_supporting_mapping_qsum']}")
    if params.get('mismatch_base_quality_threshold') is not None:
        opts.append(f"--mismatch-base-quality-threshold {params['mismatch_base_quality_threshold']}")
    if params.get('read_mismatch_limit') is not None:
        opts.append(f"--read-mismatch-limit {params['read_mismatch_limit']}")
    if params.get('read_max_mismatch_fraction') is not None:
        opts.append(f"--read-max-mismatch-fraction {params['read_max_mismatch_fraction']}")
    if params.get('read_snp_limit') is not None:
        opts.append(f"--read-snp-limit {params['read_snp_limit']}")
    if params.get('read_indel_limit') is not None:
        opts.append(f"--read-indel-limit {params['read_indel_limit']}")
    if params.get('indel_exclusion_window') is not None:
        opts.append(f"--indel-exclusion-window {params['indel_exclusion_window']}")
    if params.get('min_repeat_entropy') is not None:
        opts.append(f"--min-repeat-entropy {params['min_repeat_entropy']}")
    if params.get('min_alternate_fraction') is not None:
        opts.append(f"--min-alternate-fraction {params['min_alternate_fraction']}")
    if params.get('min_alternate_count') is not None:
        opts.append(f"--min-alternate-count {params['min_alternate_count']}")
    if params.get('min_alternate_qsum') is not None:
        opts.append(f"--min-alternate-qsum {params['min_alternate_qsum']}")
    if params.get('min_alternate_total') is not None:
        opts.append(f"--min-alternate-total {params['min_alternate_total']}")
    if params.get('min_coverage') is not None:
        opts.append(f"--min-coverage {params['min_coverage']}")
    if params.get('max_coverage') is not None:
        opts.append(f"--max-coverage {params['max_coverage']}")
    if params.get('no_partial_observations'):
        opts.append("--no-partial-observations")
    if params.get('ploidy') is not None:
        opts.append(f"--ploidy {params['ploidy']}")
    if params.get('pooled_discrete'):
        opts.append("--pooled-discrete")
    if params.get('pooled_continuous'):
        opts.append("--pooled-continuous")
    if params.get('use_duplicate_reads'):
        opts.append("--use-duplicate-reads")
    if params.get('no_mnps'):
        opts.append("--no-mnps")
    if params.get('no_complex'):
        opts.append("--no-complex")
    if params.get('no_snps'):
        opts.append("--no-snps")
    if params.get('no_indels'):
        opts.append("--no-indels")
    if params.get('max_complex_gap') is not None:
        opts.append(f"--max-complex-gap {params['max_complex_gap']}")
    if params.get('haplotype_length') is not None:
        opts.append(f"--haplotype-length {params['haplotype_length']}")
    if params.get('min_repeat_length') is not None:
        opts.append(f"--min-repeat-length {params['min_repeat_length']}")
    if params.get('min_repeat_entropy_for_detection') is not None:
        opts.append(f"--min-repeat-entropy {params['min_repeat_entropy_for_detection']}")

# Population model parameters
if 'population_model_parameters' in config:
    params = config['population_model_parameters']
    
    if params.get('theta') is not None:
        opts.append(f"--theta {params['theta']}")
    if params.get('posterior_integration_limits') is not None:
        opts.append(f"--posterior-integration-limits {params['posterior_integration_limits']}")
    if params.get('use_reference_allele'):
        opts.append("--use-reference-allele")
    if params.get('gvcf'):
        opts.append("--gvcf")
    if params.get('gvcf_chunk') is not None:
        opts.append(f"--gvcf-chunk {params['gvcf_chunk']}")

# Genotype likelihood parameters
if 'genotype_likelihood_parameters' in config:
    params = config['genotype_likelihood_parameters']
    
    if params.get('base_quality_cap') is not None:
        opts.append(f"--base-quality-cap {params['base_quality_cap']}")
    if params.get('prob_contamination') is not None:
        opts.append(f"--prob-contamination {params['prob_contamination']}")
    if params.get('legacy_gls'):
        opts.append("--legacy-gls")
    if params.get('contamination_estimates') is not None:
        opts.append(f"--contamination-estimates {params['contamination_estimates']}")

# Reporting parameters
if 'reporting_parameters' in config:
    params = config['reporting_parameters']
    
    if params.get('genotype_qualities'):
        opts.append("--genotype-qualities")
    if params.get('report_genotype_likelihood_max'):
        opts.append("--report-genotype-likelihood-max")
    if params.get('genotyping_max_iterations') is not None:
        opts.append(f"--genotyping-max-iterations {params['genotyping_max_iterations']}")
    if params.get('genotyping_max_banddepth') is not None:
        opts.append(f"--genotyping-max-banddepth {params['genotyping_max_banddepth']}")
    if params.get('exclude_unobserved_genotypes'):
        opts.append("--exclude-unobserved-genotypes")
    if params.get('genotype_variant_threshold') is not None:
        opts.append(f"--genotype-variant-threshold {params['genotype_variant_threshold']}")
    if params.get('use_mapping_quality'):
        opts.append("--use-mapping-quality")
    if params.get('harmonic_indel_quality'):
        opts.append("--harmonic-indel-quality")
    if params.get('read_dependence_factor') is not None:
        opts.append(f"--read-dependence-factor {params['read_dependence_factor']}")

# Additional options
if 'additional_options' in config:
    params = config['additional_options']
    
    if params.get('debug'):
        opts.append("--debug")
    if params.get('dd'):
        opts.append("--dd")

# Print the options
print(' '.join(opts))
PYTHON_SCRIPT
)
        
        echo "Freebayes options: \$FREEBAYES_OPTS"
    else
        echo "No freebayes config provided, using default parameters"
    fi
    
    # Run freebayes with targets file
    echo "Running freebayes on chunk ${chunk_id}..."
    echo "Using BAM arguments: ${bam_args}"
    
    freebayes \\
        --fasta-reference ${reference} \\
        --targets chunk_targets.bed \\
        ${bam_args} \\
        \$FREEBAYES_OPTS \\
        --vcf chunk_${chunk_id}.vcf
    
    # Check if freebayes produced output
    if [ -s chunk_${chunk_id}.vcf ]; then
        echo "Freebayes completed successfully for chunk ${chunk_id}"
        
        # Count variants found
        variant_count=\$(grep -v '^#' chunk_${chunk_id}.vcf | wc -l)
        echo "Found \$variant_count variants in chunk ${chunk_id}"
        
        # Show sample names in VCF
        echo "Samples in VCF:"
        grep "^#CHROM" chunk_${chunk_id}.vcf | cut -f10- | tr '\\t' '\\n' | head -10
        echo "Total samples: \$(grep "^#CHROM" chunk_${chunk_id}.vcf | cut -f10- | tr '\\t' '\\n' | wc -l)"
        
        # Compress and index the VCF
        bgzip -c chunk_${chunk_id}.vcf > chunk_${chunk_id}.vcf.gz
        tabix -p vcf chunk_${chunk_id}.vcf.gz
        
        echo "Successfully compressed and indexed chunk ${chunk_id}"
    else
        echo "No variants found in chunk ${chunk_id}, creating empty VCF"
        
        # Create minimal VCF header for empty chunks
        cat > chunk_${chunk_id}.vcf << 'VCF_HEADER'
##fileformat=VCFv4.2
##reference=${reference}
##source=freebayes
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
VCF_HEADER
        
        bgzip -c chunk_${chunk_id}.vcf > chunk_${chunk_id}.vcf.gz
        tabix -p vcf chunk_${chunk_id}.vcf.gz
    fi
    
    # Clean up intermediate files
    rm -f chunk_${chunk_id}.vcf regions_list.txt chunk_targets.bed
    
    echo "Chunk ${chunk_id} processing complete"
    """
}