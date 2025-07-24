// modules/freebayes_chunk.nf - COMPLETELY REWRITTEN for proper input handling

process FREEBAYES_CHUNK {
    tag "chunk_${chunk_id}"
    
    input:
    tuple val(chunk_id), val(regions_string), path(reference_fa), path(reference_fai), path(bam_files), path(bai_files), path(config_file)
    
    output:
    tuple val(chunk_id), path("chunk_${chunk_id}.vcf.gz")
    
    script:
    def config_name = config_file.name
    def is_config_provided = (config_name != "NO_CONFIG_FILE")
    """
    # Debug information
    echo "Processing chunk ${chunk_id}"
    echo "Reference FA: ${reference_fa}"
    echo "Reference FAI: ${reference_fai}"
    echo "Config file: ${config_file} (name: ${config_name})"
    echo "Is config provided: ${is_config_provided}"
    echo "Regions: ${regions_string}"
    
    # Count and list BAM files
    echo "BAM files received:"
    ls -la *.bam 2>/dev/null || echo "No BAM files found via ls"
    
    # Count actual files
    BAM_COUNT=0
    for f in *.bam; do
        if [[ -f "\$f" ]]; then
            echo "Found BAM: \$f"
            BAM_COUNT=\$((BAM_COUNT + 1))
        fi
    done
    echo "Total BAM files found: \$BAM_COUNT"
    
    if [ \$BAM_COUNT -eq 0 ]; then
        echo "ERROR: No BAM files found!"
        exit 1
    fi
    
    # Ensure BAM indices exist
    for bam in *.bam; do
        if [[ -f "\$bam" ]]; then
            if [[ ! -f "\${bam}.bai" ]] && [[ ! -f "\${bam%.bam}.bai" ]]; then
                echo "Creating index for \$bam"
                samtools index "\$bam"
            fi
        fi
    done
    
    # Create BED file for this chunk
    echo "${regions_string}" | tr ',' '\\n' > regions_list.txt
    
    touch chunk_targets.bed
    while IFS= read -r region; do
        if [[ -n "\$region" ]]; then
            # Parse region format: chr:start-end
            chrom=\$(echo "\$region" | cut -d':' -f1)
            positions=\$(echo "\$region" | cut -d':' -f2)
            start=\$(echo "\$positions" | cut -d'-' -f1)
            end=\$(echo "\$positions" | cut -d'-' -f2)
            
            # Convert to 0-based BED format
            bed_start=\$((start - 1))
            
            echo -e "\$chrom\\t\$bed_start\\t\$end" >> chunk_targets.bed
            echo "  Added region: \$chrom:\$start-\$end (BED: \$chrom:\$bed_start-\$end)"
        fi
    done < regions_list.txt
    
    echo "Created BED file for chunk ${chunk_id}:"
    cat chunk_targets.bed
    
    if [[ ! -s chunk_targets.bed ]]; then
        echo "ERROR: BED file is empty for chunk ${chunk_id}"
        exit 1
    fi
    
    # Handle freebayes configuration
    FREEBAYES_OPTS=""
    
    if [[ "${is_config_provided}" == "true" ]] && [[ -f "${config_file}" ]] && [[ "${config_file}" == *.json ]]; then
        echo "Loading freebayes configuration from ${config_file}"
        
        # Parse JSON config
        FREEBAYES_OPTS=\$(python3 << 'EOF'
import json
import sys

try:
    with open('${config_file}', 'r') as f:
        config = json.load(f)
    
    opts = []
    
    # Algorithm parameters
    if 'algorithm_parameters' in config:
        params = config['algorithm_parameters']
        
        for param, flag in [
            ('min_mapping_quality', '--min-mapping-quality'),
            ('min_base_quality', '--min-base-quality'),
            ('min_supporting_allele_qsum', '--min-supporting-allele-qsum'),
            ('min_supporting_mapping_qsum', '--min-supporting-mapping-qsum'),
            ('mismatch_base_quality_threshold', '--mismatch-base-quality-threshold'),
            ('read_mismatch_limit', '--read-mismatch-limit'),
            ('read_max_mismatch_fraction', '--read-max-mismatch-fraction'),
            ('read_snp_limit', '--read-snp-limit'),
            ('read_indel_limit', '--read-indel-limit'),
            ('indel_exclusion_window', '--indel-exclusion-window'),
            ('min_repeat_entropy', '--min-repeat-entropy'),
            ('min_alternate_fraction', '--min-alternate-fraction'),
            ('min_alternate_count', '--min-alternate-count'),
            ('min_alternate_qsum', '--min-alternate-qsum'),
            ('min_alternate_total', '--min-alternate-total'),
            ('min_coverage', '--min-coverage'),
            ('max_coverage', '--max-coverage'),
            ('ploidy', '--ploidy'),
            ('max_complex_gap', '--max-complex-gap'),
            ('haplotype_length', '--haplotype-length'),
            ('min_repeat_length', '--min-repeat-length'),
            ('min_repeat_entropy_for_detection', '--min-repeat-entropy')
        ]:
            if params.get(param) is not None:
                opts.append(f"{flag} {params[param]}")
        
        # Boolean flags
        for param, flag in [
            ('no_partial_observations', '--no-partial-observations'),
            ('pooled_discrete', '--pooled-discrete'),
            ('pooled_continuous', '--pooled-continuous'),
            ('use_duplicate_reads', '--use-duplicate-reads'),
            ('no_mnps', '--no-mnps'),
            ('no_complex', '--no-complex'),
            ('no_snps', '--no-snps'),
            ('no_indels', '--no-indels')
        ]:
            if params.get(param):
                opts.append(flag)
    
    # Population model parameters
    if 'population_model_parameters' in config:
        params = config['population_model_parameters']
        
        for param, flag in [
            ('theta', '--theta'),
            ('posterior_integration_limits', '--posterior-integration-limits'),
            ('gvcf_chunk', '--gvcf-chunk')
        ]:
            if params.get(param) is not None:
                opts.append(f"{flag} {params[param]}")
        
        for param, flag in [
            ('use_reference_allele', '--use-reference-allele'),
            ('gvcf', '--gvcf')
        ]:
            if params.get(param):
                opts.append(flag)
    
    # Genotype likelihood parameters
    if 'genotype_likelihood_parameters' in config:
        params = config['genotype_likelihood_parameters']
        
        for param, flag in [
            ('base_quality_cap', '--base-quality-cap'),
            ('prob_contamination', '--prob-contamination'),
            ('contamination_estimates', '--contamination-estimates')
        ]:
            if params.get(param) is not None:
                opts.append(f"{flag} {params[param]}")
        
        if params.get('legacy_gls'):
            opts.append('--legacy-gls')
    
    # Reporting parameters
    if 'reporting_parameters' in config:
        params = config['reporting_parameters']
        
        for param, flag in [
            ('genotyping_max_iterations', '--genotyping-max-iterations'),
            ('genotyping_max_banddepth', '--genotyping-max-banddepth'),
            ('genotype_variant_threshold', '--genotype-variant-threshold'),
            ('read_dependence_factor', '--read-dependence-factor')
        ]:
            if params.get(param) is not None:
                opts.append(f"{flag} {params[param]}")
        
        for param, flag in [
            ('genotype_qualities', '--genotype-qualities'),
            ('report_genotype_likelihood_max', '--report-genotype-likelihood-max'),
            ('exclude_unobserved_genotypes', '--exclude-unobserved-genotypes'),
            ('use_mapping_quality', '--use-mapping-quality'),
            ('harmonic_indel_quality', '--harmonic-indel-quality')
        ]:
            if params.get(param):
                opts.append(flag)
    
    # Additional options
    if 'additional_options' in config:
        params = config['additional_options']
        
        for param, flag in [
            ('debug', '--debug'),
            ('dd', '--dd')
        ]:
            if params.get(param):
                opts.append(flag)
    
    print(' '.join(opts))

except Exception as e:
    print(f"Error parsing config: {e}", file=sys.stderr)
    sys.exit(1)
EOF
)
        
        echo "Freebayes options from config: \$FREEBAYES_OPTS"
    else
        echo "Using default freebayes parameters"
    fi
    
    # Build BAM arguments
    BAM_ARGS=""
    for bam in *.bam; do
        if [[ -f "\$bam" ]]; then
            BAM_ARGS="\$BAM_ARGS --bam \$bam"
        fi
    done
    
    echo "BAM arguments: \$BAM_ARGS"
    
    # Run freebayes
    echo "Running freebayes on chunk ${chunk_id}..."
    
    freebayes \\
        --fasta-reference ${reference_fa} \\
        --targets chunk_targets.bed \\
        \$BAM_ARGS \\
        \$FREEBAYES_OPTS \\
        --vcf chunk_${chunk_id}.vcf
    
    # Process output
    if [[ -s chunk_${chunk_id}.vcf ]]; then
        echo "Freebayes completed successfully for chunk ${chunk_id}"
        
        variant_count=\$(grep -v '^#' chunk_${chunk_id}.vcf | wc -l)
        echo "Found \$variant_count variants in chunk ${chunk_id}"
        
        echo "Samples in VCF:"
        if grep -q "^#CHROM" chunk_${chunk_id}.vcf; then
            sample_count=\$(grep "^#CHROM" chunk_${chunk_id}.vcf | cut -f10- | tr '\\t' '\\n' | wc -l)
            grep "^#CHROM" chunk_${chunk_id}.vcf | cut -f10- | tr '\\t' '\\n' | head -10
            if [[ \$sample_count -gt 10 ]]; then
                echo "... and \$((sample_count - 10)) more samples"
            fi
            echo "Total samples in VCF: \$sample_count"
        fi
        
        bgzip -c chunk_${chunk_id}.vcf > chunk_${chunk_id}.vcf.gz
        tabix -p vcf chunk_${chunk_id}.vcf.gz
        
        echo "Successfully compressed and indexed chunk ${chunk_id}"
    else
        echo "No variants found in chunk ${chunk_id}, creating empty VCF"
        
        cat > chunk_${chunk_id}.vcf << 'VCF_HEADER'
##fileformat=VCFv4.2
##reference=${reference_fa}
##source=freebayes
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
VCF_HEADER
        
        bgzip -c chunk_${chunk_id}.vcf > chunk_${chunk_id}.vcf.gz
        tabix -p vcf chunk_${chunk_id}.vcf.gz
    fi
    
    # Cleanup
    rm -f chunk_${chunk_id}.vcf regions_list.txt chunk_targets.bed
    
    echo "Chunk ${chunk_id} processing complete"
    """
}