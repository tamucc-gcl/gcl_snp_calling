// modules/freebayes_chunk.nf - Updated to handle per-BAM ploidy via CNV map

process FREEBAYES_CHUNK {
    
    input:
    tuple val(chunk_id), val(regions_string)
    path reference
    path reference_fai
    path bams        // List of BAM files
    path bais        // List of BAI files
    path config
    path ploidy_map  // New input for ploidy map file
    
    output:
    tuple val(chunk_id), path("chunk_${chunk_id}.vcf.gz")
    
    script:
    def config_file = config.size() > 0 ? config[0] : null
    def has_config = config_file != null
    def ploidy_file = ploidy_map.name != 'NO_FILE' ? ploidy_map : null
    def has_ploidy_map = ploidy_file != null
    """
    echo "Processing chunk ${chunk_id}"
    echo "Regions: ${regions_string}"
    echo "Reference: ${reference}"
    echo "Reference FAI: ${reference_fai}"
    echo "Config provided: ${has_config}"
    echo "Ploidy map provided: ${has_ploidy_map}"
    
    # List all staged files
    echo "All files in work directory:"
    ls -la
    
    # Count BAM files
    BAM_COUNT=\$(ls -1 *.bam 2>/dev/null | wc -l)
    echo "Total BAM files: \$BAM_COUNT"
    
    if [ \$BAM_COUNT -eq 0 ]; then
        echo "ERROR: No BAM files found!"
        exit 1
    fi
    
    # Verify indices exist
    echo "Verifying BAM indices:"
    for bam in *.bam; do
        if [ -f "\$bam" ]; then
            if [ -f "\${bam}.bai" ]; then
                echo "  \$bam: index present"
            else
                echo "  \$bam: ERROR - index missing!"
                exit 1
            fi
        fi
    done
    
    # Create BED file for regions
    echo "${regions_string}" | tr ',' '\\n' > regions.txt
    
    > chunk.bed
    while read -r region; do
        if [ -n "\$region" ]; then
            chrom=\$(echo "\$region" | cut -d':' -f1)
            pos=\$(echo "\$region" | cut -d':' -f2)
            start=\$(echo "\$pos" | cut -d'-' -f1)
            end=\$(echo "\$pos" | cut -d'-' -f2)
            bed_start=\$((start - 1))
            echo -e "\$chrom\\t\$bed_start\\t\$end" >> chunk.bed
        fi
    done < regions.txt
    
    echo "BED file contents:"
    cat chunk.bed
    
    # Process ploidy map if provided to create CNV map for FreeBayes
    if [ "${has_ploidy_map}" = "true" ]; then
        echo "Creating CNV map from ploidy map..."
        
        # Create a CNV map file in FreeBayes format
        # FreeBayes CNV format: sample_name<TAB>ploidy
        > cnv_map.txt
        
        # Read the ploidy map and create CNV entries
        while read -r line; do
            # Skip comments and empty lines
            [[ "\$line" =~ ^#.*$ ]] && continue
            [[ -z "\$line" ]] && continue
            
            # Parse BAM name and ploidy
            bam_name=\$(echo "\$line" | awk '{print \$1}')
            ploidy=\$(echo "\$line" | awk '{print \$2}')
            
            # Get the sample name from the BAM file
            # FreeBayes uses the SM tag from the @RG header
            if [ -f "\$bam_name" ]; then
                # Extract sample name from BAM header
                sample_name=\$(samtools view -H "\$bam_name" | grep '^@RG' | sed -n 's/.*SM:\\([^\\t]*\\).*/\\1/p' | head -1)
                
                # If no sample name in header, use the BAM filename
                if [ -z "\$sample_name" ]; then
                    sample_name="\${bam_name%.bam}"
                    sample_name="\${sample_name%.filtered}"
                fi
                
                echo -e "\$sample_name\\t\$ploidy" >> cnv_map.txt
                echo "  Sample: \$sample_name, Ploidy: \$ploidy"
            fi
        done < ${ploidy_map}
        
        echo "CNV map created:"
        cat cnv_map.txt
        
        CNV_MAP_OPTION="--cnv-map cnv_map.txt"
    else
        CNV_MAP_OPTION=""
    fi
    
    # Build freebayes command
    FREEBAYES_CMD="freebayes --fasta-reference ${reference} --targets chunk.bed"
    
    # Add CNV map if available
    if [ -n "\$CNV_MAP_OPTION" ]; then
        FREEBAYES_CMD="\$FREEBAYES_CMD \$CNV_MAP_OPTION"
    fi
    
    # Add BAM files
    for bam in *.bam; do
        if [ -f "\$bam" ]; then
            FREEBAYES_CMD="\$FREEBAYES_CMD --bam \$bam"
        fi
    done
    
    # Add config options if provided
    if [ "${has_config}" = "true" ]; then
        CONFIG_FILE=\$(ls *.json 2>/dev/null | head -n1)
        if [ -n "\$CONFIG_FILE" ]; then
            echo "Loading config from \$CONFIG_FILE"
            
            # Parse ALL relevant config options using Python
            python3 << 'PYTHON_PARSE' > freebayes_params.sh
import json
import sys

try:
    with open('${config_file}' if '${config_file}' != '' else next((f for f in ['${config_file}'] + [f for f in __import__('glob').glob('*.json')] if f), None)) as f:
        config = json.load(f)
    
    algo = config.get('algorithm_parameters', {})
    pop = config.get('population_model_parameters', {})
    
    # Algorithm parameters
    params = []
    
    # Basic quality filters
    if 'min_mapping_quality' in algo:
        params.append(f"--min-mapping-quality {algo['min_mapping_quality']}")
    if 'min_base_quality' in algo:
        params.append(f"--min-base-quality {algo['min_base_quality']}")
    if 'min_supporting_allele_qsum' in algo:
        params.append(f"--min-supporting-allele-qsum {algo['min_supporting_allele_qsum']}")
    if 'min_supporting_mapping_qsum' in algo:
        params.append(f"--min-supporting-mapping-qsum {algo['min_supporting_mapping_qsum']}")
    
    # Allele filters
    if 'min_alternate_fraction' in algo:
        params.append(f"--min-alternate-fraction {algo['min_alternate_fraction']}")
    if 'min_alternate_count' in algo:
        params.append(f"--min-alternate-count {algo['min_alternate_count']}")
    if 'min_alternate_qsum' in algo:
        params.append(f"--min-alternate-qsum {algo['min_alternate_qsum']}")
    if 'min_alternate_total' in algo:
        params.append(f"--min-alternate-total {algo['min_alternate_total']}")
    
    # Coverage filters
    if 'min_coverage' in algo:
        params.append(f"--min-coverage {algo['min_coverage']}")
    if algo.get('max_coverage') is not None:
        params.append(f"--max-coverage {algo['max_coverage']}")
    
    # ONLY add global ploidy if no CNV map is provided
    # The CNV map overrides the global ploidy setting
    if not ${has_ploidy_map} and 'ploidy' in algo:
        params.append(f"--ploidy {algo['ploidy']}")
        print(f"echo 'Global ploidy set to: {algo['ploidy']}'")
    
    # Pooling modes still apply globally
    if algo.get('pooled_discrete'):
        params.append("--pooled-discrete")
        print("echo 'Using pooled-discrete mode for pooled sequencing'")
    elif algo.get('pooled_continuous'):
        params.append("--pooled-continuous")
        print("echo 'Using pooled-continuous mode for pooled sequencing'")
    
    # Other algorithm parameters
    if algo.get('use_duplicate_reads'):
        params.append("--use-duplicate-reads")
    if algo.get('no_partial_observations'):
        params.append("--no-partial-observations")
    if algo.get('no_mnps'):
        params.append("--no-mnps")
    if algo.get('no_complex'):
        params.append("--no-complex")
    if algo.get('no_snps'):
        params.append("--no-snps")
    if algo.get('no_indels'):
        params.append("--no-indels")
    
    if 'haplotype_length' in algo:
        params.append(f"--haplotype-length {algo['haplotype_length']}")
    if 'max_complex_gap' in algo:
        params.append(f"--max-complex-gap {algo['max_complex_gap']}")
    if 'min_repeat_entropy' in algo:
        params.append(f"--min-repeat-entropy {algo['min_repeat_entropy']}")
    
    # Population model parameters
    if 'theta' in pop:
        params.append(f"--theta {pop['theta']}")
    if pop.get('use_reference_allele'):
        params.append("--use-reference-allele")
    if pop.get('gvcf'):
        params.append("--gvcf")
    if pop.get('gvcf_chunk') is not None:
        params.append(f"--gvcf-chunk {pop['gvcf_chunk']}")
    
    print(f"FREEBAYES_PARAMS='{' '.join(params)}'")
    
except Exception as e:
    print(f"# Error parsing config: {e}", file=sys.stderr)
    print("FREEBAYES_PARAMS=''")
PYTHON_PARSE
            
            # Source the parameters
            source freebayes_params.sh
            
            # Add parsed parameters to command
            if [ -n "\$FREEBAYES_PARAMS" ]; then
                FREEBAYES_CMD="\$FREEBAYES_CMD \$FREEBAYES_PARAMS"
                echo "Using freebayes parameters: \$FREEBAYES_PARAMS"
            fi
        fi
    else
        echo "Using default freebayes parameters"
    fi
    
    # Add output
    FREEBAYES_CMD="\$FREEBAYES_CMD --vcf chunk_${chunk_id}.vcf"
    
    # Run freebayes
    echo "Full command: \$FREEBAYES_CMD"
    \$FREEBAYES_CMD
    
    # Check output and compress
    if [ -s "chunk_${chunk_id}.vcf" ]; then
        echo "Freebayes completed successfully"
        
        # Count variants
        VARIANT_COUNT=\$(grep -v '^#' chunk_${chunk_id}.vcf | wc -l)
        echo "Found \$VARIANT_COUNT variants"
        
        # Show samples
        if grep -q '^#CHROM' chunk_${chunk_id}.vcf; then
            SAMPLE_COUNT=\$(grep '^#CHROM' chunk_${chunk_id}.vcf | cut -f10- | wc -w)
            echo "Samples (\$SAMPLE_COUNT): \$(grep '^#CHROM' chunk_${chunk_id}.vcf | cut -f10- | tr '\\t' ' ')"
        fi
        
        # Compress and index
        bgzip chunk_${chunk_id}.vcf
        tabix -p vcf chunk_${chunk_id}.vcf.gz
        
    else
        echo "No variants found, creating empty VCF"
        
        # Create minimal header
        echo '##fileformat=VCFv4.2' > chunk_${chunk_id}.vcf
        echo '##reference=${reference}' >> chunk_${chunk_id}.vcf
        echo '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' >> chunk_${chunk_id}.vcf
        
        bgzip chunk_${chunk_id}.vcf
        tabix -p vcf chunk_${chunk_id}.vcf.gz
    fi
    
    echo "Chunk ${chunk_id} complete"
    """
}