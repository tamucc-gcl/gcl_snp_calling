// modules/freebayes_chunk.nf - Clean, simple approach

process FREEBAYES_CHUNK {
    
    input:
    tuple val(chunk_id), val(regions_string)  // separate the tuple components
    path reference          // reference genome FASTA
    path bams               // list of BAM files
    path config             // config file (optional, could be empty)
    
    output:
    tuple val(chunk_id), path("chunk_${chunk_id}.vcf.gz")
    
    script:
    def config_file = config.size() > 0 ? config[0] : null
    def has_config = config_file != null
    """
    echo "Processing chunk ${chunk_id}"
    echo "Regions: ${regions_string}"
    echo "Reference: ${reference}"
    echo "Config provided: ${has_config}"
    
    # List BAM files
    echo "BAM files:"
    for bam in *.bam; do
        if [ -f "\$bam" ]; then
            echo "  \$bam"
        fi
    done
    
    # Count BAM files
    BAM_COUNT=\$(ls -1 *.bam 2>/dev/null | wc -l)
    echo "Total BAM files: \$BAM_COUNT"
    
    if [ \$BAM_COUNT -eq 0 ]; then
        echo "ERROR: No BAM files found!"
        exit 1
    fi
    
    # Ensure BAM indices exist
    for bam in *.bam; do
        if [ -f "\$bam" ] && [ ! -f "\${bam}.bai" ]; then
            echo "Creating index for \$bam"
            samtools index "\$bam"
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
    
    # Build freebayes command
    FREEBAYES_CMD="freebayes --fasta-reference ${reference} --targets chunk.bed"
    
    # Add BAM files
    for bam in *.bam; do
        if [ -f "\$bam" ]; then
            FREEBAYES_CMD="\$FREEBAYES_CMD --bam \$bam"
        fi
    done
    
    # Add config options if provided
    if [ "${has_config}" = "true" ] && [ -f "${config_file}" ]; then
        echo "Loading config from ${config_file}"
        
        # Parse basic config options (simplified)
        MIN_MAPQ=\$(python3 -c "
import json
try:
    with open('${config_file}') as f:
        config = json.load(f)
    print(config.get('algorithm_parameters', {}).get('min_mapping_quality', 20))
except:
    print(20)
" 2>/dev/null)
        
        MIN_BASEQ=\$(python3 -c "
import json
try:
    with open('${config_file}') as f:
        config = json.load(f)
    print(config.get('algorithm_parameters', {}).get('min_base_quality', 20))
except:
    print(20)
" 2>/dev/null)
        
        MIN_ALT_FRAC=\$(python3 -c "
import json
try:
    with open('${config_file}') as f:
        config = json.load(f)
    print(config.get('algorithm_parameters', {}).get('min_alternate_fraction', 0.05))
except:
    print(0.05)
" 2>/dev/null)
        
        MIN_ALT_COUNT=\$(python3 -c "
import json
try:
    with open('${config_file}') as f:
        config = json.load(f)
    print(config.get('algorithm_parameters', {}).get('min_alternate_count', 2))
except:
    print(2)
" 2>/dev/null)
        
        FREEBAYES_CMD="\$FREEBAYES_CMD --min-mapping-quality \$MIN_MAPQ --min-base-quality \$MIN_BASEQ --min-alternate-fraction \$MIN_ALT_FRAC --min-alternate-count \$MIN_ALT_COUNT"
        
        echo "Using config parameters: mapq=\$MIN_MAPQ, baseq=\$MIN_BASEQ, alt_frac=\$MIN_ALT_FRAC, alt_count=\$MIN_ALT_COUNT"
    else
        echo "Using default freebayes parameters"
    fi
    
    # Add output
    FREEBAYES_CMD="\$FREEBAYES_CMD --vcf chunk_${chunk_id}.vcf"
    
    # Run freebayes
    echo "Running: \$FREEBAYES_CMD"
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