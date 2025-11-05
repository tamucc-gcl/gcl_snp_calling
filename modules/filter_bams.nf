// modules/filter_bams.nf - Fixed with deterministic config filename

process FILTER_BAMS {
    tag "${bam.simpleName}"
    
    input:
    path bam
    path bai
    path filter_config
    
    output:
    tuple path("${bam.simpleName}.filtered.bam"), path("${bam.simpleName}.filtered.bam.bai")
    
    script:
    def has_config = filter_config.name != 'NO_FILE'
    """
    echo "Filtering BAM file: ${bam}"
    
    # Copy config to fixed name for deterministic processing
    if [ "${has_config}" = "true" ]; then
        cp ${filter_config} filter_config.json
        CONFIG_FILE="filter_config.json"
    else
        CONFIG_FILE=""
    fi
    
    # Initialize filtering parameters
    MAPPING_MIN_QUALITY=20
    REMOVE_ORPHANS="yes"
    
    # Load parameters from config if provided
    if [ -n "\$CONFIG_FILE" ]; then
        python3 << 'PYTHON_SCRIPT' > filter_params.sh
import json

with open('filter_config.json', 'r') as f:
    config = json.load(f)

filters = config.get('filter_parameters', {})
print(f"MAPPING_MIN_QUALITY={filters.get('mapping_min_quality', 20)}")
print(f"REMOVE_ORPHANS='{filters.get('remove_reads_orphaned_by_filters', 'yes')}'")
PYTHON_SCRIPT
        
        source filter_params.sh
    fi
    
    # Simple filter for testing
    samtools view -@ 8 -b -q \$MAPPING_MIN_QUALITY ${bam} > temp.bam
    
    if [ "\$REMOVE_ORPHANS" = "yes" ]; then
        samtools sort -@ 8 -n temp.bam -o temp_ns.bam
        samtools fixmate -@ 8 -r temp_ns.bam temp_fix.bam
        samtools sort -@ 8 temp_fix.bam -o ${bam.simpleName}.filtered.bam
        rm temp_ns.bam temp_fix.bam temp.bam
    else
        mv temp.bam ${bam.simpleName}.filtered.bam
    fi
    
    samtools index -@ 8 ${bam.simpleName}.filtered.bam
    
    echo "Filtering complete"
    """
}