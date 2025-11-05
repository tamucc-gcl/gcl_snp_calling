// modules/filter_bams.nf - Test without Python

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
    
    # Copy config to fixed name
    if [ "${has_config}" = "true" ]; then
        cp ${filter_config} filter_config.json
    fi
    
    # Set defaults
    MAPPING_MIN_QUALITY=20
    F_FLAGS=264  # 8 + 256 (unpaired mapped + secondary)
    f_FLAGS=2    # proper pairs
    
    # Simple override from config using jq (no Python)
    if [ "${has_config}" = "true" ]; then
        MAPPING_MIN_QUALITY=\$(jq -r '.filter_parameters.mapping_min_quality // 20' filter_config.json)
    fi
    
    # Apply basic samtools filtering
    samtools view -@ ${task.cpus ?: 8} -b -q \$MAPPING_MIN_QUALITY -F \$F_FLAGS -f \$f_FLAGS ${bam} > temp_filtered.bam
    
    # Name-sort, fixmate, coordinate-sort, markdup pipeline
    samtools sort -@ ${task.cpus ?: 8} -n -o ${bam.simpleName}.nsrt.bam temp_filtered.bam
    samtools fixmate -@ ${task.cpus ?: 8} -m -r ${bam.simpleName}.nsrt.bam ${bam.simpleName}.fxmt.bam
    samtools sort -@ ${task.cpus ?: 8} -o ${bam.simpleName}.csrt.bam ${bam.simpleName}.fxmt.bam
    samtools markdup -@ ${task.cpus ?: 8} ${bam.simpleName}.csrt.bam ${bam.simpleName}.bam
    
    rm -f temp_filtered.bam ${bam.simpleName}.nsrt.bam ${bam.simpleName}.fxmt.bam ${bam.simpleName}.csrt.bam
    
    mv ${bam.simpleName}.bam ${bam.simpleName}.filtered.bam
    samtools index -@ ${task.cpus ?: 8} ${bam.simpleName}.filtered.bam
    
    echo "Filtering complete"
    """
}