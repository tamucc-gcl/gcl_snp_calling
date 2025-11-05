// modules/filter_bams.nf - MINIMAL TEST VERSION

process FILTER_BAMS {
    tag "${bam.simpleName}"
    
    input:
    path bam
    path bai
    path filter_config
    
    output:
    tuple path("${bam.simpleName}.filtered.bam"), path("${bam.simpleName}.filtered.bam.bai")
    
    script:
    """
    echo "Filtering ${bam}"
    cp ${bam} ${bam.simpleName}.filtered.bam
    cp ${bai} ${bam.simpleName}.filtered.bam.bai
    """
}