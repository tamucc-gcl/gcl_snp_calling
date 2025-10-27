process samtools_stats {
    label 'samtools_stats'
    tag "$sample_id"

    input:
        tuple val(sample_id), path(bam), path(bam_index)

    output:
        tuple val(sample_id), path("${sample_id}.stats"), path("${sample_id}.flagstats")

    script:
    """
    # Generate samtools stats (detailed statistics)
    samtools stats -@ ${task.cpus ?: 4} ${bam} > ${sample_id}.stats
    
    # Generate samtools flagstats (alignment summary)
    samtools flagstats -@ ${task.cpus ?: 4} ${bam} > ${sample_id}.flagstats
    """
}