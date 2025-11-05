#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.bams = "./data/bam/*.bam"

process TEST_PROCESS {
    tag "$sample_id"
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    
    output:
    path("${sample_id}.txt")
    
    script:
    """
    echo "Processing ${sample_id}" > ${sample_id}.txt
    samtools flagstat ${bam} >> ${sample_id}.txt
    """
}

workflow {
    // Simplest possible channel
    bam_ch = Channel.fromPath(params.bams)
        .map { bam ->
            tuple(bam.simpleName, bam, file("${bam}.bai"))
        }
    
    TEST_PROCESS(bam_ch)
}