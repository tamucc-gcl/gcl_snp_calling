#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.reference = "./genome/genome.fa"

process SIMPLE_TEST {
    tag "simple"
    
    input:
    path reference
    
    output:
    path "output.txt"
    
    script:
    """
    echo "Testing simple caching" > output.txt
    samtools faidx ${reference} > output.txt
    """
}

workflow {
    ref = Channel.fromPath(params.reference, checkIfExists: true).first()
    SIMPLE_TEST(ref)
}