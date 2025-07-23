#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { CREATE_CHUNKS } from './modules/create_chunks'
include { FREEBAYES_CHUNK } from './modules/freebayes_chunk'
include { COMBINE_VCFS } from './modules/combine_vcfs'

// Parameters
params.bams = "*.bam"
params.reference = "reference.fasta"
params.num_chunks = 10  // Number of chunks to split genome into
params.output_dir = "results"
params.output_vcf = "combined_variants.vcf.gz"
params.freebayes_config = null  // Optional: path to freebayes config JSON file

// Help message
def helpMessage() {
    log.info"""
    ================================================================
    Parallel Freebayes Variant Calling Pipeline
    ================================================================
    
    Usage:
    nextflow run pipeline.nf --bams "path/to/*.bam" --reference "reference.fasta"
    
    Required parameters:
    --bams          Path to BAM files (glob pattern, e.g., "*.bam")
    --reference     Path to reference genome FASTA file
    
    Optional parameters:
    --num_chunks        Number of chunks to split genome into (default: 10)
    --output_dir        Output directory (default: "results")
    --output_vcf        Name of final combined VCF file (default: "combined_variants.vcf.gz")
    --freebayes_config  Path to JSON file containing freebayes parameters (optional)
    
    Example:
    nextflow run pipeline.nf \\
        --bams "data/*.bam" \\
        --reference "genome.fasta" \\
        --num_chunks 50 \\
        --output_dir "variant_calls" \\
        --freebayes_config "freebayes_params.json"
    
    Note: This version creates exactly the specified number of chunks,
    with each chunk potentially spanning multiple contigs for complete
    genome coverage.
    
    """.stripIndent()
}

// Show help message if requested
if (params.help) {
    helpMessage()
    exit 0
}

// Input validation
if (!params.bams || !params.reference) {
    log.error "ERROR: Please provide both --bams and --reference parameters"
    helpMessage()
    exit 1
}

if (params.num_chunks <= 0) {
    log.error "ERROR: --num_chunks must be a positive integer"
    exit 1
}

/*
 * Main Workflow
 */
workflow {
    // Print pipeline info
    log.info "================================================================"
    log.info "Parallel Freebayes Variant Calling Pipeline"
    log.info "================================================================"
    log.info "BAM files:         ${params.bams}"
    log.info "Reference:         ${params.reference}"
    log.info "Number of chunks:  ${params.num_chunks}"
    log.info "Output directory:  ${params.output_dir}"
    log.info "Output VCF:        ${params.output_vcf}"
    log.info "Freebayes config:  ${params.freebayes_config ?: 'Using default parameters'}"
    log.info "================================================================"
    
    // Create input channels
    bam_ch = Channel.fromPath(params.bams, checkIfExists: true)
    reference_ch = Channel.fromPath(params.reference, checkIfExists: true)
    
    // Create config channel - use empty file if no config provided
    config_ch = params.freebayes_config ? 
        Channel.fromPath(params.freebayes_config, checkIfExists: true) : 
        Channel.value(file("NO_CONFIG"))
    
    // Collect BAM files and check for indices
    bam_files = bam_ch.collect()
    
    // Create channel for BAM indices (they'll be created if missing)
    bam_indices_ch = bam_ch.map { bam ->
        def bai = file("${bam}.bai")
        def bai_alt = file("${bam.baseName}.bai")
        if (bai.exists()) {
            return bai
        } else if (bai_alt.exists()) {
            return bai_alt
        } else {
            // Will be created in the freebayes process
            return file("${bam}.bai") // placeholder that will be created
        }
    }.collect()
    
    // Step 1: Create genome chunks
    chunk_outputs = CREATE_CHUNKS(reference_ch, params.num_chunks)
    chunks_bed = chunk_outputs[0]
    chunk_regions_file = chunk_outputs[1]
    
    // Step 2: Parse chunk regions file to create input channel for freebayes
    chunks_ch = chunk_regions_file
        .splitText()
        .map { line ->
            def trimmed = line.trim()
            if (trimmed) {  // Only process non-empty lines
                def parts = trimmed.split('\t')
                if (parts.size() >= 2) {
                    def chunk_id = parts[0]
                    def regions_string = parts[1]
                    return [chunk_id, regions_string]
                }
            }
            return null
        }
        .filter { it != null }  // Remove null entries
    
    // Step 3: Run freebayes on each chunk in parallel
    // Create tuples matching the new FREEBAYES_CHUNK input signature:
    // tuple val(chunk_id), val(regions_string), path(reference), path(bams), path(bam_indices), path(config)
    
    freebayes_inputs = chunks_ch
        .combine(reference_ch)
        .combine(bam_files)
        .combine(bam_indices_ch)
        .combine(config_ch)
        .map { it ->
            // 'it' is a LinkedList containing all combined elements
            // Destructure the list elements
            def chunk_id = it[0]
            def regions_string = it[1]
            def ref = it[2]
            def bams = it[3]
            def indices = it[4]
            def config = it[5]
            
            // Return tuple matching the new process input signature
            return tuple(chunk_id, regions_string, ref, bams, indices, config)
        }
    
    vcf_chunks = FREEBAYES_CHUNK(freebayes_inputs)
    
    // Step 4: Combine all VCF files
    all_vcfs = vcf_chunks.map { chunk_id, vcf -> vcf }.collect()
    
    COMBINE_VCFS(
        all_vcfs,
        params.output_vcf
    )
}

workflow.onComplete {
    log.info "================================================================"
    log.info "Pipeline Summary"
    log.info "================================================================"
    log.info "Completed at:     ${workflow.complete}"
    log.info "Duration:         ${workflow.duration}"
    log.info "Success:          ${workflow.success}"
    log.info "Exit status:      ${workflow.exitStatus}"
    log.info "Error report:     ${workflow.errorReport ?: 'None'}"
    
    if (workflow.success) {
        log.info ""
        log.info "Results are available in: ${params.output_dir}/"
        log.info "Final VCF file: ${params.output_dir}/${params.output_vcf}"
        log.info ""
        log.info "Pipeline used exactly ${params.num_chunks} chunks with complete genome coverage"
        log.info ""
        log.info "You can view the execution reports at:"
        log.info "- Timeline: ${params.output_dir}/pipeline/pipeline_timeline.html"
        log.info "- Report:   ${params.output_dir}/pipeline/pipeline_report.html"
        log.info "- Trace:    ${params.output_dir}/pipeline/pipeline_trace.txt"
    } else {
        log.info ""
        log.info "Pipeline failed. Check the error messages above."
    }
    log.info "================================================================"
}