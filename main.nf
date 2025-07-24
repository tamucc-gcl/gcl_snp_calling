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

// Subworkflow for processing chunks in parallel
workflow PROCESS_CHUNKS {
    take:
    chunks_ch        // Channel of [chunk_id, regions_string]
    reference        // Reference genome file
    reference_fai    // Reference genome index
    all_bams         // List of all BAM files
    all_bais         // List of all BAI files  
    config_file      // Config file (or empty)
    
    main:
    // Create input tuples for each chunk
    chunk_inputs = chunks_ch
        .combine(reference)
        .combine(reference_fai)  
        .combine(all_bams)
        .combine(all_bais)
        .combine(config_file)
    
    vcf_results = FREEBAYES_CHUNK(chunk_inputs)
    
    emit:
    vcf_results
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
    reference_file = Channel.fromPath(params.reference, checkIfExists: true).first()
    reference_fai = Channel.fromPath("${params.reference}.fai", checkIfExists: true).first()
    
    // Collect all BAM files and their indices
    all_bam_files = Channel.fromPath(params.bams, checkIfExists: true).collect()
    all_bai_files = Channel.fromPath(params.bams, checkIfExists: true)
        .map { bam -> file("${bam}.bai") }
        .collect()
    
    // Handle config file
    config_file = params.freebayes_config ? 
        Channel.fromPath(params.freebayes_config, checkIfExists: true).first() : 
        Channel.value(file("NO_CONFIG"))
    
    // Step 1: Create genome chunks
    chunks = CREATE_CHUNKS(reference_file, params.num_chunks)
    chunk_regions = chunks[1]  // The regions file
    
    // Step 2: Parse chunk regions into individual chunks
    chunks_ch = chunk_regions
        .splitText()
        .map { line ->
            def parts = line.trim().split('\t')
            if (parts.size() >= 2) {
                return [parts[0], parts[1]]  // [chunk_id, regions_string]
            }
            return null
        }
        .filter { it != null }
    
    // Step 3: Process all chunks in parallel using subworkflow
    vcf_chunks = PROCESS_CHUNKS(
        chunks_ch,
        reference_file,
        reference_fai,
        all_bam_files,
        all_bai_files,
        config_file
    )
    
    // Step 4: Combine all VCF files
    all_vcfs = vcf_chunks.map { chunk_id, vcf -> vcf }.collect()
    
    COMBINE_VCFS(all_vcfs, params.output_vcf)
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