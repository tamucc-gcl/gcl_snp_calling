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
    reference_ch = Channel.fromPath(params.reference, checkIfExists: true)
    bam_ch = Channel.fromPath(params.bams, checkIfExists: true)
    
    // Config file channel
    if (params.freebayes_config) {
        config_ch = Channel.fromPath(params.freebayes_config, checkIfExists: true)
    } else {
        config_ch = Channel.empty()
    }
    
    // Step 1: Create genome chunks
    chunks = CREATE_CHUNKS(reference_ch, params.num_chunks)
    chunk_regions = chunks[1]  // The regions file
    
    // Step 2: Parse chunk regions into individual chunks
    chunk_ch = chunk_regions
        .splitText()
        .map { line ->
            def parts = line.trim().split('\t')
            if (parts.size() >= 2) {
                return [parts[0], parts[1]]  // [chunk_id, regions_string]
            }
            return null
        }
        .filter { it != null }
    
    // Step 3: Run freebayes on each chunk
    // Use a simpler approach - pass channels directly to process
    vcf_chunks = FREEBAYES_CHUNK(
        chunk_ch,
        reference_ch,
        bam_ch.collect(),
        config_ch.collect().ifEmpty([])
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