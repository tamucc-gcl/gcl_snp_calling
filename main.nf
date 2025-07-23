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
    --num_chunks    Number of chunks to split genome into (default: 10)
    --output_dir    Output directory (default: "results")
    --output_vcf    Name of final combined VCF file (default: "combined_variants.vcf.gz")
    
    Example:
    nextflow run pipeline.nf \\
        --bams "data/*.bam" \\
        --reference "genome.fasta" \\
        --num_chunks 20 \\
        --output_dir "variant_calls"
    
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
    log.info "BAM files:        ${params.bams}"
    log.info "Reference:        ${params.reference}"
    log.info "Number of chunks: ${params.num_chunks}"
    log.info "Output directory: ${params.output_dir}"
    log.info "Output VCF:       ${params.output_vcf}"
    log.info "================================================================"
    
    // Create input channels
    bam_ch = Channel.fromPath(params.bams, checkIfExists: true)
    reference_ch = Channel.fromPath(params.reference, checkIfExists: true)
    
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
    chunks_bed = CREATE_CHUNKS(reference_ch, params.num_chunks)
    
    // Convert BED file to channel of regions
    chunks_ch = chunks_bed
        .splitText()
        .map { line ->
            def parts = line.trim().split('\t')
            def chrom = parts[0]
            def start = parts[1].toInteger() + 1  // Convert from 0-based to 1-based
            def end = parts[2].toInteger()
            def region = "${chrom}:${start}-${end}"
            def chunk_id = "${chrom}_${parts[1]}_${parts[2]}"
            return [chunk_id, region]
        }
    
    // Step 2: Run freebayes on each chunk
    vcf_chunks = FREEBAYES_CHUNK(
        reference_ch,
        bam_files,
        bam_indices_ch,
        chunks_ch
    )
    
    // Step 3: Combine all VCF files
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
        log.info "You can view the execution reports at:"
        log.info "- Timeline: ${params.output_dir}/timeline.html"
        log.info "- Report:   ${params.output_dir}/report.html"
        log.info "- Trace:    ${params.output_dir}/trace.txt"
    } else {
        log.info ""
        log.info "Pipeline failed. Check the error messages above."
    }
    log.info "================================================================"
}