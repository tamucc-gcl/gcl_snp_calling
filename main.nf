#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { CREATE_CHUNKS } from './modules/create_chunks'
include { FILTER_BAMS } from './modules/filter_bams'
include { FREEBAYES_CHUNK } from './modules/freebayes_chunk'
include { COMBINE_VCFS } from './modules/combine_vcfs'

// Parameters
params.bams = "*.bam"
params.reference = "reference.fasta"
params.num_chunks = 10
params.output_dir = "results"
params.output_vcf = "combined_variants.vcf.gz"
params.freebayes_config = null
params.bam_filter_config = null

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
    --bam_filter_config Path to JSON file containing BAM filter parameters (optional)
    
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

/*
 * Main Workflow
 */
workflow {
    log.info "================================================================"
    log.info "Parallel Freebayes Variant Calling Pipeline"
    log.info "================================================================"
    log.info "BAM files:         ${params.bams}"
    log.info "Reference:         ${params.reference}"
    log.info "Number of chunks:  ${params.num_chunks}"
    log.info "Output directory:  ${params.output_dir}"
    log.info "Output VCF:        ${params.output_vcf}"
    log.info "Freebayes config:  ${params.freebayes_config ?: 'Using default parameters'}"
    log.info "BAM filter config: ${params.bam_filter_config ?: 'Using default parameters'}"
    log.info "================================================================"
    
    // Create channels - use .first() to make them value channels that can be reused
    reference_ch = Channel.fromPath(params.reference, checkIfExists: true).first()
    reference_fai_ch = Channel.fromPath("${params.reference}.fai", checkIfExists: true).first()
    
    // BAM filter config
    if (params.bam_filter_config) {
        bam_filter_config_ch = Channel.fromPath(params.bam_filter_config, checkIfExists: true).first()
    } else {
        bam_filter_config_ch = Channel.value(file('NO_FILE'))
    }
    
    // Process BAM files - either filter them or use them directly
    if (params.bam_filter_config) {
        // Filter BAMs in parallel
        bam_pairs_ch = Channel.fromPath(params.bams, checkIfExists: true)
            .map { bam ->
                def bai = file("${bam}.bai")
                if (!bai.exists()) {
                    error "BAM index file not found: ${bai}. Please run 'samtools index ${bam}'"
                }
                return tuple(bam, bai)
            }
        
        // Run filtering on each BAM file
        filtered_bams = FILTER_BAMS(
            bam_pairs_ch.map { it[0] },  // just the BAM
            bam_pairs_ch.map { it[1] },  // just the BAI
            bam_filter_config_ch
        )
        
        // Collect filtered BAMs and BAIs
        bam_bai_split = filtered_bams
            .toList()
            .multiMap { pairs ->
                bams: pairs.collect { it[0] }
                bais: pairs.collect { it[1] }
            }
        
        bam_files = bam_bai_split.bams
        bai_files = bam_bai_split.bais
    } else {
        // Use original BAMs without filtering
        bam_bai_split = Channel.fromPath(params.bams, checkIfExists: true)
            .map { bam ->
                def bai = file("${bam}.bai")
                if (!bai.exists()) {
                    error "BAM index file not found: ${bai}. Please run 'samtools index ${bam}'"
                }
                return tuple(bam, bai)
            }
            .toList()
            .multiMap { pairs ->
                bams: pairs.collect { it[0] }
                bais: pairs.collect { it[1] }
            }
        
        bam_files = bam_bai_split.bams
        bai_files = bam_bai_split.bais
    }
    
    // Config file
    if (params.freebayes_config) {
        config_ch = Channel.fromPath(params.freebayes_config, checkIfExists: true).first()
    } else {
        config_ch = Channel.value([])
    }
    
    // Step 1: Create genome chunks (runs in parallel with BAM filtering)
    chunks = CREATE_CHUNKS(reference_ch, params.num_chunks)
    chunk_regions = chunks[1]
    
    // Step 2: Parse chunks into individual emissions
    chunk_ch = chunk_regions
        .splitText()
        .map { line ->
            def parts = line.trim().split('\t')
            if (parts.size() >= 2) {
                return tuple(parts[0], parts[1])  // [chunk_id, regions_string]
            }
            return null
        }
        .filter { it != null }
    
    // Step 3: For each chunk, run freebayes with ALL inputs
    // This will wait for both BAM filtering and chunking to complete
    vcf_chunks = FREEBAYES_CHUNK(
        chunk_ch,
        reference_ch,
        reference_fai_ch,
        bam_files,
        bai_files,
        config_ch
    )
    
    // Step 4: Combine VCFs
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