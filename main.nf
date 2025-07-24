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
    
    // Create input channels separately and explicitly
    reference_file = Channel.fromPath(params.reference, checkIfExists: true)
    
    // Get reference FAI file (should be alongside the reference)
    reference_fai = Channel.fromPath("${params.reference}.fai", checkIfExists: true)
    
    // Collect BAM files
    bam_files_list = Channel.fromPath(params.bams, checkIfExists: true)
        .collect()
        .map { bams ->
            log.info "DEBUG: Collected ${bams.size()} BAM files: ${bams.collect{it.name}.join(', ')}"
            return bams
        }
    
    // Collect corresponding BAM index files
    bam_index_files = Channel.fromPath(params.bams, checkIfExists: true)
        .map { bam ->
            def bai_file = file("${bam}.bai")
            if (bai_file.exists()) {
                return bai_file
            } else {
                def alt_bai = file("${bam.toString().replaceAll('\\.bam$', '.bai')}")
                if (alt_bai.exists()) {
                    return alt_bai
                } else {
                    // Return placeholder - will be created if needed
                    return file("${bam}.bai")
                }
            }
        }
        .collect()
    
    // Handle config file
    if (params.freebayes_config && params.freebayes_config != "null" && params.freebayes_config != null) {
        config_file_ch = Channel.fromPath(params.freebayes_config, checkIfExists: true)
        log.info "DEBUG: Using config file: ${params.freebayes_config}"
    } else {
        // Create a dummy file for "no config"
        config_file_ch = Channel.value(file("NO_CONFIG_FILE"))
        log.info "DEBUG: No config file provided"
    }
    
    // Step 1: Create genome chunks
    chunk_results = CREATE_CHUNKS(reference_file, params.num_chunks)
    chunks_bed_file = chunk_results[0]
    chunk_regions_file = chunk_results[1]
    
    // Step 2: Parse chunk regions
    chunks_channel = chunk_regions_file
        .splitText()
        .map { line ->
            def trimmed = line.trim()
            if (trimmed && !trimmed.startsWith('#')) {
                def parts = trimmed.split('\t')
                if (parts.size() >= 2) {
                    return [parts[0], parts[1]]  // [chunk_id, regions_string]
                }
            }
            return null
        }
        .filter { it != null }
    
    // Step 3: Create inputs for FREEBAYES_CHUNK
    // Process signature: tuple val(chunk_id), val(regions_string), path(reference_fa), path(reference_fai), path(bam_files), path(bai_files), path(config_file)
    
    freebayes_inputs = chunks_channel
        .combine(reference_file)
        .combine(reference_fai)
        .combine(bam_files_list)
        .combine(bam_index_files)
        .combine(config_file_ch)
        .map { it ->
            // it is [chunk_info, ref_fa, ref_fai, bams, bais, config]
            // where chunk_info is [chunk_id, regions_string]
            def chunk_info = it[0]
            def ref_fa = it[1]
            def ref_fai = it[2]
            def bams = it[3]
            def bais = it[4]
            def config = it[5]
            return [chunk_info[0], chunk_info[1], ref_fa, ref_fai, bams, bais, config]
        }
    
    // Run freebayes
    vcf_results = FREEBAYES_CHUNK(freebayes_inputs)
    
    // Step 4: Combine VCFs
    all_vcf_files = vcf_results
        .map { chunk_id, vcf_file -> vcf_file }
        .collect()
    
    COMBINE_VCFS(all_vcf_files, params.output_vcf)
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