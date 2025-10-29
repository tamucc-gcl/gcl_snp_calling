#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { CREATE_CHUNKS } from './modules/create_chunks'
include { FILTER_BAMS } from './modules/filter_bams'
include { FREEBAYES_CHUNK } from './modules/freebayes_chunk'
include { COMBINE_VCFS } from './modules/combine_vcfs'
include { SUMMARIZE_VCFS } from './modules/summarize_vcfs'
include { samtools_stats as SAMTOOLS_STATS_RAW } from './modules/samtools_stats'
include { samtools_stats as SAMTOOLS_STATS_FILTERED } from './modules/samtools_stats'
include { multiqc as MULTIQC_RAW_BAMS } from './modules/multiqc'
include { multiqc as MULTIQC_FILTERED_BAMS } from './modules/multiqc'

// Parameters
params.bams = "*.bam"
params.reference = "reference.fasta"
params.num_chunks = 10
params.output_dir = "results"
params.output_vcf = "raw_variants.vcf.gz"
params.freebayes_config = null
params.bam_filter_config = null
params.ploidy_map = null  // NEW PARAMETER for per-BAM ploidy

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
    --output_vcf        Name of final combined VCF file (default: "raw_variants.vcf.gz")
    --freebayes_config  Path to JSON file containing freebayes parameters (optional)
    --bam_filter_config Path to JSON file containing BAM filter parameters (optional)
    --ploidy_map        Path to file mapping BAM files to ploidy values (optional)
                        Format: bam_filename<TAB>ploidy (e.g., sample1.bam 40)
                        For pooled samples: ploidy = num_individuals × 2 (diploid)
    
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
    log.info "Ploidy map:        ${params.ploidy_map ?: 'Using global ploidy from config or default'}"
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
    
    // Ploidy map channel (NEW)
    if (params.ploidy_map) {
        ploidy_map_ch = Channel.fromPath(params.ploidy_map, checkIfExists: true).first()
    } else {
        ploidy_map_ch = Channel.value(file('NO_FILE'))
    }
    
    // CREATE A SINGLE BAM CHANNEL AND SPLIT IT PROPERLY
    // This ensures consistent ordering and allows parallel execution
    Channel.fromPath(params.bams, checkIfExists: true)
        .toSortedList()  // Ensure consistent ordering
        .flatMap { bams -> bams }
        .map { bam ->
            def bai = file("${bam}.bai")
            if (!bai.exists()) {
                error "BAM index file not found: ${bai}. Please run 'samtools index ${bam}'"
            }
            return tuple(bam.simpleName, bam, bai)
        }
        .into { raw_bam_pairs_ch; bam_pairs_for_filter_ch }  // Split into two channels
    
    // Run samtools stats on raw BAMs (uses first channel)
    raw_bam_stats = SAMTOOLS_STATS_RAW(raw_bam_pairs_ch)
    
    // Run MultiQC on raw BAM stats
    MULTIQC_RAW_BAMS(
        raw_bam_stats
            .map { sid, stats, flagstats -> [stats, flagstats] }
            .flatten()
            .collect(),
        Channel.value('raw_bams')
    )

    // Process BAM files for filtering or direct use (uses second channel)
    if (params.bam_filter_config) {
        // Convert channel for filtering
        bams_for_filter = bam_pairs_for_filter_ch.map { sid, bam, bai -> tuple(bam, bai) }
        
        // Run filtering on each BAM file
        filtered_bams = FILTER_BAMS(
            bams_for_filter.map { it[0] },  // just the BAM
            bams_for_filter.map { it[1] },  // just the BAI
            bam_filter_config_ch
        )
        
        // Run samtools stats on filtered BAMs
        filtered_bam_stats = SAMTOOLS_STATS_FILTERED(
            filtered_bams.map { bam, bai -> 
                tuple(bam.simpleName.replace('.filtered', ''), bam, bai)
            }
        )
        
        // Run MultiQC on filtered BAM stats
        MULTIQC_FILTERED_BAMS(
            filtered_bam_stats
                .map { sid, stats, flagstats -> [stats, flagstats] }
                .flatten()
                .collect(),
            Channel.value('filtered_bams')
        )

        // Collect all filtered BAMs into lists for FreeBayes
        // IMPORTANT: Use collect() to ensure all files are available
        bam_files_ch = filtered_bams.map { bam, bai -> bam }.collect()
        bai_files_ch = filtered_bams.map { bam, bai -> bai }.collect()
        
    } else {
        // Use original BAMs without filtering
        // IMPORTANT: Use collect() to ensure all files are available
        bam_files_ch = bam_pairs_for_filter_ch.map { sid, bam, bai -> bam }.collect()
        bai_files_ch = bam_pairs_for_filter_ch.map { sid, bam, bai -> bai }.collect()
    }
    
    // Config file
    if (params.freebayes_config) {
        config_ch = Channel.fromPath(params.freebayes_config, checkIfExists: true).first()
    } else {
        config_ch = Channel.value([])
    }
    
    // Step 1: Create genome chunks (runs in parallel with BAM processing)
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
    
    // Step 3: For each chunk, run freebayes with ALL BAM files
    // Use combine to ensure each chunk gets all BAM files
    vcf_chunks = FREEBAYES_CHUNK(
        chunk_ch,
        reference_ch,
        reference_fai_ch,
        bam_files_ch,  // This is now a value channel with all BAMs
        bai_files_ch,  // This is now a value channel with all BAIs
        config_ch,
        ploidy_map_ch
    )
    
    // Step 4: Combine VCFs
    all_vcfs = vcf_chunks.map { chunk_id, vcf -> vcf }.collect()
    COMBINE_VCFS(all_vcfs, params.output_vcf)

    SUMMARIZE_VCFS(COMBINE_VCFS.out.vcf)
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
        log.info "QC Reports:"
        log.info "- Raw BAMs MultiQC: ${params.output_dir}/multiqc_reports/multiqc_raw_bams.html"
        if (params.bam_filter_config) {
            log.info "- Filtered BAMs MultiQC: ${params.output_dir}/multiqc_reports/multiqc_filtered_bams.html"
        }
        log.info ""
        log.info "Pipeline used exactly ${params.num_chunks} chunks with complete genome coverage"
        if (params.ploidy_map) {
            log.info "Used per-BAM ploidy values from: ${params.ploidy_map}"
        }
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