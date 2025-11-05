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
params.ploidy_map = null
params.skip_qc = false  // New parameter to skip QC if needed

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
    --skip_qc           Skip QC stats and MultiQC (default: false)
    
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
    log.info "Skip QC:           ${params.skip_qc}"
    log.info "================================================================"
    
    // Create value channels for reference files
    reference_ch = Channel.fromPath(params.reference, checkIfExists: true).first()
    reference_fai_ch = Channel.fromPath("${params.reference}.fai", checkIfExists: true).first()
    
    // Config channels
    if (params.bam_filter_config) {
        bam_filter_config_ch = Channel.fromPath(params.bam_filter_config, checkIfExists: true).first()
    } else {
        bam_filter_config_ch = Channel.value(file('NO_FILE'))
    }
    
    if (params.ploidy_map) {
        ploidy_map_ch = Channel.fromPath(params.ploidy_map, checkIfExists: true).first()
    } else {
        ploidy_map_ch = Channel.value(file('NO_FILE'))
    }
    
    if (params.freebayes_config) {
        config_ch = Channel.fromPath(params.freebayes_config, checkIfExists: true).first()
    } else {
        config_ch = Channel.value([])
    }
    
    // ===== MAIN VARIANT CALLING WORKFLOW =====
    // This is kept separate from QC to maintain caching
    
    // Step 1: Create genome chunks (independent of BAMs)
    chunks = CREATE_CHUNKS(reference_ch, params.num_chunks)
    chunk_regions = chunks[1]
    
    // Step 2: Process BAM files for variant calling
    if (params.bam_filter_config) {
        // Create channel for filtering
        bam_channel = Channel.fromPath(params.bams, checkIfExists: true)
            .map { bam ->
                def bai = file("${bam}.bai")
                if (!bai.exists()) {
                    error "BAM index file not found: ${bai}. Please run 'samtools index ${bam}'"
                }
                return tuple(bam, bai)
            }
        
        // Filter BAMs
        filtered_bams = FILTER_BAMS(
            bam_channel.map { it[0] },
            bam_channel.map { it[1] },
            bam_filter_config_ch
        )
        
        // Collect for FreeBayes
        bam_files_ch = filtered_bams.map { bam, bai -> bam }.collect()
        bai_files_ch = filtered_bams.map { bam, bai -> bai }.collect()
        
    } else {
        // Use original BAMs without filtering
        bam_channel = Channel.fromPath(params.bams, checkIfExists: true)
            .map { bam ->
                def bai = file("${bam}.bai")
                if (!bai.exists()) {
                    error "BAM index file not found: ${bai}. Please run 'samtools index ${bam}'"
                }
                return tuple(bam, bai)
            }
        
        bam_files_ch = bam_channel.map { bam, bai -> bam }.collect()
        bai_files_ch = bam_channel.map { bam, bai -> bai }.collect()
    }
    
    // Step 3: Parse chunks into individual emissions
    chunk_ch = chunk_regions
        .splitText()
        .map { line ->
            def parts = line.trim().split('\t')
            if (parts.size() >= 2) {
                return tuple(parts[0], parts[1])
            }
            return null
        }
        .filter { it != null }
    
    // Step 4: Run freebayes on each chunk
    vcf_chunks = FREEBAYES_CHUNK(
        chunk_ch,
        reference_ch,
        reference_fai_ch,
        bam_files_ch,
        bai_files_ch,
        config_ch,
        ploidy_map_ch
    )
    
    // Step 5: Combine VCFs
    all_vcfs = vcf_chunks.map { chunk_id, vcf -> vcf }.collect()
    COMBINE_VCFS(all_vcfs, params.output_vcf)
    
    // Step 6: Summarize final VCF
    SUMMARIZE_VCFS(COMBINE_VCFS.out.vcf)
    
    // ===== QC WORKFLOW (INDEPENDENT) =====
    // Run QC in parallel, not affecting the main pipeline
    if (!params.skip_qc) {
        
        // Create separate channels for QC (independent of main workflow)
        qc_bam_channel = Channel.fromPath(params.bams, checkIfExists: true)
            .map { bam ->
                def bai = file("${bam}.bai")
                if (!bai.exists()) {
                    error "BAM index file not found: ${bai}. Please run 'samtools index ${bam}'"
                }
                return tuple(bam.simpleName, bam, bai)
            }
        
        // Run samtools stats on raw BAMs
        raw_bam_stats = SAMTOOLS_STATS_RAW(qc_bam_channel)
        
        // Collect stats for MultiQC
        raw_stats_files = raw_bam_stats
            .map { sid, stats, flagstats -> [stats, flagstats] }
            .flatten()
            .collect()
        
        // Run MultiQC on raw BAM stats
        MULTIQC_RAW_BAMS(
            raw_stats_files,
            Channel.value('raw_bams')
        )
        
        // If filtering was done, also run QC on filtered BAMs
        if (params.bam_filter_config) {
            // Create another independent channel for filtered BAM QC
            qc_filter_channel = Channel.fromPath(params.bams, checkIfExists: true)
                .map { bam ->
                    def bai = file("${bam}.bai")
                    return tuple(bam, bai)
                }
            
            // Filter BAMs independently for QC
            filtered_bams_qc = FILTER_BAMS(
                qc_filter_channel.map { it[0] },
                qc_filter_channel.map { it[1] },
                bam_filter_config_ch
            )
            
            // Run samtools stats on filtered BAMs
            filtered_bam_stats = SAMTOOLS_STATS_FILTERED(
                filtered_bams_qc.map { bam, bai -> 
                    tuple(bam.simpleName.replace('.filtered', ''), bam, bai)
                }
            )
            
            // Collect filtered stats for MultiQC
            filtered_stats_files = filtered_bam_stats
                .map { sid, stats, flagstats -> [stats, flagstats] }
                .flatten()
                .collect()
            
            // Run MultiQC on filtered BAM stats
            MULTIQC_FILTERED_BAMS(
                filtered_stats_files,
                Channel.value('filtered_bams')
            )
        }
    }
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
        if (!params.skip_qc) {
            log.info "QC Reports:"
            log.info "- Raw BAMs MultiQC: ${params.output_dir}/multiqc_reports/multiqc_raw_bams.html"
            if (params.bam_filter_config) {
                log.info "- Filtered BAMs MultiQC: ${params.output_dir}/multiqc_reports/multiqc_filtered_bams.html"
            }
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