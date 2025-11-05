#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { CREATE_CHUNKS } from './modules/create_chunks'
include { FILTER_BAMS } from './modules/filter_bams'
include { FREEBAYES_CHUNK } from './modules/freebayes_chunk'
include { COMBINE_VCFS } from './modules/combine_vcfs'
include { SUMMARIZE_VCFS } from './modules/summarize_vcfs'

// QC modules - import but make optional
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
params.run_qc = false  // QC disabled by default to not interfere with caching

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
    --run_qc            Run QC stats and MultiQC (default: false - enable only after main pipeline works)
    
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

// Main variant calling workflow - keep this clean and simple
workflow variant_calling {
    
    // Create value channels for reference files
    reference_ch = Channel.fromPath(params.reference, checkIfExists: true).first()
    reference_fai_ch = Channel.fromPath("${params.reference}.fai", checkIfExists: true).first()
    
    // Config channels
    bam_filter_config_ch = params.bam_filter_config 
        ? Channel.fromPath(params.bam_filter_config, checkIfExists: true).first()
        : Channel.value(file('NO_FILE'))
    
    ploidy_map_ch = params.ploidy_map 
        ? Channel.fromPath(params.ploidy_map, checkIfExists: true).first()
        : Channel.value(file('NO_FILE'))
    
    config_ch = params.freebayes_config 
        ? Channel.fromPath(params.freebayes_config, checkIfExists: true).first()
        : Channel.value([])
    
    // Step 1: Create genome chunks
    chunks = CREATE_CHUNKS(reference_ch, params.num_chunks)
    chunk_regions = chunks[1]
    
    // Step 2: Parse chunks into individual emissions
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
    
    // Step 3: Process BAM files for filtering (if config provided)
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
        
        // Collect filtered BAMs for FreeBayes
        bam_files_ch = filtered_bams.map { bam, bai -> bam }.collect()
        bai_files_ch = filtered_bams.map { bam, bai -> bai }.collect()
        
    } else {
        // Use original BAMs without filtering
        bam_files_ch = Channel.fromPath(params.bams, checkIfExists: true).collect()
        bai_files_ch = Channel.fromPath(params.bams, checkIfExists: true)
            .map { bam -> file("${bam}.bai") }
            .collect()
    }
    
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
    SUMMARIZE_VCFS(COMBINE_VCFS.out.vcf, ploidy_map_ch)
}

// Separate QC workflow - only run if explicitly requested
workflow qc_workflow {
    
    // QC on raw BAMs
    raw_bam_channel = Channel.fromPath(params.bams, checkIfExists: true)
        .map { bam ->
            def bai = file("${bam}.bai")
            if (!bai.exists()) {
                error "BAM index file not found: ${bai}. Please run 'samtools index ${bam}'"
            }
            return tuple(bam.simpleName, bam, bai)
        }
    
    raw_bam_stats = SAMTOOLS_STATS_RAW(raw_bam_channel)
    
    raw_stats_files = raw_bam_stats
        .map { sid, stats, flagstats -> [stats, flagstats] }
        .flatten()
        .collect()
    
    MULTIQC_RAW_BAMS(
        raw_stats_files,
        Channel.value('raw_bams')
    )
    
    // QC on filtered BAMs if they exist
    if (params.bam_filter_config) {
        // Look for filtered BAMs in the output directory
        filtered_bam_channel = Channel.fromPath("${params.output_dir}/filtered_bams/*.filtered.bam", checkIfExists: false)
            .map { bam ->
                def bai = file("${bam}.bai")
                if (bai.exists()) {
                    return tuple(bam.simpleName.replace('.filtered', ''), bam, bai)
                }
                return null
            }
            .filter { it != null }
            .ifEmpty { 
                log.warn "No filtered BAMs found for QC. Run main pipeline first."
                return Channel.empty()
            }
        
        if (filtered_bam_channel) {
            filtered_bam_stats = SAMTOOLS_STATS_FILTERED(filtered_bam_channel)
            
            filtered_stats_files = filtered_bam_stats
                .map { sid, stats, flagstats -> [stats, flagstats] }
                .flatten()
                .collect()
            
            MULTIQC_FILTERED_BAMS(
                filtered_stats_files,
                Channel.value('filtered_bams')
            )
        }
    }
}

// Main workflow
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
    log.info "Run QC:            ${params.run_qc}"
    log.info "================================================================"
    
    // Always run the main variant calling pipeline
    variant_calling()
    
    // Only run QC if explicitly requested (to avoid cache interference)
    if (params.run_qc) {
        log.info "Running QC workflow..."
        qc_workflow()
    } else {
        log.info "QC workflow skipped. Use --run_qc to enable QC stats and MultiQC reports."
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
        
        if (params.run_qc) {
            log.info "QC Reports:"
            log.info "- Raw BAMs MultiQC: ${params.output_dir}/multiqc_reports/multiqc_raw_bams.html"
            if (params.bam_filter_config) {
                log.info "- Filtered BAMs MultiQC: ${params.output_dir}/multiqc_reports/multiqc_filtered_bams.html"
            }
        } else {
            log.info "To generate QC reports, run with --run_qc"
        }
        
        log.info ""
        log.info "Pipeline used exactly ${params.num_chunks} chunks with complete genome coverage"
        if (params.ploidy_map) {
            log.info "Used per-BAM ploidy values from: ${params.ploidy_map}"
        }
        log.info ""
        log.info "You can view the execution reports at:"
        log.info "- Timeline: ${params.output_dir}/pipeline/snp_pipeline_timeline.html"
        log.info "- Report:   ${params.output_dir}/pipeline/snp_pipeline_report.html"
        log.info "- Trace:    ${params.output_dir}/pipeline/snp_pipeline_trace.txt"
    } else {
        log.info ""
        log.info "Pipeline failed. Check the error messages above."
    }
    log.info "================================================================"
}