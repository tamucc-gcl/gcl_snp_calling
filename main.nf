#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import original modules
include { CREATE_CHUNKS } from './modules/create_chunks'
include { FILTER_BAMS } from './modules/filter_bams'
include { FREEBAYES_CHUNK } from './modules/freebayes_chunk'
include { COMBINE_VCFS } from './modules/combine_vcfs'
include { SUMMARIZE_VCFS } from './modules/summarize_vcfs'

// Import comprehensive QC modules
include { VARIANT_STATS } from './modules/variant_stats'
include { SAMPLE_QC } from './modules/sample_qc'
include { LOCUS_QC } from './modules/locus_qc'
include { CONTAMINATION_CHECK } from './modules/contamination_check'
include { GENERATE_PLOTS } from './modules/generate_plots'
include { FILTER_RECOMMENDATIONS } from './modules/filter_recommendations'
include { FILTER_VALIDATION } from './modules/filter_validation'
include { PUBLICATION_TABLES } from './modules/publication_tables'
include { GENERATE_REPORT } from './modules/generate_report'

// Parameters
params.bams = "*.bam"
params.reference = "reference.fasta"
params.num_chunks = 10
params.output_dir = "results"
params.output_vcf = "combined_variants.vcf.gz"
params.freebayes_config = null
params.bam_filter_config = null
params.ploidy_map = null
params.generate_report = true
params.check_contamination = true
params.validate_filters = true
params.generate_publication_tables = false

// Help message
def helpMessage() {
    log.info"""
    ================================================================
    Parallel Freebayes Variant Calling Pipeline with Comprehensive QC
    ================================================================
    
    Usage:
    nextflow run main.nf --bams "path/to/*.bam" --reference "reference.fasta"
    
    Required parameters:
    --bams                      Path to BAM files (glob pattern)
    --reference                 Path to reference genome FASTA file
    
    Optional parameters:
    --num_chunks                Number of chunks to split genome (default: 10)
    --output_dir                Output directory (default: "results")
    --output_vcf                Name of final VCF file (default: "combined_variants.vcf.gz")
    --freebayes_config          FreeBayes parameters JSON file
    --bam_filter_config         BAM filter parameters JSON file
    --ploidy_map                Ploidy map for pooled samples (TSV format)
    
    QC and Reporting options:
    --generate_report           Generate comprehensive QC report (default: true)
    --check_contamination       Check for sample contamination (default: true)
    --validate_filters          Validate filter effectiveness (default: true)
    --generate_publication_tables  Generate publication-ready tables (default: false)
    
    Examples:
    # Standard diploid samples
    nextflow run main.nf \\
        --bams "data/*.bam" \\
        --reference "ref.fa" \\
        --freebayes_config config/diploid.json
    
    # Pooled samples with ploidy map
    nextflow run main.nf \\
        --bams "pools/*.bam" \\
        --reference "ref.fa" \\
        --ploidy_map pools.txt \\
        --freebayes_config config/pooled.json
    
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
 * Main Workflow with Comprehensive QC Pipeline
 */
workflow {
    log.info """
    ╔════════════════════════════════════════════════════════════════╗
    ║  Parallel FreeBayes Variant Calling Pipeline with QC Suite    ║
    ╠════════════════════════════════════════════════════════════════╣
    ║ Input Configuration:                                           ║
    ║   BAM files:         ${params.bams}                           
    ║   Reference:         ${params.reference}                      
    ║   Chunks:            ${params.num_chunks}                     
    ║   Output:            ${params.output_dir}/${params.output_vcf}
    ║                                                                ║
    ║ Analysis Options:                                              ║
    ║   FreeBayes config:  ${params.freebayes_config ?: 'Default'}  
    ║   BAM filters:       ${params.bam_filter_config ?: 'None'}    
    ║   Ploidy map:        ${params.ploidy_map ?: 'Diploid'}        
    ║                                                                ║
    ║ QC & Reporting:                                                ║
    ║   Generate report:   ${params.generate_report}                
    ║   Contamination:     ${params.check_contamination}            
    ║   Validate filters:  ${params.validate_filters}               
    ║   Publication tables: ${params.generate_publication_tables}    
    ╚════════════════════════════════════════════════════════════════╝
    """.stripIndent()
    
    // ========== SETUP CHANNELS ==========
    
    reference_ch = Channel.fromPath(params.reference, checkIfExists: true).first()
    reference_fai_ch = Channel.fromPath("${params.reference}.fai", checkIfExists: true).first()
    
    // BAM filter config
    if (params.bam_filter_config) {
        bam_filter_config_ch = Channel.fromPath(params.bam_filter_config, checkIfExists: true).first()
    } else {
        bam_filter_config_ch = Channel.value(file('NO_FILE'))
    }
    
    // Ploidy map channel
    if (params.ploidy_map) {
        ploidy_map_ch = Channel.fromPath(params.ploidy_map, checkIfExists: true).first()
    } else {
        ploidy_map_ch = Channel.value(file('NO_FILE'))
    }
    
    // FreeBayes config
    if (params.freebayes_config) {
        config_ch = Channel.fromPath(params.freebayes_config, checkIfExists: true).first()
    } else {
        config_ch = Channel.value([])
    }
    
    // ========== VARIANT CALLING PIPELINE ==========
    
    log.info "Phase 1: BAM Processing and Filtering"
    
    // Process BAM files
    if (params.bam_filter_config) {
        bam_pairs_ch = Channel.fromPath(params.bams, checkIfExists: true)
            .map { bam ->
                def bai = file("${bam}.bai")
                if (!bai.exists()) {
                    error "BAM index file not found: ${bai}. Please run 'samtools index ${bam}'"
                }
                return tuple(bam, bai)
            }
        
        filtered_bams = FILTER_BAMS(
            bam_pairs_ch.map { it[0] },
            bam_pairs_ch.map { it[1] },
            bam_filter_config_ch
        )
        
        bam_bai_split = filtered_bams
            .toList()
            .multiMap { pairs ->
                bams: pairs.collect { it[0] }
                bais: pairs.collect { it[1] }
            }
        
        bam_files = bam_bai_split.bams
        bai_files = bam_bai_split.bais
    } else {
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
    
    log.info "Phase 2: Parallel Variant Calling"
    
    // Create genome chunks
    chunks = CREATE_CHUNKS(reference_ch, params.num_chunks)
    chunk_regions = chunks[1]
    
    // Parse chunks
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
    
    // Run FreeBayes on each chunk
    vcf_chunks = FREEBAYES_CHUNK(
        chunk_ch,
        reference_ch,
        reference_fai_ch,
        bam_files,
        bai_files,
        config_ch,
        ploidy_map_ch
    )
    
    // Combine VCFs
    all_vcfs = vcf_chunks.map { chunk_id, vcf -> vcf }.collect()
    combined_vcf = COMBINE_VCFS(all_vcfs, params.output_vcf)
    
    // Basic summary
    vcf_summary = SUMMARIZE_VCFS(combined_vcf.vcf)
    
    // ========== COMPREHENSIVE QC PIPELINE ==========
    
    if (params.generate_report) {
        
        log.info "Phase 3: Quality Control Analysis"
        
        // Extract variant statistics
        variant_stats_out = VARIANT_STATS(
            combined_vcf.vcf,
            combined_vcf.index
        )
        
        // Sample-level QC
        sample_qc_out = SAMPLE_QC(
            combined_vcf.vcf,
            combined_vcf.index,
            ploidy_map_ch
        )
        
        // Locus-level QC
        locus_qc_out = LOCUS_QC(
            combined_vcf.vcf,
            combined_vcf.index
        )
        
        // Contamination check (optional)
        if (params.check_contamination) {
            contamination_out = CONTAMINATION_CHECK(
                combined_vcf.vcf,
                combined_vcf.index,
                sample_qc_out.heterozygosity
            )
            contamination_report = contamination_out.report
        } else {
            contamination_report = Channel.value(file('NO_FILE'))
        }
        
        log.info "Phase 4: Visualization and Reporting"
        
        // Generate QC plots
        plots_out = GENERATE_PLOTS(
            variant_stats_out.stats,
            variant_stats_out.sample_stats,
            variant_stats_out.chr_stats,
            variant_stats_out.af_dist,
            variant_stats_out.depth_dist,
            variant_stats_out.qual_dist,
            sample_qc_out.metrics,
            sample_qc_out.heterozygosity,
            sample_qc_out.depth_stats,
            sample_qc_out.missingness,
            sample_qc_out.relatedness,
            locus_qc_out.metrics,
            sample_qc_out.pool_metrics
        )
        
        // Generate filter recommendations
        filter_rec_out = FILTER_RECOMMENDATIONS(
            variant_stats_out.stats,
            locus_qc_out.metrics,
            sample_qc_out.metrics,
            variant_stats_out.depth_dist,
            variant_stats_out.qual_dist,
            combined_vcf.vcf
        )
        
        // Validate filters (optional)
        if (params.validate_filters) {
            filter_validation_out = FILTER_VALIDATION(
                combined_vcf.vcf,
                combined_vcf.index,
                filter_rec_out.thresholds
            )
            filter_validation_report = filter_validation_out.report
        } else {
            filter_validation_report = Channel.value(file('NO_FILE'))
        }
        
        // Generate publication tables (optional)
        if (params.generate_publication_tables) {
            publication_out = PUBLICATION_TABLES(
                variant_stats_out.summary,
                variant_stats_out.sample_stats,
                variant_stats_out.chr_stats,
                filter_validation_report,
                contamination_report
            )
        }
        
        // Collect config files for report
        config_files_ch = Channel.empty()
        if (params.freebayes_config) {
            config_files_ch = config_files_ch.mix(
                Channel.fromPath(params.freebayes_config)
            )
        }
        if (params.bam_filter_config) {
            config_files_ch = config_files_ch.mix(
                Channel.fromPath(params.bam_filter_config)
            )
        }
        config_files_collected = config_files_ch.collect().ifEmpty([])
        
        // Pipeline trace
        pipeline_trace_ch = file("${params.output_dir}/pipeline/pipeline_trace.txt").exists() ?
            Channel.fromPath("${params.output_dir}/pipeline/pipeline_trace.txt") :
            Channel.value(file('NO_FILE'))
        
        // Generate comprehensive report
        final_report = GENERATE_REPORT(
            variant_stats_out.summary,
            filter_rec_out.recommendations,
            filter_rec_out.thresholds,
            filter_rec_out.comparison,
            sample_qc_out.metrics,
            plots_out.plots.first(),
            pipeline_trace_ch,
            config_files_collected,
            params
        )
        
        log.info "Phase 5: Report Generation Complete"
    }
}

workflow.onComplete {
    log.info """
    ╔════════════════════════════════════════════════════════════════╗
    ║                    Pipeline Execution Summary                  ║
    ╠════════════════════════════════════════════════════════════════╣
    ║ Status:      ${workflow.success ? '✅ SUCCESS' : '❌ FAILED'}
    ║ Completed:   ${workflow.complete}
    ║ Duration:    ${workflow.duration}
    ║ Exit status: ${workflow.exitStatus}
    ║ Error:       ${workflow.errorReport ?: 'None'}
    ╚════════════════════════════════════════════════════════════════╝
    """
    
    if (workflow.success) {
        log.info """
        📁 Output Files:
        ├── 📄 ${params.output_dir}/${params.output_vcf}
        """
        
        if (params.generate_report) {
            log.info """├── 📊 ${params.output_dir}/snp_calling_report.html
        ├── 📁 ${params.output_dir}/qc/
        │   ├── 📄 filter_commands.sh
        │   ├── 📄 filtering_recommendations.txt
        │   ├── 📁 plots/
        │   │   ├── sample_qc/
        │   │   ├── variant_qc/
        │   │   └── locus_qc/"""
            
            if (params.check_contamination) {
                log.info """│   └── 📁 contamination/
        │       └── 📄 contamination_report.tsv"""
            }
            
            if (params.validate_filters) {
                log.info """│   └── 📁 filter_validation/
        │       └── 📄 filter_validation_report.tsv"""
            }
            
            if (params.generate_publication_tables) {
                log.info """└── 📁 ${params.output_dir}/publication/
            ├── 📄 Table1_sample_summary.csv
            ├── 📄 Table2_variant_summary.csv
            └── 📄 publication_tables_latex.tex"""
            }
            
            log.info """
        
        🎯 Recommended Next Steps:
        1. Review the HTML report: ${params.output_dir}/snp_calling_report.html
        2. Apply filters: bash ${params.output_dir}/qc/filter_commands.sh
        3. Check contamination suspects if any were identified
        4. Annotate filtered variants with SnpEff/VEP
        """
        }
        
        log.info """
        📈 Pipeline Performance:
        ├── Chunks processed: ${params.num_chunks}
        ├── Execution time: ${workflow.duration}
        └── Resources: See ${params.output_dir}/pipeline/pipeline_report.html
        """
    } else {
        log.info """
        ❌ Pipeline failed. Please check:
        1. Error messages above
        2. Log files in work/ directory
        3. Ensure all input files are valid
        4. Check available system resources
        """
    }
}