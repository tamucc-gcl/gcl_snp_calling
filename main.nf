#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Parameters
params.bams = "*.bam"
params.reference = "reference.fasta"
params.num_chunks = 10
params.output_dir = "results"
params.output_vcf = "raw_variants.vcf.gz"
params.freebayes_config = null
params.angsd_config = null
params.bam_filter_config = null
params.ploidy_map = null
params.sites_file = null
params.genotyper = "freebayes"  // "freebayes" or "angsd"

// ##### CHANGE 1: data type, decomposition and diagnostic parameters #####
// data_type controls the Ts/Tv sweep grids and the PCA site-spacing default.
// Reduced-representation and whole-genome data need different values: ddRAD
// has structural locus dropout, so call rates below 50% are normal, while WGS
// at high sample counts should not, and its dominant false-positive class is
// collapsed-repeat pileups rather than dropout.
params.data_type     = "ddrad"   // "ddrad" or "wgs"
params.atom_overlaps = "."       // "." = missing GT on overlap, "*" = star allele
params.tstv_ns_grid   = null     // null = use the data_type default grid
params.tstv_maf_grid  = null
params.tstv_min_sites = 200
params.pca_ns_frac     = 0.50
params.pca_min_spacing = null    // null = 5000 for wgs, 500 for ddrad
// ##### END CHANGE 1 #####

// PCAngsd parameters
params.pcangsd_minMaf = null //null to skip and default to all loci which pass ANGSD SNP calling filters, or set to a value (e.g., 0.05) to filter by MAF
params.pcangsd_maf_iter = 10000
params.pcangsd_iter = 10000
params.pcangsd_eigenvalues = null //null for automatic selection with MAP values

// Import modules
include { CREATE_CHUNKS } from './modules/create_chunks'
include { FILTER_BAMS } from './modules/filter_bams'
include { FREEBAYES_CHUNK } from './modules/freebayes_chunk'
include { ANGSD_CHUNK } from './modules/angsd_chunk'
include { COMBINE_VCFS } from './modules/combine_vcfs'
include { COMBINE_ANGSD } from './modules/combine_angsd'
include { PCANGSD } from './modules/pcangsd'
include { SUMMARIZE_VCFS } from './modules/summarize_vcfs'
include { samtools_stats as SAMTOOLS_STATS_RAW } from './modules/samtools_stats'
include { samtools_stats as SAMTOOLS_STATS_FILTERED } from './modules/samtools_stats'
include { multiqc as MULTIQC_RAW_BAMS } from './modules/multiqc'
include { multiqc as MULTIQC_FILTERED_BAMS } from './modules/multiqc'
include { REPORT_SNP_CALLING_SUMMARY } from './modules/create_report.nf'

// ##### CHANGE 2: new module imports #####
include { DECOMPOSE_VCF } from './modules/decompose_vcf'
include { TSTV_SWEEP }    from './modules/tstv_sweep'
// ##### END CHANGE 2 #####


def helpMessage() {
    log.info"""
    ================================================================
    Parallel Genotyping Pipeline (FreeBayes/ANGSD)
    ================================================================
    
    Usage:
    nextflow run main.nf --bams "path/to/*.bam" --reference "reference.fasta"
    
    Required parameters:
    --bams              Path to BAM files (glob pattern, e.g., "*.bam")
    --reference         Path to reference genome FASTA file
    
    Genotyping method:
    --genotyper         Choose genotyping method: "freebayes" or "angsd" (default: "freebayes")
    
    Optional parameters:
    --num_chunks        Number of chunks to split genome into (default: 10)
    --output_dir        Output directory (default: "results")
    --output_vcf        Name of output file (default: "raw_variants.vcf.gz")
                        - FreeBayes: Creates this exact VCF file
                        - ANGSD: Uses as base name (e.g., "out.vcf.gz" → "out.beagle.gz", "out.mafs.gz")
    --bam_filter_config Path to JSON file containing BAM filter parameters (optional)
    --data_type         "ddrad" or "wgs" (default: "ddrad"). Sets the Ts/Tv sweep
                        grids and the PCA site-spacing default.

    FreeBayes-specific:
    --freebayes_config  Path to JSON file containing FreeBayes parameters (optional)
    --ploidy_map        Path to file mapping BAM files to ploidy values (optional)
    --atom_overlaps     "." or "*" (default: "."). How bcftools norm handles
                        overlapping atomized alleles during decomposition.

    Diagnostics (FreeBayes only):
    --tstv_ns_grid      Comma-separated call-rate fractions for the Ts/Tv sweep
                        (default: set by --data_type)
    --tstv_maf_grid     Comma-separated MAF floors for the Ts/Tv sweep
                        (default: set by --data_type)
    --tstv_min_sites    Grid cells below this site count are flagged unreliable (default: 200)
    --pca_ns_frac       PCA site call-rate floor as a fraction of n_samples (default: 0.50)
    --pca_min_spacing   Minimum bp between PCA sites; 0 disables (default: by --data_type)

    ANGSD-specific:
    --angsd_config      Path to JSON file containing ANGSD parameters (optional)
    --sites_file        Path to sites file for ANGSD (optional, for targeted analysis)
    
    Examples:
    # Run with FreeBayes (default)
    nextflow run main.nf --bams "*.bam" --reference genome.fa
    
    # Run with ANGSD for genotype likelihoods
    nextflow run main.nf --bams "*.bam" --reference genome.fa --genotyper angsd
    
    # Run ANGSD with custom config and consistent naming
    nextflow run main.nf --bams "*.bam" --reference genome.fa \\
        --genotyper angsd \\
        --angsd_config angsd_parameters.json \\
        --output_vcf "my_project.vcf.gz"
    # Creates: my_project.beagle.gz, my_project.mafs.gz, etc.
    
    # Run FreeBayes with same naming
    nextflow run main.nf --bams "*.bam" --reference genome.fa \\
        --genotyper freebayes \\
        --output_vcf "my_project.vcf.gz"
    # Creates: my_project.vcf.gz and my_project.decomposed.vcf.gz

    # Whole-genome data instead of ddRAD
    nextflow run main.nf --bams "*.bam" --reference genome.fa \\
        --data_type wgs
    
    # ANGSD with VCF output only
    nextflow run main.nf --bams "*.bam" --reference genome.fa \\
        --genotyper angsd \\
        --angsd_output_format vcf \\
        --output_vcf "only_vcf.vcf.gz"
    # Creates: only_vcf.vcf.gz (VCF format only, no beagle)
    
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    exit 0
}

if (!params.bams || !params.reference) {
    log.error "ERROR: Please provide both --bams and --reference parameters"
    helpMessage()
    exit 1
}

// ##### CHANGE 3: validate data_type #####
if (!(params.data_type in ['ddrad', 'wgs'])) {
    log.error "ERROR: Invalid data_type '${params.data_type}'. Must be 'ddrad' or 'wgs'"
    exit 1
}
// ##### END CHANGE 3 #####

workflow {
    // Derive output prefix from output_vcf parameter
    // Remove .vcf.gz or .vcf extensions to get base name
    def output_prefix = params.output_vcf.replaceAll(/\.(vcf|VCF)(\.gz)?$/, "")
    
    log.info "================================================================"
    log.info "Parallel Genotyping Pipeline"
    log.info "================================================================"
    log.info "Genotyper:         ${params.genotyper}"
    log.info "Data type:         ${params.data_type}"   // ##### CHANGE 4 #####
    log.info "BAM files:         ${params.bams}"
    log.info "Reference:         ${params.reference}"
    log.info "Number of chunks:  ${params.num_chunks}"
    log.info "Output directory:  ${params.output_dir}"
    
    if (params.genotyper == "freebayes") {
        log.info "Output VCF:        ${params.output_vcf}"
        log.info "Decomposed VCF:    ${output_prefix}.decomposed.vcf.gz"
        log.info "FreeBayes config:  ${params.freebayes_config ?: 'Using default parameters'}"
        log.info "Ploidy map:        ${params.ploidy_map ?: 'Using global ploidy from config or default'}"
    } else if (params.genotyper == "angsd") {
        log.info "Output base name:  ${output_prefix}"
        log.info "ANGSD config:      ${params.angsd_config ?: 'Using default parameters'}"
        log.info "Sites file:        ${params.sites_file ?: 'None - will discover sites'}"
        log.info "Output files:      Beagle + VCF formats"
    } else {
        log.error "ERROR: Invalid genotyper '${params.genotyper}'. Must be 'freebayes' or 'angsd'"
        exit 1
    }
    
    log.info "BAM filter config: ${params.bam_filter_config ?: 'Using default parameters'}"
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
    
    // Simple channel pattern that caches properly
    raw_bam_ch = Channel.fromPath(params.bams, checkIfExists: true)
        .map { bam ->
            def bai = file("${bam}.bai")
            if (!bai.exists()) {
                error "BAM index file not found: ${bai}"
            }
            tuple(bam.simpleName, bam, bai)
        }
    
    // Run samtools stats on raw BAMs
    raw_bam_stats = SAMTOOLS_STATS_RAW(raw_bam_ch)
    
    // Collect stats for MultiQC
    raw_stats_files = raw_bam_stats
        .map { sid, stats, flagstats -> [stats, flagstats] }
        .flatten()
        .collect()
    
    MULTIQC_RAW_BAMS(raw_stats_files, Channel.value('raw_bams'))
    
    // Process BAM files for filtering or direct use
    if (params.bam_filter_config) {
        // Separate channel for filtering
        filter_ch = Channel.fromPath(params.bams, checkIfExists: true)
            .map { bam -> tuple(bam, file("${bam}.bai")) }
        
        filtered_bams = FILTER_BAMS(
            filter_ch.map { bam, bai -> bam },
            filter_ch.map { bam, bai -> bai },
            bam_filter_config_ch
        )
        
        filtered_bam_stats = SAMTOOLS_STATS_FILTERED(
            filtered_bams.map { bam, bai -> 
                tuple(bam.simpleName.replace('.filtered', ''), bam, bai)
            }
        )
        
        filtered_stats_files = filtered_bam_stats
            .map { sid, stats, flagstats -> [stats, flagstats] }
            .flatten()
            .collect()
        
        MULTIQC_FILTERED_BAMS(filtered_stats_files, Channel.value('filtered_bams'))
        
        // Collect for genotyping
        bam_files_ch = filtered_bams.map { bam, bai -> bam }.collect()
        bai_files_ch = filtered_bams.map { bam, bai -> bai }.collect()
        
    } else {
        // Use raw BAMs
        raw_for_genotyping = Channel.fromPath(params.bams, checkIfExists: true)
        bam_files_ch = raw_for_genotyping.collect()
        bai_files_ch = raw_for_genotyping.map { bam -> file("${bam}.bai") }.collect()
    }
    
    // Step 1: Create genome chunks
    chunks = CREATE_CHUNKS(reference_ch, params.num_chunks)
    chunk_regions = chunks[1]
    
    // Step 2: Parse chunks
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
    
    // Step 3: Run genotyping based on selected method
    if (params.genotyper == "freebayes") {
        // FreeBayes workflow
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
        
        vcf_chunks = FREEBAYES_CHUNK(
            chunk_ch,
            reference_ch,
            reference_fai_ch,
            bam_files_ch,
            bai_files_ch,
            config_ch,
            ploidy_map_ch
        )
        
        // ##### CHANGE 5 starts here #####

        // Step 4: Combine VCFs. This publishes the FreeBayes-native callset
        // unchanged, retaining multiallelic, MNP and complex records. It is the
        // only place the collapsed-paralog signal survives per-site, and the
        // only callset that can be re-decomposed under different settings.
        all_vcfs = vcf_chunks.map { chunk_id, vcf -> vcf }.collect()
        COMBINE_VCFS(all_vcfs, params.output_vcf)

        // Step 5: Decompose once on the combined callset.
        // Measured at 8 minutes for 5.7M records x 874 samples, so there is no
        // case for fanning this out per chunk. Running after the concat also
        // avoids the boundary shifts that norm's left-alignment creates at
        // chunk edges: 15% of records were realigned on this dataset.
        decomposed_name = "${output_prefix}.decomposed.vcf.gz"

        DECOMPOSE_VCF(
            COMBINE_VCFS.out.vcf,
            COMBINE_VCFS.out.index,
            reference_ch,
            reference_fai_ch,
            decomposed_name
        )

        // Step 6: Ts/Tv versus stringency diagnostic.
        // Emits no filtered VCF on purpose. A MAF floor suits structure and PCA
        // and is destructive to anything reading the site frequency spectrum,
        // so the thresholds cannot be chosen once at pipeline level.
        TSTV_SWEEP(
            DECOMPOSE_VCF.out.vcf,
            DECOMPOSE_VCF.out.index,
            Channel.value(params.data_type)
        )

        // Step 7: Summarize the decomposed callset
        SUMMARIZE_VCFS(
            DECOMPOSE_VCF.out.vcf,
            ploidy_map_ch,
            Channel.value('freebayes')
        )

        report_raw_stats_ch = COMBINE_VCFS.out.raw_stats
        report_sweep_ch     = TSTV_SWEEP.out.sweep
        report_precomp_ch   = DECOMPOSE_VCF.out.summary
        report_decomp_link  = Channel.value(decomposed_name)

        // ##### CHANGE 5 ends here #####
        
    } else if (params.genotyper == "angsd") {
        // ANGSD workflow
        if (params.angsd_config) {
            angsd_config_ch = Channel.fromPath(params.angsd_config, checkIfExists: true).first()
        } else {
            angsd_config_ch = Channel.value([])
        }
        
        if (params.sites_file) {
            sites_file_ch = Channel.fromPath(params.sites_file, checkIfExists: true).first()
        } else {
            sites_file_ch = Channel.value(file('NO_FILE'))
        }
        
        angsd_chunks = ANGSD_CHUNK(
            chunk_ch,
            reference_ch,
            reference_fai_ch,
            bam_files_ch,
            bai_files_ch,
            angsd_config_ch,
            sites_file_ch
        )
        
        // Collect all chunk output files
        all_angsd_files = angsd_chunks.chunk_files
            .map { chunk_id, files -> files }
            .flatten()
            .collect()
        
        // Combine ANGSD outputs
        COMBINE_ANGSD(
            all_angsd_files,
            Channel.value(output_prefix)
        )

        // Run PCAngsd on combined Beagle file
        PCANGSD(
            COMBINE_ANGSD.out.beagle,
            Channel.value(output_prefix)
        )
        
        // Optional: Run VCF summarize
        if (params.ploidy_map) {
            ploidy_map_ch = Channel.fromPath(params.ploidy_map, checkIfExists: true).first()
        } else {
            ploidy_map_ch = Channel.value(file('NO_FILE'))
        }
        //SUMMARIZE_VCFS(COMBINE_ANGSD.out.vcf, ploidy_map_ch)
        SUMMARIZE_VCFS(COMBINE_ANGSD.out.vcf, ploidy_map_ch, Channel.value('angsd'))

        // ##### CHANGE 6: placeholders so the report call is uniform #####
        // ANGSD emits simple biallelic calls, so decomposition is a no-op and
        // the Ts/Tv sweep has no raw-versus-decomposed comparison to make. The
        // report detects these sentinels and omits the corresponding sections.
        report_raw_stats_ch = Channel.value(file('NO_FILE'))
        report_sweep_ch     = Channel.value(file('NO_FILE'))
        report_precomp_ch   = Channel.value(file('NO_FILE'))
        report_decomp_link  = Channel.value('NO_DECOMPOSED_VCF')
        // ##### END CHANGE 6 #####
    }

    // ##### CHANGE 7: report call now takes seven inputs #####
    REPORT_SNP_CALLING_SUMMARY(
        Channel.value(params.genotyper),
        Channel.value("${output_prefix}.vcf.gz"),
        report_decomp_link,
        report_raw_stats_ch,
        report_sweep_ch,
        report_precomp_ch,
        SUMMARIZE_VCFS.out.report_inputs
    )
    // ##### END CHANGE 7 #####
}

workflow.onComplete {
    // Derive output prefix for completion message
    def output_prefix = params.output_vcf.replaceAll(/\.(vcf|VCF)(\.gz)?$/, "")
    
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
        
        // ##### CHANGE 8 starts here #####
        if (params.genotyper == "freebayes") {
            log.info "Callsets (both published):"
            log.info "- Raw, FreeBayes-native:  ${params.output_dir}/${params.output_vcf}"
            log.info "    Multiallelic and MNP records intact. Use for the"
            log.info "    collapsed-paralog signal and for re-decomposition."
            log.info "- Decomposed, analysis:   ${params.output_dir}/${output_prefix}.decomposed.vcf.gz"
            log.info "    Split, atomized, deduplicated, tags refilled."
            log.info "    Filter from this one."
            log.info ""
            log.info "Read this before choosing filtering thresholds:"
            log.info "  ${params.output_dir}/snp_qc/${output_prefix}.tstv_sweep.txt"
        } else if (params.genotyper == "angsd") {
        // ##### CHANGE 8 ends here #####
            log.info "ANGSD outputs:"
            log.info "- Genotype likelihoods: ${params.output_dir}/${output_prefix}.beagle.gz"
            log.info "- Allele frequencies: ${params.output_dir}/${output_prefix}.mafs.gz"
            log.info "- VCF file: ${params.output_dir}/${output_prefix}.vcf.gz"
            log.info "- Summary: ${params.output_dir}/${output_prefix}_summary.txt"
        }
        
        log.info ""
        log.info "QC Reports:"
        log.info "- Raw BAMs MultiQC: ${params.output_dir}/multiqc_reports/multiqc_raw_bams.html"
        if (params.bam_filter_config) {
            log.info "- Filtered BAMs MultiQC: ${params.output_dir}/multiqc_reports/multiqc_filtered_bams.html"
        }
        
        log.info ""
        log.info "Pipeline used exactly ${params.num_chunks} chunks with complete genome coverage"
        
    } else {
        log.info ""
        log.info "Pipeline failed. Check the error messages above."
    }
    log.info "================================================================"
}
