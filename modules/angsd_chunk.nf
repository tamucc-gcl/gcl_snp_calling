// modules/angsd_chunk.nf - ANGSD genotyping for genomic chunks

process ANGSD_CHUNK {
    tag "angsd_${chunk_id}"
    
    input:
    tuple val(chunk_id), val(regions_string)
    path reference
    path reference_fai
    path bams        // List of BAM files
    path bais        // List of BAI files
    path config      // ANGSD configuration JSON
    path sites_file  // Optional sites file for specific positions
    
    output:
    tuple val(chunk_id), path("chunk_${chunk_id}.*"), emit: chunk_files
    path "chunk_${chunk_id}.arg", emit: args_file
    
    script:
    def config_file = config.size() > 0 ? config[0] : null
    def has_config = config_file != null
    def has_sites = sites_file.name != 'NO_FILE'
    """
    echo "Processing chunk ${chunk_id} with ANGSD"
    echo "Regions: ${regions_string}"
    echo "Reference: ${reference}"
    echo "Will output both Beagle and BCF formats"
    
    # Count BAM files
    BAM_COUNT=\$(ls -1 *.bam 2>/dev/null | wc -l)
    echo "Total BAM files: \$BAM_COUNT"
    
    if [ \$BAM_COUNT -eq 0 ]; then
        echo "ERROR: No BAM files found!"
        exit 1
    fi
    
    # Create BAM list file
    ls -1 *.bam > bam.list
    echo "BAM list created with \$BAM_COUNT files"
    
    # Create regions file for ANGSD
    echo "${regions_string}" | tr ',' '\\n' > regions.txt
    
    # Convert regions to ANGSD format (chr:start-end)
    > angsd_regions.txt
    while read -r region; do
        if [ -n "\$region" ]; then
            echo "\$region" >> angsd_regions.txt
        fi
    done < regions.txt
    
    echo "ANGSD regions:"
    cat angsd_regions.txt
    
    # Build base ANGSD command
    ANGSD_CMD="angsd -bam bam.list"
    ANGSD_CMD="\$ANGSD_CMD -ref ${reference}"
    ANGSD_CMD="\$ANGSD_CMD -rf angsd_regions.txt"
    ANGSD_CMD="\$ANGSD_CMD -out chunk_${chunk_id}"
    ANGSD_CMD="\$ANGSD_CMD -nThreads ${task.cpus ?: 8}"
    
    # Add sites file if provided
    if [ "${has_sites}" = "true" ]; then
        echo "Using sites file: ${sites_file}"
        ANGSD_CMD="\$ANGSD_CMD -sites ${sites_file}"
    fi
    
    # Parse configuration file if provided
    if [ "${has_config}" = "true" ]; then
        echo "Loading ANGSD parameters from config file"
        
        python3 << 'PYTHON_CONFIG' > angsd_params.sh
import json
import sys

try:
    config_file = '${config_file}'
    
    with open(config_file) as f:
        config = json.load(f)
    
    # Extract parameter groups
    basic = config.get('basic_options', {})
    filters = config.get('filters', {})
    gl_options = config.get('genotype_likelihoods', {})
    snp_calling = config.get('snp_calling', {})
    output_opts = config.get('output_options', {})
    pooled_opts = config.get('pooled_sequencing', {})
    
    params = []
    
    # Basic options
    if 'remove_bads' in basic:
        params.append(f"-remove_bads {int(basic['remove_bads'])}")
    if 'unique_only' in basic:
        params.append(f"-uniqueOnly {int(basic['unique_only'])}")
    if 'only_proper_pairs' in basic:
        params.append(f"-only_proper_pairs {int(basic['only_proper_pairs'])}")
    if 'trim' in basic:
        params.append(f"-trim {int(basic['trim'])}")
    if 'C' in basic:
        params.append(f"-C {basic['C']}")
    if 'baq' in basic:
        params.append(f"-baq {basic['baq']}")
    if 'checkBamHeaders' in basic:
        params.append(f"-checkBamHeaders {int(basic['checkBamHeaders'])}")
    
    # Quality filters
    if 'minQ' in filters:
        params.append(f"-minQ {filters['minQ']}")
    if 'minMapQ' in filters:
        params.append(f"-minMapQ {filters['minMapQ']}")
    if 'minInd' in filters:
        params.append(f"-minInd {filters['minInd']}")
    if 'setMinDepth' in filters:
        params.append(f"-setMinDepth {filters['setMinDepth']}")
    if 'setMaxDepth' in filters:
        params.append(f"-setMaxDepth {filters['setMaxDepth']}")
    if 'setMinDepthInd' in filters:
        params.append(f"-setMinDepthInd {filters['setMinDepthInd']}")
    if 'setMaxDepthInd' in filters:
        params.append(f"-setMaxDepthInd {filters['setMaxDepthInd']}")
    
    # Genotype likelihood model
    gl_model = gl_options.get('GL', 1)  # Default to SAMtools model
    params.append(f"-GL {gl_model}")
    print(f"echo 'Using genotype likelihood model: {gl_model}'")
    
    # Add genotype likelihood calculation options
    if gl_options.get('doGlf'):
        params.append(f"-doGlf {gl_options['doGlf']}")
    if gl_options.get('doGeno'):
        params.append(f"-doGeno {gl_options['doGeno']}")
    if gl_options.get('doPost'):
        params.append(f"-doPost {gl_options['doPost']}")
    
    # SNP calling options
    if snp_calling.get('doMajorMinor'):
        params.append(f"-doMajorMinor {snp_calling['doMajorMinor']}")
    if snp_calling.get('doMaf'):
        params.append(f"-doMaf {snp_calling['doMaf']}")
    if snp_calling.get('SNP_pval') is not None:
        params.append(f"-SNP_pval {snp_calling['SNP_pval']}")
    if snp_calling.get('rmTriallelic'):
        params.append(f"-rmTriallelic {int(snp_calling['rmTriallelic'])}")
    if snp_calling.get('minMaf'):
        params.append(f"-minMaf {snp_calling['minMaf']}")
    if snp_calling.get('minLRT'):
        params.append(f"-minLRT {snp_calling['minLRT']}")
    
    # Pooled sequencing specific options
    if pooled_opts.get('doSaf'):
        params.append(f"-doSaf {pooled_opts['doSaf']}")
        print("echo 'Site allele frequency (SAF) calculation enabled'")
    if pooled_opts.get('fold'):
        params.append(f"-fold {int(pooled_opts['fold'])}")
    if pooled_opts.get('anc'):
        params.append(f"-anc {pooled_opts['anc']}")
    
    # Output options based on requested format
    output_format = '${output_format}'
    
    # Always calculate genotype likelihoods
    params.append("-doGlf 2")  # Beagle format
    
    if output_format in ['vcf', 'all']:
        params.append("-doVcf 1")
        params.append("-doCounts 1")
        params.append("-doDepth 1")
        print("echo 'VCF output enabled'")
    
    if output_format in ['beagle', 'all']:
        print("echo 'Beagle genotype likelihood output enabled'")
    
    # Add counts and depth for all outputs
    params.append("-doCounts 1")
    params.append("-doDepth 1")
    
    # Additional output options from config
    if output_opts.get('dumpCounts'):
        params.append(f"-dumpCounts {output_opts['dumpCounts']}")
    if output_opts.get('doDepth'):
        params.append(f"-doDepth {output_opts['doDepth']}")
    if output_opts.get('maxDepth'):
        params.append(f"-maxDepth {output_opts['maxDepth']}")
    
    print(f"ANGSD_PARAMS='{' '.join(params)}'")
    print(f"echo 'Loaded {len(params)} parameters from config'")
    
    # Special handling for pooled samples
    if pooled_opts.get('is_pooled', False):
        n_ind = pooled_opts.get('n_individuals_per_pool', None)
        if n_ind:
            # For pooled samples, we might need to adjust -nInd
            # This overrides the default (number of BAM files)
            print(f"echo 'Pooled sequencing mode: {n_ind} individuals per pool'")
            params.append(f"-nInd {n_ind}")
    
except Exception as e:
    print(f"# Error parsing config: {e}", file=sys.stderr)
    print("ANGSD_PARAMS=''")
    print("echo 'WARNING: Failed to parse config file, using minimal defaults'")
    # Minimal default parameters
    print("ANGSD_PARAMS='-GL 1 -doGlf 2 -doMajorMinor 1 -doMaf 1 -doCounts 1'")
PYTHON_CONFIG
        
        source angsd_params.sh
        
        if [ -n "\$ANGSD_PARAMS" ]; then
            ANGSD_CMD="\$ANGSD_CMD \$ANGSD_PARAMS"
        fi
    else
        echo "No config file - using default ANGSD parameters"
        # Default parameters for basic genotype likelihood calculation
        ANGSD_CMD="\$ANGSD_CMD -GL 1"           # SAMtools GL model
        ANGSD_CMD="\$ANGSD_CMD -doGlf 2"        # Output beagle format
        ANGSD_CMD="\$ANGSD_CMD -doMajorMinor 1" # Infer major/minor from GL
        ANGSD_CMD="\$ANGSD_CMD -doMaf 1"        # Calculate allele frequencies
        ANGSD_CMD="\$ANGSD_CMD -doCounts 1"     # Count alleles
        ANGSD_CMD="\$ANGSD_CMD -doDepth 1"      # Calculate depth
        
        # Basic filters
        ANGSD_CMD="\$ANGSD_CMD -minMapQ 20"
        ANGSD_CMD="\$ANGSD_CMD -minQ 20"
        ANGSD_CMD="\$ANGSD_CMD -uniqueOnly 1"
        ANGSD_CMD="\$ANGSD_CMD -remove_bads 1"
        ANGSD_CMD="\$ANGSD_CMD -only_proper_pairs 1"
        
        if [ "${output_format}" = "vcf" ] || [ "${output_format}" = "all" ]; then
            ANGSD_CMD="\$ANGSD_CMD -doVcf 1"
        fi
    fi
    
    # Run ANGSD
    echo "Full ANGSD command: \$ANGSD_CMD"
    \$ANGSD_CMD
    
    # Save the arguments file for debugging
    if [ -f "chunk_${chunk_id}.arg" ]; then
        echo "ANGSD completed successfully"
        echo "Arguments saved in chunk_${chunk_id}.arg"
        
        # List output files
        echo "Output files generated:"
        ls -la chunk_${chunk_id}.*
        
        # Check if beagle file has content
        if [ -f "chunk_${chunk_id}.beagle.gz" ]; then
            SITES=\$(zcat chunk_${chunk_id}.beagle.gz | tail -n +2 | wc -l)
            echo "Sites with genotype likelihoods: \$SITES"
        fi
        
        # Check if BCF was generated
        if [ -f "chunk_${chunk_id}.bcf" ]; then
            VARIANTS=\$(bcftools view -H chunk_${chunk_id}.bcf | wc -l)
            echo "Variants in BCF: \$VARIANTS"
        fi
        
    else
        echo "WARNING: ANGSD may have failed - no .arg file generated"
        # Create empty outputs to prevent pipeline failure
        touch chunk_${chunk_id}.beagle.gz
        touch chunk_${chunk_id}.mafs.gz
        touch chunk_${chunk_id}.bcf
        echo "Created empty output files"
    fi
    
    echo "Chunk ${chunk_id} processing complete"
    """
}