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
    tuple val(chunk_id), path("${chunk_id}.*"), emit: chunk_files
    path "${chunk_id}.arg", emit: args_file
    
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
    
    # Verify indices exist
    echo "Verifying BAM indices:"
    for bam in *.bam; do
        if [ -f "\$bam" ]; then
            if [ -f "\${bam}.bai" ]; then
                echo "  \$bam: index present"
            else
                echo "  \$bam: ERROR - index missing!"
                exit 1
            fi
        fi
    done
    
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
    ANGSD_CMD="\$ANGSD_CMD -out ${chunk_id}"
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
    # Read the config file
    config_file = '${config_file}' if '${config_file}' else None
    if not config_file:
        # Try to find a JSON file in current directory
        import glob
        json_files = glob.glob('*.json')
        if json_files:
            config_file = json_files[0]
    
    if not config_file:
        raise FileNotFoundError("No config file found")
    
    with open(config_file) as f:
        config = json.load(f)
    
    # Extract sections based on structure
    basic = config.get('basic_options', {})
    filters = config.get('filters', {})
    gl_options = config.get('genotype_likelihoods', {})
    snp_calling = config.get('snp_calling', {})
    output_opts = config.get('output_options', {})
    pooled_opts = config.get('pooled_sequencing', {})
    
    params = []
    
    # ========== BASIC OPTIONS ==========
    if 'doCounts' in basic:
        params.append(f"-doCounts {int(basic['doCounts'])}")
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
    
    # ========== FILTERS ==========
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
    
    # ========== GENOTYPE LIKELIHOODS ==========
    gl_model = gl_options.get('GL', 1)
    params.append(f"-GL {gl_model}")
    print(f"echo 'Using genotype likelihood model: {gl_model}'")
    
    # Always output beagle format
    params.append("-doGlf 2")
    
    if gl_options.get('doGeno'):
        params.append(f"-doGeno {gl_options['doGeno']}")
    if gl_options.get('doPost'):
        params.append(f"-doPost {gl_options['doPost']}")
    if gl_options.get('geno_minDepth'):
        params.append(f"-geno_minDepth {gl_options['geno_minDepth']}")
    if gl_options.get('geno_minMM'):
        params.append(f"-geno_minMM {gl_options['geno_minMM']}")
    
    # ========== SNP CALLING ==========
    if snp_calling.get('doMajorMinor'):
        params.append(f"-doMajorMinor {snp_calling['doMajorMinor']}")
    if snp_calling.get('doMaf'):
        params.append(f"-doMaf {snp_calling['doMaf']}")
    if snp_calling.get('SNP_pval') is not None:
        params.append(f"-SNP_pval {snp_calling['SNP_pval']}")
    if snp_calling.get('rmTriallelic'):
        params.append(f"-rmTriallelic {int(snp_calling['rmTriallelic'])}")
    if snp_calling.get('skipTriallelic'):
        params.append(f"-skipTriallelic {int(snp_calling['skipTriallelic'])}")
    if snp_calling.get('minMaf'):
        params.append(f"-minMaf {snp_calling['minMaf']}")
    if snp_calling.get('minLRT'):
        params.append(f"-minLRT {snp_calling['minLRT']}")
    
    # ========== POOLED SEQUENCING ==========
    if pooled_opts.get('doSaf'):
        params.append(f"-doSaf {pooled_opts['doSaf']}")
        print("echo 'Site allele frequency (SAF) calculation enabled'")
    if pooled_opts.get('fold'):
        params.append(f"-fold {int(pooled_opts['fold'])}")
    if pooled_opts.get('anc'):
        params.append(f"-anc {pooled_opts['anc']}")
    
    # Always enable BCF output if genotype calling is enabled
    if gl_options.get('doGeno', 0) > 0:
        params.append("-dobcf 1")
        print("echo 'BCF output enabled (genotype calling is active)'")
    
    # Always add counts and depth
    if 'doCounts' not in basic:
        params.append("-doCounts 1")
    if 'doDepth' not in output_opts:
        params.append("-doDepth 1")
    
    # ========== OUTPUT OPTIONS ==========
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
            print(f"echo 'Pooled sequencing mode: {n_ind} individuals per pool'")
            params.append(f"-nInd {n_ind}")
    
except Exception as e:
    print(f"# Error parsing config: {e}", file=sys.stderr)
    print("ANGSD_PARAMS=''")
    print("echo 'WARNING: Failed to parse config file, using minimal defaults'")
    # Minimal default parameters - always output both formats
    print("ANGSD_PARAMS='-GL 1 -doGlf 2 -doMajorMinor 1 -doMaf 1 -doCounts 1 -doGeno 2 -doPost 1 -dobcf 1'")
PYTHON_CONFIG
        
        source angsd_params.sh
        
        if [ -n "\$ANGSD_PARAMS" ]; then
            ANGSD_CMD="\$ANGSD_CMD \$ANGSD_PARAMS"
        fi
    else
        echo "No config file - using default ANGSD parameters"
        # Default parameters - output both Beagle and BCF
        ANGSD_CMD="\$ANGSD_CMD -GL 1"           # SAMtools GL model
        ANGSD_CMD="\$ANGSD_CMD -doGlf 2"        # Output beagle format
        ANGSD_CMD="\$ANGSD_CMD -doMajorMinor 1" # Infer major/minor from GL
        ANGSD_CMD="\$ANGSD_CMD -doMaf 1"        # Calculate allele frequencies
        ANGSD_CMD="\$ANGSD_CMD -doCounts 1"     # Count alleles
        ANGSD_CMD="\$ANGSD_CMD -doDepth 1"      # Calculate depth
        ANGSD_CMD="\$ANGSD_CMD -doGeno 2"       # Call genotypes
        ANGSD_CMD="\$ANGSD_CMD -doPost 1"       # Calculate posteriors
        ANGSD_CMD="\$ANGSD_CMD -dobcf 1"        # Output BCF format
        
        # Basic filters
        ANGSD_CMD="\$ANGSD_CMD -minMapQ 20"
        ANGSD_CMD="\$ANGSD_CMD -minQ 20"
        ANGSD_CMD="\$ANGSD_CMD -uniqueOnly 1"
        ANGSD_CMD="\$ANGSD_CMD -remove_bads 1"
        ANGSD_CMD="\$ANGSD_CMD -only_proper_pairs 1"
    fi
    
    # Run ANGSD
    echo "Full ANGSD command: \$ANGSD_CMD"
    echo "----------------------------------------"
    
    # Actually run ANGSD and capture exit status
    \$ANGSD_CMD
    ANGSD_EXIT=\$?
    
    echo "----------------------------------------"
    echo "ANGSD exit status: \$ANGSD_EXIT"
    
    # Check if ANGSD ran successfully
    if [ \$ANGSD_EXIT -eq 0 ] && [ -f "${chunk_id}.arg" ]; then
        echo "ANGSD completed successfully"
        echo "Arguments saved in ${chunk_id}.arg"
        
        # List output files
        echo "Output files generated:"
        ls -la ${chunk_id}.* 2>/dev/null || echo "No output files found"
        
        # Check if beagle file has content
        if [ -f "${chunk_id}.beagle.gz" ]; then
            SITES=\$(zcat ${chunk_id}.beagle.gz 2>/dev/null | tail -n +2 | wc -l)
            echo "Sites with genotype likelihoods: \$SITES"
            
            if [ \$SITES -eq 0 ]; then
                echo "WARNING: Beagle file exists but contains no sites"
            fi
        else
            echo "WARNING: No beagle file generated"
        fi
        
        # Check if BCF was generated
        if [ -f "${chunk_id}.bcf" ]; then
            VARIANTS=\$(bcftools view -H ${chunk_id}.bcf 2>/dev/null | wc -l || echo "0")
            echo "Variants in BCF: \$VARIANTS"
            
            if [ \$VARIANTS -eq 0 ]; then
                echo "WARNING: BCF file exists but contains no variants"
            fi
        else
            echo "WARNING: No BCF file generated"
        fi
        
        # If no sites were found, this might be okay for some chunks
        if [ ! -f "${chunk_id}.beagle.gz" ] || [ ! -f "${chunk_id}.mafs.gz" ]; then
            echo "NOTE: This chunk may have no variable sites, creating empty files"
            touch ${chunk_id}.beagle.gz
            touch ${chunk_id}.mafs.gz
            
            if [ ! -f "${chunk_id}.bcf" ]; then
                # Create minimal BCF
                echo "Creating minimal BCF file"
                echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" > ${chunk_id}.vcf
                bcftools view -Ob ${chunk_id}.vcf > ${chunk_id}.bcf
                rm ${chunk_id}.vcf
            fi
        fi
        
    else
        echo "ERROR: ANGSD failed with exit status \$ANGSD_EXIT"
        echo "No .arg file was generated"
        echo ""
        echo "Checking for error indicators:"
        
        # Check if BAMs are readable
        echo "BAM file check:"
        while read bam; do
            if [ -f "\$bam" ]; then
                echo "  \$bam: exists"
                samtools quickcheck "\$bam" && echo "    - valid" || echo "    - CORRUPT"
            else
                echo "  \$bam: NOT FOUND"
            fi
        done < bam.list
        
        # Check reference
        echo ""
        echo "Reference check:"
        if [ -f "${reference}" ]; then
            echo "  ${reference}: exists"
        else
            echo "  ${reference}: NOT FOUND"
        fi
        
        # Check regions
        echo ""
        echo "Regions to process:"
        cat angsd_regions.txt
        
        # Create empty outputs to allow pipeline to continue
        echo ""
        echo "Creating empty output files to allow pipeline to continue..."
        touch ${chunk_id}.beagle.gz
        touch ${chunk_id}.mafs.gz
        
        # Create empty but valid BCF
        echo '##fileformat=VCFv4.2' > ${chunk_id}.vcf
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> ${chunk_id}.vcf
        bcftools view -Ob ${chunk_id}.vcf > ${chunk_id}.bcf 2>/dev/null || touch ${chunk_id}.bcf
        rm -f ${chunk_id}.vcf
        
        # Exit with error
        exit 1
    fi
    
    echo "Chunk ${chunk_id} processing complete"
    """
}