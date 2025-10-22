// modules/freebayes_chunk.nf - Complete version with fixed ploidy map processing

process FREEBAYES_CHUNK {
    
    input:
    tuple val(chunk_id), val(regions_string)
    path reference
    path reference_fai
    path bams        // List of BAM files
    path bais        // List of BAI files
    path config
    path ploidy_map  // Optional ploidy map file
    
    output:
    tuple val(chunk_id), path("chunk_${chunk_id}.vcf.gz")
    
    script:
    def config_file = config.size() > 0 ? config[0] : null
    def has_config = config_file != null
    def ploidy_file = ploidy_map.name != 'NO_FILE' ? ploidy_map : null
    def has_ploidy_map = ploidy_file != null
    """
    echo "Processing chunk ${chunk_id}"
    echo "Regions: ${regions_string}"
    echo "Reference: ${reference}"
    echo "Reference FAI: ${reference_fai}"
    echo "Config provided: ${has_config}"
    echo "Ploidy map provided: ${has_ploidy_map}"
    
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
    
    # Create BED file for regions
    echo "${regions_string}" | tr ',' '\\n' > regions.txt
    
    > chunk.bed
    while read -r region; do
        if [ -n "\$region" ]; then
            chrom=\$(echo "\$region" | cut -d':' -f1)
            pos=\$(echo "\$region" | cut -d':' -f2)
            start=\$(echo "\$pos" | cut -d'-' -f1)
            end=\$(echo "\$pos" | cut -d'-' -f2)
            bed_start=\$((start - 1))
            echo -e "\$chrom\\t\$bed_start\\t\$end" >> chunk.bed
        fi
    done < regions.txt
    
    echo "BED file contents:"
    cat chunk.bed
    
    # Process ploidy map if provided to create CNV map for FreeBayes
    if [ "${has_ploidy_map}" = "true" ]; then
        echo "================================================================"
        echo "Creating CNV map from ploidy map using Python"
        echo "================================================================"
        
        # Use Python for reliable parsing
        python3 << 'PYTHON_CNV' > cnv_map.txt
import sys
import subprocess

# Get sample names from BAM files
bam_samples = {}
sample_to_bam = {}

# List all BAM files
import glob
bam_files = glob.glob('*.bam')

print("Step 1: Inventorying BAM files and sample names...", file=sys.stderr)
print("------------------------------------------------", file=sys.stderr)

for bam in bam_files:
    # Extract sample name from BAM header
    try:
        result = subprocess.run(['samtools', 'view', '-H', bam], 
                              capture_output=True, text=True)
        for line in result.stdout.split('\\n'):
            if line.startswith('@RG'):
                # Look for SM tag
                parts = line.split('\\t')
                for part in parts:
                    if part.startswith('SM:'):
                        sm_tag = part[3:]
                        break
                else:
                    # No SM tag found, derive from filename
                    sm_tag = bam.replace('.filtered.bam', '').replace('.bam', '')
                
                bam_samples[bam] = sm_tag
                sample_to_bam[sm_tag] = bam
                print(f"  BAM: {bam} → Sample: {sm_tag}", file=sys.stderr)
                break
    except Exception as e:
        print(f"  ERROR processing {bam}: {e}", file=sys.stderr)
        continue

print(f"  Found {len(bam_samples)} BAM files", file=sys.stderr)
print("", file=sys.stderr)

# Process ploidy map
print("Step 2: Processing ploidy map and matching to BAMs...", file=sys.stderr)
print("------------------------------------------------------", file=sys.stderr)

matched_count = 0
unmatched_count = 0

try:
    with open('${ploidy_map}', 'r') as f:
        for line in f:
            # Skip comments and empty lines
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            # Parse the line (handle tabs or spaces)
            parts = line.split()
            if len(parts) < 2:
                continue
                
            sample_id = parts[0]
            try:
                ploidy_value = int(parts[1])
            except ValueError:
                print(f"  ERROR: Invalid ploidy value '{parts[1]}' for sample '{sample_id}'", file=sys.stderr)
                continue
            
            # Check if we have a BAM for this sample
            if sample_id in sample_to_bam:
                # Found matching BAM
                print(f"{sample_id}\\t{ploidy_value}")
                print(f"  ✔ Matched: {sample_id} (ploidy={ploidy_value}) → BAM: {sample_to_bam[sample_id]}", file=sys.stderr)
                matched_count += 1
            else:
                # No matching BAM found
                print(f"  ✗ Unmatched: {sample_id} (ploidy={ploidy_value}) - No corresponding BAM found", file=sys.stderr)
                
                # Try to find similar sample names
                print("    Checking for similar sample names...", file=sys.stderr)
                for bam_sample in sample_to_bam.keys():
                    if sample_id in bam_sample or bam_sample in sample_id:
                        print(f"    → Possible match: '{bam_sample}' in {sample_to_bam[bam_sample]}", file=sys.stderr)
                unmatched_count += 1

except Exception as e:
    print(f"ERROR reading ploidy map: {e}", file=sys.stderr)
    sys.exit(1)

print("", file=sys.stderr)
print("Step 3: Summary", file=sys.stderr)
print("---------------", file=sys.stderr)
print(f"  Ploidy map entries processed: {matched_count + unmatched_count}", file=sys.stderr)
print(f"  Successfully matched: {matched_count}", file=sys.stderr)
print(f"  Unmatched entries: {unmatched_count}", file=sys.stderr)
print(f"  Total BAMs available: {len(bam_samples)}", file=sys.stderr)

# Check for BAMs without ploidy information
print("", file=sys.stderr)
print("Step 4: Checking for BAMs without ploidy information...", file=sys.stderr)
print("--------------------------------------------------------", file=sys.stderr)

# Get list of matched samples
matched_samples = set()
with open('${ploidy_map}', 'r') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split()
        if len(parts) >= 2 and parts[0] in sample_to_bam:
            matched_samples.add(parts[0])

bams_without_ploidy = 0
for bam, sample in bam_samples.items():
    if sample not in matched_samples:
        print(f"  WARNING: BAM '{bam}' (sample: {sample}) has no ploidy information", file=sys.stderr)
        bams_without_ploidy += 1

if bams_without_ploidy > 0:
    print(f"  → {bams_without_ploidy} BAM(s) will use default ploidy or global setting", file=sys.stderr)
else:
    print("  → All BAMs have ploidy information", file=sys.stderr)

print("", file=sys.stderr)
print("================================================================", file=sys.stderr)

if matched_count == 0:
    print("ERROR: CNV map is empty! No samples could be matched.", file=sys.stderr)
    print("Please check that sample names in your ploidy map match the SM tags in your BAM files.", file=sys.stderr)
else:
    print(f"SUCCESS: CNV map created with {matched_count} sample(s)", file=sys.stderr)
PYTHON_CNV
        
        # Check if CNV map was created successfully
        if [ -s cnv_map.txt ]; then
            echo "CNV map created successfully:"
            echo "-----------------------"
            cat cnv_map.txt | while IFS=\$'\\t' read -r sample ploidy; do
                echo "  \$sample: ploidy=\$ploidy"
            done
            echo ""
            CNV_MAP_OPTION="--cnv-map cnv_map.txt"
        else
            echo "WARNING: CNV map is empty or was not created"
            echo "FreeBayes will use default ploidy for all samples"
            CNV_MAP_OPTION=""
        fi
    else
        echo "No ploidy map provided - using global ploidy from config or FreeBayes default"
        CNV_MAP_OPTION=""
    fi
    
    # Build freebayes command
    FREEBAYES_CMD="freebayes --fasta-reference ${reference} --targets chunk.bed"
    
    # Add CNV map if available
    if [ -n "\$CNV_MAP_OPTION" ]; then
        FREEBAYES_CMD="\$FREEBAYES_CMD \$CNV_MAP_OPTION"
    fi
    
    # Add BAM files
    for bam in *.bam; do
        if [ -f "\$bam" ]; then
            FREEBAYES_CMD="\$FREEBAYES_CMD --bam \$bam"
        fi
    done
    
    # Add config options if provided - COMPLETE PARAMETER EXTRACTION
    if [ "${has_config}" = "true" ]; then
        CONFIG_FILE=\$(ls *.json 2>/dev/null | head -n1)
        if [ -n "\$CONFIG_FILE" ]; then
            echo "Loading config from \$CONFIG_FILE"
            
            # Parse ALL config options using Python
            python3 << 'PYTHON_PARSE' > freebayes_params.sh
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
    
    # Extract sections based on new structure
    input_filters = config.get('input_filters', {})
    pop_model = config.get('population_model', {})
    allele_scope = config.get('allele_scope', {})
    priors = config.get('priors', {})
    genotype_lh = config.get('genotype_likelihoods', {})
    algorithmic = config.get('algorithmic_features', {})
    output_opts = config.get('output_options', {})
    indel_realign = config.get('indel_realignment', {})
    debug_opts = config.get('debugging', {})
    
    params = []
    
    # ========== INPUT FILTERS ==========
    if input_filters.get('use_duplicate_reads'):
        params.append("--use-duplicate-reads")
    if 'min_mapping_quality' in input_filters:
        params.append(f"--min-mapping-quality {input_filters['min_mapping_quality']}")
    if 'min_base_quality' in input_filters:
        params.append(f"--min-base-quality {input_filters['min_base_quality']}")
    if 'min_supporting_allele_qsum' in input_filters:
        params.append(f"--min-supporting-allele-qsum {input_filters['min_supporting_allele_qsum']}")
    if 'min_supporting_mapping_qsum' in input_filters:
        params.append(f"--min-supporting-mapping-qsum {input_filters['min_supporting_mapping_qsum']}")
    if 'mismatch_base_quality_threshold' in input_filters:
        params.append(f"--mismatch-base-quality-threshold {input_filters['mismatch_base_quality_threshold']}")
    if input_filters.get('read_mismatch_limit') is not None:
        params.append(f"--read-mismatch-limit {input_filters['read_mismatch_limit']}")
    if 'read_max_mismatch_fraction' in input_filters:
        params.append(f"--read-max-mismatch-fraction {input_filters['read_max_mismatch_fraction']}")
    if input_filters.get('read_snp_limit') is not None:
        params.append(f"--read-snp-limit {input_filters['read_snp_limit']}")
    if input_filters.get('read_indel_limit') is not None:
        params.append(f"--read-indel-limit {input_filters['read_indel_limit']}")
    if input_filters.get('standard_filters'):
        params.append("--standard-filters")
    if 'min_alternate_fraction' in input_filters:
        params.append(f"--min-alternate-fraction {input_filters['min_alternate_fraction']}")
    if 'min_alternate_count' in input_filters:
        params.append(f"--min-alternate-count {input_filters['min_alternate_count']}")
    if 'min_alternate_qsum' in input_filters:
        params.append(f"--min-alternate-qsum {input_filters['min_alternate_qsum']}")
    if 'min_alternate_total' in input_filters:
        params.append(f"--min-alternate-total {input_filters['min_alternate_total']}")
    if 'min_coverage' in input_filters:
        params.append(f"--min-coverage {input_filters['min_coverage']}")
    if input_filters.get('limit_coverage') is not None:
        params.append(f"--limit-coverage {input_filters['limit_coverage']}")
    if input_filters.get('skip_coverage') is not None:
        params.append(f"--skip-coverage {input_filters['skip_coverage']}")
    if input_filters.get('trim_complex_tail'):
        params.append("--trim-complex-tail")
    
    # ========== POPULATION MODEL ==========
    if 'theta' in pop_model:
        params.append(f"--theta {pop_model['theta']}")
    
    # Convert bash boolean string to Python boolean
    has_ploidy_map = '${has_ploidy_map}' == 'true'
    
    # Only add global ploidy if no ploidy map is provided
    if not has_ploidy_map and 'ploidy' in pop_model:
        params.append(f"--ploidy {pop_model['ploidy']}")
        print(f"echo 'Global ploidy set to: {pop_model['ploidy']}'")
    elif has_ploidy_map:
        print("echo 'Using per-sample ploidy from CNV map'")
    
    if pop_model.get('pooled_discrete'):
        params.append("--pooled-discrete")
        print("echo 'Using pooled-discrete mode'")
    if pop_model.get('pooled_continuous'):
        params.append("--pooled-continuous")
        print("echo 'Using pooled-continuous mode'")
    if pop_model.get('use_reference_allele'):
        params.append("--use-reference-allele")
    if pop_model.get('reference_quality'):
        params.append(f"--reference-quality {pop_model['reference_quality']}")
    
    # ========== ALLELE SCOPE ==========
    if 'use_best_n_alleles' in allele_scope and allele_scope['use_best_n_alleles'] != 0:
        params.append(f"--use-best-n-alleles {allele_scope['use_best_n_alleles']}")
    if 'max_complex_gap' in allele_scope:
        params.append(f"--max-complex-gap {allele_scope['max_complex_gap']}")
    if 'haplotype_length' in allele_scope:
        params.append(f"--haplotype-length {allele_scope['haplotype_length']}")
    if 'min_repeat_size' in allele_scope:
        params.append(f"--min-repeat-size {allele_scope['min_repeat_size']}")
    if 'min_repeat_entropy' in allele_scope:
        params.append(f"--min-repeat-entropy {allele_scope['min_repeat_entropy']}")
    if allele_scope.get('no_partial_observations'):
        params.append("--no-partial-observations")
    if allele_scope.get('throw_away_snp_obs'):
        params.append("--throw-away-snp-obs")
    if allele_scope.get('throw_away_indel_obs'):
        params.append("--throw-away-indel-obs")
    if allele_scope.get('throw_away_mnps_obs'):
        params.append("--throw-away-mnps-obs")
    if allele_scope.get('throw_away_complex_obs'):
        params.append("--throw-away-complex-obs")
    
    # ========== PRIORS ==========
    if priors.get('no_population_priors'):
        params.append("--no-population-priors")
    if priors.get('hwe_priors_off'):
        params.append("--hwe-priors-off")
    if priors.get('binomial_obs_priors_off'):
        params.append("--binomial-obs-priors-off")
    if priors.get('allele_balance_priors_off'):
        params.append("--allele-balance-priors-off")
    
    # ========== GENOTYPE LIKELIHOODS ==========
    if genotype_lh.get('observation_bias'):
        params.append(f"--observation-bias {genotype_lh['observation_bias']}")
    if genotype_lh.get('base_quality_cap') is not None:
        params.append(f"--base-quality-cap {genotype_lh['base_quality_cap']}")
    if 'prob_contamination' in genotype_lh:
        params.append(f"--prob-contamination {genotype_lh['prob_contamination']}")
    if genotype_lh.get('legacy_gls'):
        params.append("--legacy-gls")
    if genotype_lh.get('contamination_estimates'):
        params.append(f"--contamination-estimates {genotype_lh['contamination_estimates']}")
    
    # ========== ALGORITHMIC FEATURES ==========
    if algorithmic.get('report_genotype_likelihood_max'):
        params.append("--report-genotype-likelihood-max")
    if 'genotyping_max_iterations' in algorithmic:
        params.append(f"--genotyping-max-iterations {algorithmic['genotyping_max_iterations']}")
    if 'genotyping_max_banddepth' in algorithmic:
        params.append(f"--genotyping-max-banddepth {algorithmic['genotyping_max_banddepth']}")
    if 'posterior_integration_limits' in algorithmic:
        params.append(f"--posterior-integration-limits {algorithmic['posterior_integration_limits']}")
    if algorithmic.get('exclude_unobserved_genotypes'):
        params.append("--exclude-unobserved-genotypes")
    if algorithmic.get('genotype_variant_threshold') is not None:
        params.append(f"--genotype-variant-threshold {algorithmic['genotype_variant_threshold']}")
    if algorithmic.get('use_mapping_quality'):
        params.append("--use-mapping-quality")
    if algorithmic.get('harmonic_indel_quality'):
        params.append("--harmonic-indel-quality")
    if 'read_dependence_factor' in algorithmic:
        params.append(f"--read-dependence-factor {algorithmic['read_dependence_factor']}")
    if algorithmic.get('genotype_qualities'):
        params.append("--genotype-qualities")
    
    # ========== OUTPUT OPTIONS ==========
    if output_opts.get('gvcf'):
        params.append("--gvcf")
    if output_opts.get('gvcf_chunk') is not None:
        params.append(f"--gvcf-chunk {output_opts['gvcf_chunk']}")
    if output_opts.get('gvcf_dont_use_chunk'):
        params.append("--gvcf-dont-use-chunk true")
    if output_opts.get('variant_input'):
        params.append(f"--variant-input {output_opts['variant_input']}")
    if output_opts.get('only_use_input_alleles'):
        params.append("--only-use-input-alleles")
    if output_opts.get('haplotype_basis_alleles'):
        params.append(f"--haplotype-basis-alleles {output_opts['haplotype_basis_alleles']}")
    if output_opts.get('report_all_haplotype_alleles'):
        params.append("--report-all-haplotype-alleles")
    if output_opts.get('report_monomorphic'):
        params.append("--report-monomorphic")
    if 'pvar' in output_opts and output_opts['pvar'] > 0:
        params.append(f"--pvar {output_opts['pvar']}")
    if output_opts.get('strict_vcf'):
        params.append("--strict-vcf")
    
    # ========== INDEL REALIGNMENT ==========
    if indel_realign.get('dont_left_align_indels'):
        params.append("--dont-left-align-indels")
    
    # ========== DEBUGGING ==========
    if debug_opts.get('debug'):
        params.append("--debug")
    if debug_opts.get('dd'):
        params.append("-dd")
    
    # Output the parameters
    print(f"FREEBAYES_PARAMS='{' '.join(params)}'")
    
    # Print summary of what was configured
    param_count = len(params)
    print(f"echo 'Loaded {param_count} parameters from config file'")
    
    # List key parameters for confirmation
    if pop_model.get('pooled_discrete'):
        print("echo '  - Pooled discrete mode: ENABLED'")
    if pop_model.get('pooled_continuous'):
        print("echo '  - Pooled continuous mode: ENABLED'")
    if 'min_alternate_fraction' in input_filters:
        print(f"echo '  - Min alternate fraction: {input_filters['min_alternate_fraction']}'")
    if 'min_coverage' in input_filters:
        print(f"echo '  - Min coverage: {input_filters['min_coverage']}'")
    if 'limit_coverage' in input_filters:
        print(f"echo '  - Limit coverage: {input_filters['limit_coverage']}'")
    
except Exception as e:
    print(f"# Error parsing config: {e}", file=sys.stderr)
    print("FREEBAYES_PARAMS=''")
    print("echo 'WARNING: Failed to parse config file, using defaults'")
    print(f"echo 'Error details: {e}'")
PYTHON_PARSE
            
            # Source the parameters
            source freebayes_params.sh
            
            # Add parsed parameters to command
            if [ -n "\$FREEBAYES_PARAMS" ]; then
                FREEBAYES_CMD="\$FREEBAYES_CMD \$FREEBAYES_PARAMS"
                echo "Using freebayes parameters from config"
            fi
        fi
    else
        echo "No config file provided - using FreeBayes defaults"
    fi
    
    # Add output
    FREEBAYES_CMD="\$FREEBAYES_CMD --vcf chunk_${chunk_id}.vcf"
    
    # Run freebayes
    echo "Full command: \$FREEBAYES_CMD"
    \$FREEBAYES_CMD
    
    # Check output and compress
    if [ -s "chunk_${chunk_id}.vcf" ]; then
        echo "Freebayes completed successfully"
        
        # Count variants
        VARIANT_COUNT=\$(grep -v '^#' chunk_${chunk_id}.vcf | wc -l)
        echo "Found \$VARIANT_COUNT variants"
        
        # Show samples
        if grep -q '^#CHROM' chunk_${chunk_id}.vcf; then
            SAMPLE_COUNT=\$(grep '^#CHROM' chunk_${chunk_id}.vcf | cut -f10- | wc -w)
            echo "Samples (\$SAMPLE_COUNT): \$(grep '^#CHROM' chunk_${chunk_id}.vcf | cut -f10- | tr '\\t' ' ')"
        fi
        
        # Compress and index
        bgzip chunk_${chunk_id}.vcf
        tabix -p vcf chunk_${chunk_id}.vcf.gz
        
    else
        echo "No variants found, creating empty VCF"
        
        # Create minimal header
        echo '##fileformat=VCFv4.2' > chunk_${chunk_id}.vcf
        echo '##reference=${reference}' >> chunk_${chunk_id}.vcf
        echo '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' >> chunk_${chunk_id}.vcf
        
        bgzip chunk_${chunk_id}.vcf
        tabix -p vcf chunk_${chunk_id}.vcf.gz
    fi
    
    echo "Chunk ${chunk_id} complete"
    """
}