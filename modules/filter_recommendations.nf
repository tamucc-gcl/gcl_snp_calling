// modules/filter_recommendations.nf - Analyze QC metrics and recommend filters

process FILTER_RECOMMENDATIONS {
    tag "filter_recommendations"
    publishDir "${params.output_dir}/qc", mode: 'copy'
    
    input:
    path variant_stats
    path locus_metrics
    path sample_metrics
    path depth_dist
    path qual_dist
    path vcf
    
    output:
    path "filtering_recommendations.txt", emit: recommendations
    path "filter_commands.sh", emit: commands
    path "filter_thresholds.json", emit: thresholds
    path "pre_post_filter_comparison.tsv", emit: comparison
    
    script:
    """
    echo "Analyzing QC metrics to generate filtering recommendations..."
    
    # Create Python script to analyze distributions and suggest thresholds
    cat > analyze_thresholds.py << 'PYTHON_SCRIPT'
import json
import numpy as np
import pandas as pd
from collections import defaultdict

# Read QC data
locus_df = pd.read_csv('${locus_metrics}', sep='\\t')
depth_df = pd.read_csv('${depth_dist}', sep='\\t')
qual_df = pd.read_csv('${qual_dist}', sep='\\t')

# Calculate distribution statistics
depth_values = []
qual_values = []

# Expand binned distributions to approximate values
depth_bins = {
    "0-5": 2.5, "5-10": 7.5, "10-20": 15, "20-30": 25, 
    "30-50": 40, "50-100": 75, "100-200": 150, "200+": 250
}

qual_bins = {
    "0-10": 5, "10-20": 15, "20-30": 25, "30-50": 40,
    "50-100": 75, "100-200": 150, "200-500": 350, "500+": 750
}

# Approximate depth distribution
for _, row in depth_df.iterrows():
    if row['depth_bin'] in depth_bins:
        depth_values.extend([depth_bins[row['depth_bin']]] * int(row['count']))

# Approximate quality distribution
for _, row in qual_df.iterrows():
    if row['quality_bin'] in qual_bins:
        qual_values.extend([qual_bins[row['quality_bin']]] * int(row['count']))

# Calculate percentiles
depth_percentiles = {
    'p1': np.percentile(depth_values, 1) if depth_values else 0,
    'p5': np.percentile(depth_values, 5) if depth_values else 0,
    'p10': np.percentile(depth_values, 10) if depth_values else 0,
    'median': np.percentile(depth_values, 50) if depth_values else 0,
    'p90': np.percentile(depth_values, 90) if depth_values else 0,
    'p95': np.percentile(depth_values, 95) if depth_values else 0,
    'p99': np.percentile(depth_values, 99) if depth_values else 0
}

qual_percentiles = {
    'p1': np.percentile(qual_values, 1) if qual_values else 0,
    'p5': np.percentile(qual_values, 5) if qual_values else 0,
    'p10': np.percentile(qual_values, 10) if qual_values else 0,
    'median': np.percentile(qual_values, 50) if qual_values else 0,
    'p90': np.percentile(qual_values, 90) if qual_values else 0,
    'p95': np.percentile(qual_values, 95) if qual_values else 0,
    'p99': np.percentile(qual_values, 99) if qual_values else 0
}

# Analyze missingness
miss_vals = locus_df['missingness_rate'].dropna()
miss_percentiles = {
    'p50': miss_vals.quantile(0.5) if len(miss_vals) > 0 else 0,
    'p90': miss_vals.quantile(0.9) if len(miss_vals) > 0 else 0,
    'p95': miss_vals.quantile(0.95) if len(miss_vals) > 0 else 0,
    'p99': miss_vals.quantile(0.99) if len(miss_vals) > 0 else 0
}

# Determine recommended thresholds
thresholds = {
    'min_depth': max(5, int(depth_percentiles['p5'])),
    'max_depth': min(int(depth_percentiles['p99']), int(depth_percentiles['median'] * 3)),
    'min_quality': max(20, int(qual_percentiles['p5'])),
    'max_missingness': min(0.2, round(miss_percentiles['p95'], 2)),
    'min_maf': 0.01  # Default for typical use case
}

# Adjust for pooled samples if detected (high ploidy)
if 'ploidy' in locus_df.columns and locus_df['ploidy'].max() > 2:
    thresholds['min_maf'] = 0.001  # Lower MAF for pools
    thresholds['min_depth'] = max(10, thresholds['min_depth'])  # Higher depth for pools

# Count variants that would pass filters
total_vars = len(locus_df)
passing_vars = locus_df[
    (locus_df['depth'] >= thresholds['min_depth']) &
    (locus_df['depth'] <= thresholds['max_depth']) &
    (locus_df['quality'] >= thresholds['min_quality']) &
    (locus_df['missingness_rate'] <= thresholds['max_missingness'])
]
n_passing = len(passing_vars)
pct_passing = (n_passing / total_vars * 100) if total_vars > 0 else 0

# Save thresholds
with open('filter_thresholds.json', 'w') as f:
    json.dump({
        'recommended_thresholds': thresholds,
        'depth_percentiles': depth_percentiles,
        'quality_percentiles': qual_percentiles,
        'missingness_percentiles': miss_percentiles,
        'filtering_impact': {
            'total_variants': total_vars,
            'passing_variants': n_passing,
            'percent_passing': round(pct_passing, 2)
        }
    }, f, indent=2)

print(f"Recommended thresholds saved to filter_thresholds.json")
print(f"Would retain {n_passing}/{total_vars} variants ({pct_passing:.1f}%)")
PYTHON_SCRIPT
    
    # Run analysis
    python3 analyze_thresholds.py
    
    # Read the generated thresholds
    MIN_DEPTH=\$(python3 -c "import json; print(json.load(open('filter_thresholds.json'))['recommended_thresholds']['min_depth'])")
    MAX_DEPTH=\$(python3 -c "import json; print(json.load(open('filter_thresholds.json'))['recommended_thresholds']['max_depth'])")
    MIN_QUAL=\$(python3 -c "import json; print(json.load(open('filter_thresholds.json'))['recommended_thresholds']['min_quality'])")
    MAX_MISS=\$(python3 -c "import json; print(json.load(open('filter_thresholds.json'))['recommended_thresholds']['max_missingness'])")
    MIN_MAF=\$(python3 -c "import json; print(json.load(open('filter_thresholds.json'))['recommended_thresholds']['min_maf'])")
    
    # Generate filtering recommendations report
    cat > filtering_recommendations.txt << RECOMMENDATIONS
================================================================================
                        VARIANT FILTERING RECOMMENDATIONS
================================================================================

Based on the QC analysis of your variant dataset, we recommend the following
filtering strategy to balance quality and retention of informative variants.

--------------------------------------------------------------------------------
RECOMMENDED THRESHOLDS
--------------------------------------------------------------------------------

1. DEPTH FILTERING
   - Minimum depth: \${MIN_DEPTH}X
   - Maximum depth: \${MAX_DEPTH}X
   - Rationale: Removes low-confidence calls and potential mapping artifacts
   
2. QUALITY FILTERING
   - Minimum quality score (QUAL): \${MIN_QUAL}
   - Rationale: Ensures high-confidence variant calls
   
3. MISSINGNESS FILTERING
   - Maximum missingness per variant: \${MAX_MISS} ($(echo "\${MAX_MISS} * 100" | bc | cut -d. -f1)%)
   - Rationale: Removes poorly called variants across samples
   
4. ALLELE FREQUENCY FILTERING
   - Minimum minor allele frequency: \${MIN_MAF}
   - Rationale: Removes very rare variants that may be errors
   
5. SAMPLE-LEVEL FILTERING
   - Remove samples with >10% missing genotypes
   - Remove samples with extreme heterozygosity (>3 SD from mean)
   - Rationale: Identifies potentially contaminated or poor-quality samples

--------------------------------------------------------------------------------
DISTRIBUTION STATISTICS
--------------------------------------------------------------------------------

\$(python3 -c "
import json
data = json.load(open('filter_thresholds.json'))
print('DEPTH DISTRIBUTION:')
for k, v in data['depth_percentiles'].items():
    print(f'  {k:10s}: {v:>8.1f}X')
print()
print('QUALITY DISTRIBUTION:')
for k, v in data['quality_percentiles'].items():
    print(f'  {k:10s}: {v:>8.1f}')
print()
print('MISSINGNESS DISTRIBUTION:')
for k, v in data['missingness_percentiles'].items():
    print(f'  {k:10s}: {v:>8.3f}')
")

--------------------------------------------------------------------------------
EXPECTED FILTERING IMPACT
--------------------------------------------------------------------------------

\$(python3 -c "
import json
data = json.load(open('filter_thresholds.json'))
impact = data['filtering_impact']
print(f'Total variants before filtering:     {impact[\"total_variants\"]:,}')
print(f'Expected variants after filtering:   {impact[\"passing_variants\"]:,}')
print(f'Retention rate:                      {impact[\"percent_passing\"]:.1f}%')
")

--------------------------------------------------------------------------------
ADDITIONAL CONSIDERATIONS
--------------------------------------------------------------------------------

1. Hardy-WEINBERG EQUILIBRIUM
   - Consider filtering variants with HWE p-value < 1e-6 in control samples
   - Be cautious with HWE filtering in case-control studies or pooled samples

2. TECHNICAL ARTIFACTS
   - Review variants in repetitive regions (use RepeatMasker annotations)
   - Check for strand bias and mapping quality bias
   - Consider removing variants in segmental duplications

3. POPULATION-SPECIFIC FILTERING
   - For pooled samples: adjust MAF threshold based on pool size
   - For family studies: check Mendelian inheritance patterns
   - For case-control: consider differential missingness between groups

4. POST-FILTERING VALIDATION
   - Verify Ti/Tv ratio remains >2.0 for WGS (>3.0 for WES)
   - Check that filtering doesn't create sample stratification
   - Validate a subset of filtered variants by alternative methods

--------------------------------------------------------------------------------
NEXT STEPS
--------------------------------------------------------------------------------

1. Apply the recommended filters using the provided shell script
2. Review the filtered dataset statistics
3. Perform functional annotation with SnpEff or VEP
4. Consider additional filtering based on functional impact
5. Validate filtering strategy with known variants or external datasets

================================================================================
RECOMMENDATIONS

    # Generate filter commands script
    cat > filter_commands.sh << 'FILTER_SCRIPT'
#!/bin/bash

# Variant Filtering Commands
# Generated based on QC analysis

# Set variables (modify as needed)
INPUT_VCF="${vcf}"
OUTPUT_DIR="filtered_variants"
MIN_DEPTH=\${MIN_DEPTH}
MAX_DEPTH=\${MAX_DEPTH}
MIN_QUAL=\${MIN_QUAL}
MAX_MISS=\${MAX_MISS}
MIN_MAF=\${MIN_MAF}

# Create output directory
mkdir -p \$OUTPUT_DIR

echo "Starting variant filtering pipeline..."
echo "Input VCF: \$INPUT_VCF"

# Step 1: Basic quality filtering
echo "Step 1: Applying basic quality filters..."
bcftools filter \\
    -e "QUAL < \$MIN_QUAL || INFO/DP < \$MIN_DEPTH || INFO/DP > \$MAX_DEPTH" \\
    -o \$OUTPUT_DIR/step1_quality_filtered.vcf.gz \\
    -O z \\
    \$INPUT_VCF

bcftools index \$OUTPUT_DIR/step1_quality_filtered.vcf.gz

# Step 2: Missingness filtering
echo "Step 2: Filtering by missingness..."
bcftools filter \\
    -e "F_MISSING > \$MAX_MISS" \\
    -o \$OUTPUT_DIR/step2_missingness_filtered.vcf.gz \\
    -O z \\
    \$OUTPUT_DIR/step1_quality_filtered.vcf.gz

bcftools index \$OUTPUT_DIR/step2_missingness_filtered.vcf.gz

# Step 3: Allele frequency filtering
echo "Step 3: Filtering by minor allele frequency..."
bcftools filter \\
    -e "MAF < \$MIN_MAF" \\
    -o \$OUTPUT_DIR/step3_maf_filtered.vcf.gz \\
    -O z \\
    \$OUTPUT_DIR/step2_missingness_filtered.vcf.gz

bcftools index \$OUTPUT_DIR/step3_maf_filtered.vcf.gz

# Step 4: Hardy-Weinberg filtering (optional, commented out)
# echo "Step 4: Filtering by Hardy-Weinberg equilibrium..."
# bcftools filter \\
#     -e "INFO/HWE < 1e-6" \\
#     -o \$OUTPUT_DIR/step4_hwe_filtered.vcf.gz \\
#     -O z \\
#     \$OUTPUT_DIR/step3_maf_filtered.vcf.gz

# Final filtered file
cp \$OUTPUT_DIR/step3_maf_filtered.vcf.gz \$OUTPUT_DIR/final_filtered.vcf.gz
cp \$OUTPUT_DIR/step3_maf_filtered.vcf.gz.csi \$OUTPUT_DIR/final_filtered.vcf.gz.csi

# Generate filtering statistics
echo ""
echo "Filtering Statistics:"
echo "===================="
echo -n "Original variants: "
bcftools view -H \$INPUT_VCF | wc -l

echo -n "After quality filtering: "
bcftools view -H \$OUTPUT_DIR/step1_quality_filtered.vcf.gz | wc -l

echo -n "After missingness filtering: "
bcftools view -H \$OUTPUT_DIR/step2_missingness_filtered.vcf.gz | wc -l

echo -n "After MAF filtering: "
bcftools view -H \$OUTPUT_DIR/final_filtered.vcf.gz | wc -l

# Calculate Ti/Tv ratio
echo ""
echo "Ti/Tv ratio check:"
echo -n "Before filtering: "
bcftools stats \$INPUT_VCF | grep "TSTV" | tail -1 | cut -f5

echo -n "After filtering: "
bcftools stats \$OUTPUT_DIR/final_filtered.vcf.gz | grep "TSTV" | tail -1 | cut -f5

echo ""
echo "Filtering complete!"
echo "Final filtered VCF: \$OUTPUT_DIR/final_filtered.vcf.gz"
FILTER_SCRIPT

    chmod +x filter_commands.sh
    
    # Generate pre/post filtering comparison
    echo "stage\ttotal_variants\tsnps\tindels\tmean_depth\tmean_quality\tti_tv_ratio" > pre_post_filter_comparison.tsv
    
    # Pre-filter stats
    pre_total=\$(bcftools view -H ${vcf} | wc -l)
    pre_snps=\$(bcftools view -H -v snps ${vcf} | wc -l)
    pre_indels=\$(bcftools view -H -v indels ${vcf} | wc -l)
    pre_titv=\$(bcftools stats ${vcf} | grep "TSTV" | tail -1 | cut -f5)
    
    echo -e "pre_filter\\t\$pre_total\\t\$pre_snps\\t\$pre_indels\\tNA\\tNA\\t\$pre_titv" >> pre_post_filter_comparison.tsv
    
    # Simulate post-filter stats (would be calculated after actual filtering)
    post_total=\$(python3 -c "import json; print(json.load(open('filter_thresholds.json'))['filtering_impact']['passing_variants'])")
    post_snps=\$(echo "\$post_total * \$pre_snps / \$pre_total" | bc)
    post_indels=\$(echo "\$post_total * \$pre_indels / \$pre_total" | bc)
    
    echo -e "post_filter_expected\\t\$post_total\\t\$post_snps\\t\$post_indels\\tNA\\tNA\\tNA" >> pre_post_filter_comparison.tsv
    
    echo "Filter recommendations generated successfully!"
    """
}