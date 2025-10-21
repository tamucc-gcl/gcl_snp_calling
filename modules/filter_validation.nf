// modules/filter_validation.nf - Test and validate filtering strategies

process FILTER_VALIDATION {
    tag "filter_validation"
    publishDir "${params.output_dir}/qc/filter_validation", mode: 'copy'
    
    input:
    path vcf
    path vcf_index
    path filter_thresholds
    
    output:
    path "filter_validation_report.tsv", emit: report
    path "filter_impact_summary.json", emit: summary
    path "titv_by_filter.tsv", emit: titv
    path "maf_spectrum_comparison.tsv", emit: maf_spectrum
    path "filter_validation_plots.pdf", emit: plots
    
    script:
    """
    echo "Validating filter effectiveness..."
    
    # Read filter thresholds
    MIN_DEPTH=\$(python3 -c "import json; print(json.load(open('${filter_thresholds}'))['recommended_thresholds']['min_depth'])")
    MAX_DEPTH=\$(python3 -c "import json; print(json.load(open('${filter_thresholds}'))['recommended_thresholds']['max_depth'])")
    MIN_QUAL=\$(python3 -c "import json; print(json.load(open('${filter_thresholds}'))['recommended_thresholds']['min_quality'])")
    MAX_MISS=\$(python3 -c "import json; print(json.load(open('${filter_thresholds}'))['recommended_thresholds']['max_missingness'])")
    MIN_MAF=\$(python3 -c "import json; print(json.load(open('${filter_thresholds}'))['recommended_thresholds']['min_maf'])")
    
    # Test different filtering strategies
    echo "filter_name\ttotal_variants\tsnps\tindels\tti_tv_ratio\tmean_depth\tmean_qual\tretention_rate" > filter_validation_report.tsv
    
    # 1. No filtering (baseline)
    echo "Analyzing unfiltered data..."
    total_orig=\$(bcftools view -H ${vcf} | wc -l)
    snps_orig=\$(bcftools view -H -v snps ${vcf} | wc -l)
    indels_orig=\$(bcftools view -H -v indels ${vcf} | wc -l)
    titv_orig=\$(bcftools stats ${vcf} | grep "TSTV" | tail -1 | cut -f5)
    mean_depth_orig=\$(bcftools query -f '%INFO/DP\\n' ${vcf} | awk '{sum+=\$1; n++} END {if(n>0) print sum/n; else print 0}')
    mean_qual_orig=\$(bcftools query -f '%QUAL\\n' ${vcf} | awk '{sum+=\$1; n++} END {if(n>0) print sum/n; else print 0}')
    
    echo -e "no_filter\\t\$total_orig\\t\$snps_orig\\t\$indels_orig\\t\$titv_orig\\t\$mean_depth_orig\\t\$mean_qual_orig\\t100.0" >> filter_validation_report.tsv
    
    # 2. Quality filter only
    echo "Testing quality filter only..."
    bcftools filter -e "QUAL < \$MIN_QUAL" ${vcf} -o temp_qual.vcf
    total_qual=\$(grep -v "^#" temp_qual.vcf | wc -l)
    snps_qual=\$(bcftools view -H -v snps temp_qual.vcf | wc -l)
    indels_qual=\$(bcftools view -H -v indels temp_qual.vcf | wc -l)
    titv_qual=\$(bcftools stats temp_qual.vcf | grep "TSTV" | tail -1 | cut -f5)
    retention_qual=\$(echo "scale=2; \$total_qual * 100 / \$total_orig" | bc)
    
    echo -e "qual_only\\t\$total_qual\\t\$snps_qual\\t\$indels_qual\\t\$titv_qual\\tNA\\tNA\\t\$retention_qual" >> filter_validation_report.tsv
    
    # 3. Depth filter only
    echo "Testing depth filter only..."
    bcftools filter -e "INFO/DP < \$MIN_DEPTH || INFO/DP > \$MAX_DEPTH" ${vcf} -o temp_depth.vcf
    total_depth=\$(grep -v "^#" temp_depth.vcf | wc -l)
    snps_depth=\$(bcftools view -H -v snps temp_depth.vcf | wc -l)
    indels_depth=\$(bcftools view -H -v indels temp_depth.vcf | wc -l)
    titv_depth=\$(bcftools stats temp_depth.vcf | grep "TSTV" | tail -1 | cut -f5)
    retention_depth=\$(echo "scale=2; \$total_depth * 100 / \$total_orig" | bc)
    
    echo -e "depth_only\\t\$total_depth\\t\$snps_depth\\t\$indels_depth\\t\$titv_depth\\tNA\\tNA\\t\$retention_depth" >> filter_validation_report.tsv
    
    # 4. Combined basic filters (qual + depth)
    echo "Testing combined basic filters..."
    bcftools filter -e "QUAL < \$MIN_QUAL || INFO/DP < \$MIN_DEPTH || INFO/DP > \$MAX_DEPTH" ${vcf} -o temp_basic.vcf
    total_basic=\$(grep -v "^#" temp_basic.vcf | wc -l)
    snps_basic=\$(bcftools view -H -v snps temp_basic.vcf | wc -l)
    indels_basic=\$(bcftools view -H -v indels temp_basic.vcf | wc -l)
    titv_basic=\$(bcftools stats temp_basic.vcf | grep "TSTV" | tail -1 | cut -f5)
    retention_basic=\$(echo "scale=2; \$total_basic * 100 / \$total_orig" | bc)
    
    echo -e "basic_filters\\t\$total_basic\\t\$snps_basic\\t\$indels_basic\\t\$titv_basic\\tNA\\tNA\\t\$retention_basic" >> filter_validation_report.tsv
    
    # 5. Stringent filters (all recommended)
    echo "Testing all recommended filters..."
    # Note: Some filters may need to be adjusted based on available INFO fields
    bcftools filter -e "QUAL < \$MIN_QUAL || INFO/DP < \$MIN_DEPTH || INFO/DP > \$MAX_DEPTH" ${vcf} | \
    bcftools filter -e "F_MISSING > \$MAX_MISS" -o temp_stringent.vcf
    
    total_stringent=\$(grep -v "^#" temp_stringent.vcf | wc -l)
    snps_stringent=\$(bcftools view -H -v snps temp_stringent.vcf | wc -l)
    indels_stringent=\$(bcftools view -H -v indels temp_stringent.vcf | wc -l)
    titv_stringent=\$(bcftools stats temp_stringent.vcf | grep "TSTV" | tail -1 | cut -f5)
    retention_stringent=\$(echo "scale=2; \$total_stringent * 100 / \$total_orig" | bc)
    
    echo -e "stringent\\t\$total_stringent\\t\$snps_stringent\\t\$indels_stringent\\t\$titv_stringent\\tNA\\tNA\\t\$retention_stringent" >> filter_validation_report.tsv
    
    # Calculate Ti/Tv ratio at different quality thresholds
    echo "quality_threshold\tti_tv_ratio\tvariant_count" > titv_by_filter.tsv
    
    for qual_thresh in 0 10 20 30 50 100; do
        echo "Testing quality threshold: \$qual_thresh"
        bcftools filter -e "QUAL < \$qual_thresh" ${vcf} -o temp_q\${qual_thresh}.vcf
        
        titv=\$(bcftools stats temp_q\${qual_thresh}.vcf | grep "TSTV" | tail -1 | cut -f5)
        count=\$(grep -v "^#" temp_q\${qual_thresh}.vcf | wc -l)
        
        echo -e "\$qual_thresh\\t\$titv\\t\$count" >> titv_by_filter.tsv
        rm temp_q\${qual_thresh}.vcf
    done
    
    # Compare MAF spectrum before and after filtering
    echo "maf_bin\tunfiltered_count\tfiltered_count\tfilter_type" > maf_spectrum_comparison.tsv
    
    # Extract MAF for unfiltered
    bcftools query -f '%AF\\n' ${vcf} | \
    awk '{
        if(\$1 < 0.01) bin="0-0.01"
        else if(\$1 < 0.05) bin="0.01-0.05"
        else if(\$1 < 0.1) bin="0.05-0.1"
        else if(\$1 < 0.25) bin="0.1-0.25"
        else if(\$1 < 0.5) bin="0.25-0.5"
        else bin="0.5+"
        count[bin]++
    }
    END {
        for(b in count) print b"\t"count[b]
    }' | sort > unfiltered_maf.txt
    
    # Extract MAF for filtered
    bcftools query -f '%AF\\n' temp_stringent.vcf | \
    awk '{
        if(\$1 < 0.01) bin="0-0.01"
        else if(\$1 < 0.05) bin="0.01-0.05"
        else if(\$1 < 0.1) bin="0.05-0.1"
        else if(\$1 < 0.25) bin="0.1-0.25"
        else if(\$1 < 0.5) bin="0.25-0.5"
        else bin="0.5+"
        count[bin]++
    }
    END {
        for(b in count) print b"\t"count[b]
    }' | sort > filtered_maf.txt
    
    # Combine MAF distributions
    python3 << 'PYTHON_MAF'
import pandas as pd

bins = ["0-0.01", "0.01-0.05", "0.05-0.1", "0.1-0.25", "0.25-0.5", "0.5+"]

# Read unfiltered
unfiltered = {}
try:
    with open('unfiltered_maf.txt', 'r') as f:
        for line in f:
            parts = line.strip().split('\\t')
            if len(parts) == 2:
                unfiltered[parts[0]] = int(parts[1])
except:
    pass

# Read filtered
filtered = {}
try:
    with open('filtered_maf.txt', 'r') as f:
        for line in f:
            parts = line.strip().split('\\t')
            if len(parts) == 2:
                filtered[parts[0]] = int(parts[1])
except:
    pass

# Write comparison
with open('maf_spectrum_comparison.tsv', 'w') as f:
    f.write("maf_bin\\tunfiltered_count\\tfiltered_count\\n")
    for bin in bins:
        unf_count = unfiltered.get(bin, 0)
        filt_count = filtered.get(bin, 0)
        f.write(f"{bin}\\t{unf_count}\\t{filt_count}\\n")
PYTHON_MAF
    
    # Generate impact summary JSON
    python3 << 'PYTHON_SUMMARY'
import json

# Read validation report
validation_data = []
with open('filter_validation_report.tsv', 'r') as f:
    header = f.readline().strip().split('\\t')
    for line in f:
        parts = line.strip().split('\\t')
        validation_data.append(dict(zip(header, parts)))

# Calculate improvements
baseline = validation_data[0]  # no_filter
best_titv = max(float(d['ti_tv_ratio']) for d in validation_data if d['ti_tv_ratio'] != 'NA')
best_retention = max(float(d['retention_rate']) for d in validation_data[1:])

# Find optimal filter
optimal = None
for d in validation_data:
    if d['filter_name'] == 'stringent':
        optimal = d
        break

summary = {
    'baseline': {
        'total_variants': int(baseline['total_variants']),
        'ti_tv_ratio': float(baseline['ti_tv_ratio'])
    },
    'improvements': {
        'best_ti_tv_ratio': best_titv,
        'ti_tv_improvement': round(best_titv - float(baseline['ti_tv_ratio']), 3)
    },
    'optimal_filter': {
        'name': optimal['filter_name'] if optimal else 'basic_filters',
        'variants_retained': int(optimal['total_variants']) if optimal else 0,
        'retention_rate': float(optimal['retention_rate']) if optimal else 0,
        'ti_tv_ratio': float(optimal['ti_tv_ratio']) if optimal and optimal['ti_tv_ratio'] != 'NA' else 0
    },
    'recommendations': [
        'Quality filtering shows greatest Ti/Tv improvement' if best_titv > float(baseline['ti_tv_ratio']) + 0.2 else 'Minimal Ti/Tv improvement from filtering',
        f"Recommended filters retain {optimal['retention_rate']}% of variants" if optimal else "Apply basic filters",
        'Consider less stringent filtering if retention < 50%' if optimal and float(optimal['retention_rate']) < 50 else 'Good retention rate'
    ]
}

with open('filter_impact_summary.json', 'w') as f:
    json.dump(summary, f, indent=2)

print("Filter validation summary generated")
PYTHON_SUMMARY
    
    # Generate validation plots
    cat > plot_validation.R << 'R_SCRIPT'
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(scales)

# Read data
validation_df <- read_tsv("filter_validation_report.tsv", show_col_types = FALSE)
titv_df <- read_tsv("titv_by_filter.tsv", show_col_types = FALSE)
maf_df <- read_tsv("maf_spectrum_comparison.tsv", show_col_types = FALSE)

# Set theme
theme_set(theme_bw(base_size = 12))

# Plot 1: Retention vs Ti/Tv trade-off
p1 <- ggplot(validation_df %>% filter(ti_tv_ratio != "NA"), 
             aes(x = retention_rate, y = as.numeric(ti_tv_ratio))) +
    geom_point(size = 4, color = "#2E7D32") +
    geom_text(aes(label = filter_name), hjust = -0.1, vjust = 0.5, size = 3) +
    geom_hline(yintercept = 2.0, linetype = "dashed", color = "red", alpha = 0.5) +
    scale_x_continuous(limits = c(0, 105)) +
    labs(title = "Filter Trade-off: Retention vs Quality",
         subtitle = "Higher Ti/Tv indicates better quality",
         x = "Variant Retention Rate (%)",
         y = "Ti/Tv Ratio")

# Plot 2: Ti/Tv by quality threshold
p2 <- ggplot(titv_df, aes(x = quality_threshold, y = ti_tv_ratio)) +
    geom_line(color = "#1976D2", size = 1.5) +
    geom_point(size = 3, color = "#1976D2") +
    geom_hline(yintercept = 2.0, linetype = "dashed", color = "red", alpha = 0.5) +
    labs(title = "Ti/Tv Ratio by Quality Threshold",
         subtitle = "Optimal threshold maximizes Ti/Tv",
         x = "Minimum Quality Score",
         y = "Ti/Tv Ratio")

# Plot 3: Variant counts by filter
validation_long <- validation_df %>%
    select(filter_name, snps, indels) %>%
    pivot_longer(cols = c(snps, indels), names_to = "type", values_to = "count")

p3 <- ggplot(validation_long, aes(x = filter_name, y = count, fill = type)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("snps" = "#4CAF50", "indels" = "#FF9800")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Variant Counts by Filter Strategy",
         x = "Filter Strategy",
         y = "Number of Variants",
         fill = "Variant Type")

# Plot 4: MAF spectrum comparison
maf_long <- maf_df %>%
    pivot_longer(cols = c(unfiltered_count, filtered_count), 
                 names_to = "filter_status", 
                 values_to = "count")

p4 <- ggplot(maf_long, aes(x = maf_bin, y = count, fill = filter_status)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("unfiltered_count" = "#9E9E9E", 
                                "filtered_count" = "#2196F3"),
                      labels = c("Filtered", "Unfiltered")) +
    scale_y_log10(labels = comma_format()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "MAF Spectrum: Before and After Filtering",
         x = "Minor Allele Frequency",
         y = "Count (log10)",
         fill = "Status")

# Combine plots
combined_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

# Save plots
ggsave("filter_validation_plots.pdf", combined_plot, width = 14, height = 12)

cat("Filter validation plots generated successfully\\n")
R_SCRIPT
    
    Rscript plot_validation.R
    
    # Clean up temp files
    rm -f temp_*.vcf unfiltered_maf.txt filtered_maf.txt
    
    # Print summary
    echo ""
    echo "Filter Validation Summary:"
    echo "========================="
    
    echo ""
    echo "Filter Performance:"
    cat filter_validation_report.tsv | column -t -s \$'\\t'
    
    echo ""
    echo "Key Findings:"
    best_titv=\$(awk 'NR>1 && \$5!="NA" {print \$5}' filter_validation_report.tsv | sort -rn | head -1)
    echo "  - Best Ti/Tv ratio achieved: \$best_titv"
    echo "  - Original Ti/Tv ratio: \$titv_orig"
    
    if [ "\$(echo "\$best_titv > \$titv_orig + 0.1" | bc)" -eq 1 ]; then
        echo "  ✅ Filtering improves variant quality"
    else
        echo "  ⚠️ Minimal quality improvement from filtering"
    fi
    
    echo ""
    echo "Filter validation complete"
    """
}