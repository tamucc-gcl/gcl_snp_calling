// modules/contamination_check.nf - Detect potential contamination in samples

process CONTAMINATION_CHECK {
    tag "contamination_check"
    publishDir "${params.output_dir}/qc/contamination", mode: 'copy'
    
    input:
    path vcf
    path vcf_index
    path sample_heterozygosity
    
    output:
    path "contamination_report.tsv", emit: report
    path "contamination_suspects.txt", emit: suspects
    path "allele_balance_per_sample.tsv", emit: allele_balance
    path "contamination_plots.pdf", emit: plots
    
    script:
    """
    echo "Checking for potential contamination..."
    
    # Extract allele balance information for heterozygous sites
    echo "sample\tmean_allele_balance\tsd_allele_balance\tskew\tmedian_ab\toutlier_sites" > allele_balance_per_sample.tsv
    
    # Get sample list
    bcftools query -l ${vcf} > samples.txt
    
    # For each sample, calculate allele balance statistics
    while read sample; do
        echo "Processing sample: \$sample"
        
        # Extract allele depths for heterozygous sites (requires AD field)
        bcftools query -s \$sample -f '[%GT\\t%AD\\n]' ${vcf} | \
        awk -F'\\t' '\$1 ~ /0\\/1|1\\/0/ && \$2 != "." {
            split(\$2, depths, ",")
            if(depths[1] > 0 && depths[2] > 0) {
                total = depths[1] + depths[2]
                ab = depths[2] / total
                print ab
            }
        }' > \${sample}_allele_balance.txt
        
        # Calculate statistics if AD field exists and has data
        if [ -s \${sample}_allele_balance.txt ]; then
            python3 << PYTHON_SCRIPT
import numpy as np
import sys

sample = "\$sample"
ab_values = []

with open(f"{sample}_allele_balance.txt", 'r') as f:
    for line in f:
        try:
            ab_values.append(float(line.strip()))
        except:
            continue

if len(ab_values) > 0:
    mean_ab = np.mean(ab_values)
    sd_ab = np.std(ab_values)
    median_ab = np.median(ab_values)
    
    # Calculate skewness (simplified)
    if sd_ab > 0:
        skew = 3 * (mean_ab - median_ab) / sd_ab
    else:
        skew = 0
    
    # Count outlier sites (AB < 0.2 or AB > 0.8)
    outliers = sum(1 for ab in ab_values if ab < 0.2 or ab > 0.8)
    
    print(f"{sample}\\t{mean_ab:.3f}\\t{sd_ab:.3f}\\t{skew:.3f}\\t{median_ab:.3f}\\t{outliers}")
else:
    print(f"{sample}\\t0\\t0\\t0\\t0\\t0")
PYTHON_SCRIPT
        else
            # No AD field or no heterozygous sites
            echo -e "\$sample\\t0.5\\t0\\t0\\t0.5\\t0"
        fi
        
        rm -f \${sample}_allele_balance.txt
    done < samples.txt >> allele_balance_per_sample.tsv
    
    # Analyze contamination indicators
    python3 << 'CONTAMINATION_ANALYSIS'
import pandas as pd
import numpy as np
import json

# Read heterozygosity data
het_df = pd.read_csv("${sample_heterozygosity}", sep='\\t')

# Read allele balance data
ab_df = pd.read_csv("allele_balance_per_sample.tsv", sep='\\t')

# Merge data
merged_df = het_df.merge(ab_df, on='sample', how='left')

# Calculate contamination indicators
contamination_scores = []

for _, row in merged_df.iterrows():
    score = 0
    flags = []
    
    # 1. Excessive heterozygosity (>3 SD from mean)
    het_mean = het_df['heterozygosity_rate'].mean()
    het_sd = het_df['heterozygosity_rate'].std()
    if row['heterozygosity_rate'] > het_mean + 3 * het_sd:
        score += 3
        flags.append("high_heterozygosity")
    elif row['heterozygosity_rate'] > het_mean + 2 * het_sd:
        score += 1
        flags.append("elevated_heterozygosity")
    
    # 2. Skewed allele balance (deviation from 0.5)
    if abs(row['mean_allele_balance'] - 0.5) > 0.1:
        score += 2
        flags.append("skewed_allele_balance")
    
    # 3. High variation in allele balance
    if row['sd_allele_balance'] > 0.2:
        score += 2
        flags.append("variable_allele_balance")
    
    # 4. Many outlier allele balance sites
    if row['outlier_sites'] > 100:
        score += 2
        flags.append("many_ab_outliers")
    
    # 5. Negative inbreeding coefficient (excess heterozygosity)
    if row['inbreeding_coefficient'] < -0.1:
        score += 2
        flags.append("negative_inbreeding")
    
    contamination_scores.append({
        'sample': row['sample'],
        'contamination_score': score,
        'risk_level': 'HIGH' if score >= 5 else 'MEDIUM' if score >= 3 else 'LOW',
        'heterozygosity_rate': row['heterozygosity_rate'],
        'mean_allele_balance': row['mean_allele_balance'],
        'sd_allele_balance': row['sd_allele_balance'],
        'inbreeding_coefficient': row['inbreeding_coefficient'],
        'flags': '|'.join(flags) if flags else 'none'
    })

# Create contamination report
report_df = pd.DataFrame(contamination_scores)
report_df = report_df.sort_values('contamination_score', ascending=False)
report_df.to_csv('contamination_report.tsv', sep='\\t', index=False)

# Identify high-risk samples
suspects = report_df[report_df['risk_level'].isin(['HIGH', 'MEDIUM'])]['sample'].tolist()
with open('contamination_suspects.txt', 'w') as f:
    if suspects:
        f.write('\\n'.join(suspects))
        print(f"Identified {len(suspects)} samples with potential contamination")
    else:
        f.write('No samples with high contamination risk detected')
        print("No samples with high contamination risk detected")

# Print summary
print(f"\\nContamination Risk Summary:")
print(f"HIGH risk:   {len(report_df[report_df['risk_level'] == 'HIGH'])} samples")
print(f"MEDIUM risk: {len(report_df[report_df['risk_level'] == 'MEDIUM'])} samples")
print(f"LOW risk:    {len(report_df[report_df['risk_level'] == 'LOW'])} samples")
CONTAMINATION_ANALYSIS
    
    # Generate contamination plots using R
    cat > plot_contamination.R << 'R_SCRIPT'
library(ggplot2)
library(tidyverse)
library(gridExtra)

# Read data
contamination_df <- read_tsv("contamination_report.tsv", show_col_types = FALSE)
ab_df <- read_tsv("allele_balance_per_sample.tsv", show_col_types = FALSE)

# Set theme
theme_set(theme_bw(base_size = 12))

# Color scheme for risk levels
risk_colors <- c("LOW" = "#4CAF50", "MEDIUM" = "#FF9800", "HIGH" = "#F44336")

# Plot 1: Contamination score distribution
p1 <- ggplot(contamination_df, aes(x = reorder(sample, contamination_score), 
                                   y = contamination_score, 
                                   fill = risk_level)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = risk_colors) +
    coord_flip() +
    labs(title = "Contamination Risk Scores by Sample",
         x = "Sample",
         y = "Contamination Score",
         fill = "Risk Level") +
    theme(axis.text.y = element_text(size = 6))

# Plot 2: Heterozygosity vs Allele Balance
p2 <- ggplot(contamination_df, aes(x = heterozygosity_rate, 
                                   y = mean_allele_balance,
                                   color = risk_level)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_hline(yintercept = 0.5, linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = mean(contamination_df$heterozygosity_rate) + 
                           2*sd(contamination_df$heterozygosity_rate), 
               linetype = "dashed", alpha = 0.5) +
    scale_color_manual(values = risk_colors) +
    labs(title = "Contamination Indicators",
         subtitle = "Expected: heterozygosity near mean, allele balance near 0.5",
         x = "Heterozygosity Rate",
         y = "Mean Allele Balance",
         color = "Risk Level")

# Plot 3: Allele balance standard deviation
p3 <- ggplot(contamination_df, aes(x = sd_allele_balance, fill = risk_level)) +
    geom_histogram(bins = 30, alpha = 0.7, position = "stack") +
    scale_fill_manual(values = risk_colors) +
    geom_vline(xintercept = 0.15, linetype = "dashed", color = "red", alpha = 0.5) +
    labs(title = "Allele Balance Variability",
         subtitle = "High SD indicates potential contamination",
         x = "Standard Deviation of Allele Balance",
         y = "Count",
         fill = "Risk Level")

# Plot 4: Inbreeding coefficient vs heterozygosity
p4 <- ggplot(contamination_df, aes(x = inbreeding_coefficient, 
                                   y = heterozygosity_rate,
                                   color = risk_level)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_color_manual(values = risk_colors) +
    labs(title = "Inbreeding vs Heterozygosity",
         subtitle = "Negative F with high heterozygosity suggests contamination",
         x = "Inbreeding Coefficient (F)",
         y = "Heterozygosity Rate",
         color = "Risk Level")

# Combine plots
combined_plot <- grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

# Save plots
ggsave("contamination_plots.pdf", combined_plot, width = 14, height = 12)

cat("Contamination plots generated successfully\\n")
R_SCRIPT
    
    Rscript plot_contamination.R
    
    # Generate summary statistics
    echo ""
    echo "Contamination Check Summary:"
    echo "============================"
    
    n_high=\$(grep -c "HIGH" contamination_report.tsv || echo "0")
    n_medium=\$(grep -c "MEDIUM" contamination_report.tsv || echo "0")
    n_low=\$(grep -c "LOW" contamination_report.tsv || echo "0")
    
    echo "Risk Level Distribution:"
    echo "  HIGH:   \$n_high samples"
    echo "  MEDIUM: \$n_medium samples"
    echo "  LOW:    \$n_low samples"
    
    if [ \$n_high -gt 0 ] || [ \$n_medium -gt 0 ]; then
        echo ""
        echo "⚠️ WARNING: Potential contamination detected!"
        echo "Review contamination_report.tsv for details"
        echo "Consider excluding samples listed in contamination_suspects.txt"
    else
        echo ""
        echo "✅ No significant contamination detected"
    fi
    
    echo ""
    echo "Contamination check complete"
    """
}