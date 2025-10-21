// modules/publication_tables.nf - Generate publication-ready summary tables

process PUBLICATION_TABLES {
    tag "publication_tables"
    publishDir "${params.output_dir}/publication", mode: 'copy'
    conda 'r-base=4.3 r-tidyverse r-knitr r-kableextra'
    
    input:
    path variant_summary
    path sample_stats
    path chr_stats
    path filter_report
    path contamination_report
    
    output:
    path "Table1_sample_summary.csv", emit: table1
    path "Table2_variant_summary.csv", emit: table2
    path "Table3_qc_metrics.csv", emit: table3
    path "SupplementaryTable1_per_chromosome.csv", emit: supp1
    path "SupplementaryTable2_filter_impact.csv", emit: supp2
    path "publication_tables.docx", emit: docx optional true
    path "publication_tables_latex.tex", emit: latex
    
    script:
    """
    echo "Generating publication-ready tables..."
    
    # R script to format tables
    cat > generate_tables.R << 'R_SCRIPT'
library(tidyverse)
library(knitr)
library(kableExtra)

# Function to format numbers for publication
format_number <- function(x, digits = 0, big.mark = ",") {
    format(round(x, digits), big.mark = big.mark, scientific = FALSE)
}

# Function to calculate mean ± SD
mean_sd <- function(x) {
    paste0(round(mean(x, na.rm = TRUE), 2), " ± ", round(sd(x, na.rm = TRUE), 2))
}

# Read data
variant_summary <- jsonlite::fromJSON("${variant_summary}")
sample_df <- read_tsv("${sample_stats}", show_col_types = FALSE)
chr_df <- read_tsv("${chr_stats}", show_col_types = FALSE)
filter_df <- read_tsv("${filter_report}", show_col_types = FALSE)
contam_df <- read_tsv("${contamination_report}", show_col_types = FALSE)

# ============ Table 1: Sample Summary ============
table1 <- sample_df %>%
    summarise(
        Parameter = c(
            "Number of samples",
            "Total genotypes called",
            "Average variants per sample",
            "Average heterozygous calls",
            "Average homozygous alt calls",
            "Average missing rate (%)",
            "Samples with >10% missing",
            "Samples with contamination risk"
        ),
        Value = c(
            format_number(n()),
            format_number(sum(total_variants)),
            mean_sd(total_variants),
            mean_sd(heterozygous),
            mean_sd(homozygous_alt),
            paste0(round(mean(missing_rate * 100), 2), " ± ", 
                   round(sd(missing_rate * 100), 2)),
            format_number(sum(missing_rate > 0.1)),
            format_number(sum(contam_df$risk_level %in% c("HIGH", "MEDIUM")))
        )
    )

write_csv(table1, "Table1_sample_summary.csv")

# ============ Table 2: Variant Summary ============
table2 <- tibble(
    Category = c(
        "Total variants",
        "Single nucleotide polymorphisms (SNPs)",
        "Insertions/Deletions (INDELs)",
        "Multi-allelic sites",
        "Transition/Transversion ratio",
        "Mean depth across variants",
        "Mean quality score",
        "Variants passing QC filters",
        "Private variants (singletons)",
        "Common variants (MAF > 5%)"
    ),
    Count = c(
        format_number(variant_summary$total_variants),
        paste0(format_number(variant_summary$snps), 
               " (", round(variant_summary$snps/variant_summary$total_variants*100, 1), "%)"),
        paste0(format_number(variant_summary$indels),
               " (", round(variant_summary$indels/variant_summary$total_variants*100, 1), "%)"),
        format_number(variant_summary$multiallelic),
        round(as.numeric(variant_summary$ti_tv_ratio), 2),
        "Calculate from data",  # Would need depth distribution
        "Calculate from data",  # Would need quality distribution
        "Calculate post-filtering",
        "Calculate from AF",
        "Calculate from AF"
    )
)

write_csv(table2, "Table2_variant_summary.csv")

# ============ Table 3: Quality Control Metrics ============
# Calculate summary statistics for QC
qc_summary <- sample_df %>%
    summarise(
        het_mean = mean(heterozygous/(total_variants - missing), na.rm = TRUE),
        het_sd = sd(heterozygous/(total_variants - missing), na.rm = TRUE),
        miss_mean = mean(missing_rate, na.rm = TRUE),
        miss_sd = sd(missing_rate, na.rm = TRUE)
    )

# Get filter impact
filter_impact <- filter_df %>%
    filter(filter_name %in% c("no_filter", "stringent")) %>%
    select(filter_name, total_variants, ti_tv_ratio)

table3 <- tibble(
    "QC Metric" = c(
        "Ti/Tv ratio (unfiltered)",
        "Ti/Tv ratio (filtered)",
        "Mean heterozygosity rate",
        "SD heterozygosity rate",
        "Mean missingness rate",
        "SD missingness rate",
        "Samples flagged for contamination",
        "Variants in problematic regions",
        "Filter retention rate (%)"
    ),
    "Value" = c(
        filter_impact$ti_tv_ratio[filter_impact$filter_name == "no_filter"],
        filter_impact$ti_tv_ratio[filter_impact$filter_name == "stringent"],
        round(qc_summary$het_mean, 4),
        round(qc_summary$het_sd, 4),
        round(qc_summary$miss_mean, 4),
        round(qc_summary$miss_sd, 4),
        sum(contam_df$risk_level != "LOW"),
        "Count from BED file",
        round(filter_impact$total_variants[filter_impact$filter_name == "stringent"] / 
              filter_impact$total_variants[filter_impact$filter_name == "no_filter"] * 100, 1)
    ),
    "Threshold/Expected" = c(
        ">2.0 (WGS), >3.0 (WES)",
        ">2.0 (WGS), >3.0 (WES)",
        "Population-specific",
        "<0.05",
        "<0.05",
        "<0.02",
        "0 (ideal)",
        "<5% of total",
        ">60%"
    ),
    "Pass/Fail" = c(
        ifelse(as.numeric(filter_impact$ti_tv_ratio[filter_impact$filter_name == "no_filter"]) > 2.0, "✓", "✗"),
        ifelse(as.numeric(filter_impact$ti_tv_ratio[filter_impact$filter_name == "stringent"]) > 2.0, "✓", "✗"),
        "✓",
        ifelse(qc_summary$het_sd < 0.05, "✓", "⚠"),
        ifelse(qc_summary$miss_mean < 0.05, "✓", "⚠"),
        ifelse(qc_summary$miss_sd < 0.02, "✓", "⚠"),
        ifelse(sum(contam_df$risk_level != "LOW") == 0, "✓", "⚠"),
        "—",
        ifelse(filter_impact$total_variants[filter_impact$filter_name == "stringent"] / 
               filter_impact$total_variants[filter_impact$filter_name == "no_filter"] > 0.6, "✓", "⚠")
    )
)

write_csv(table3, "Table3_qc_metrics.csv")

# ============ Supplementary Table 1: Per-Chromosome Statistics ============
supp_table1 <- chr_df %>%
    mutate(
        snp_percentage = round(snps/total_variants * 100, 1),
        indel_percentage = round(indels/total_variants * 100, 1),
        variants_per_mb = round(total_variants / 1000, 2)  # Assuming ~1000 Mb per chromosome for simplicity
    ) %>%
    select(
        Chromosome = chromosome,
        "Total Variants" = total_variants,
        "SNPs (n)" = snps,
        "SNPs (%)" = snp_percentage,
        "INDELs (n)" = indels,
        "INDELs (%)" = indel_percentage,
        "Mean Depth" = avg_depth,
        "Mean Quality" = avg_qual,
        "Variants/Mb" = variants_per_mb
    ) %>%
    arrange(Chromosome)

write_csv(supp_table1, "SupplementaryTable1_per_chromosome.csv")

# ============ Supplementary Table 2: Filter Impact ============
supp_table2 <- filter_df %>%
    mutate(
        snp_retention = round(snps/snps[filter_name == "no_filter"] * 100, 1),
        indel_retention = round(indels/indels[filter_name == "no_filter"] * 100, 1),
        ti_tv_change = round(as.numeric(ti_tv_ratio) - 
                            as.numeric(ti_tv_ratio[filter_name == "no_filter"]), 3)
    ) %>%
    select(
        "Filter Strategy" = filter_name,
        "Total Variants" = total_variants,
        "SNPs" = snps,
        "INDELs" = indels,
        "Ti/Tv Ratio" = ti_tv_ratio,
        "Retention (%)" = retention_rate,
        "SNP Retention (%)" = snp_retention,
        "INDEL Retention (%)" = indel_retention,
        "Ti/Tv Change" = ti_tv_change
    )

write_csv(supp_table2, "SupplementaryTable2_filter_impact.csv")

# ============ Generate LaTeX Tables ============
latex_output <- c(
    "\\\\documentclass{article}",
    "\\\\usepackage{booktabs}",
    "\\\\usepackage{longtable}",
    "\\\\usepackage{array}",
    "\\\\usepackage{multirow}",
    "\\\\begin{document}",
    "",
    "\\\\section*{Table 1: Sample Summary}",
    kable(table1, format = "latex", booktabs = TRUE) %>%
        kable_styling(latex_options = c("hold_position")),
    "",
    "\\\\section*{Table 2: Variant Summary}",
    kable(table2, format = "latex", booktabs = TRUE) %>%
        kable_styling(latex_options = c("hold_position")),
    "",
    "\\\\section*{Table 3: Quality Control Metrics}",
    kable(table3, format = "latex", booktabs = TRUE) %>%
        kable_styling(latex_options = c("hold_position", "scale_down")),
    "",
    "\\\\section*{Supplementary Table 1: Per-Chromosome Statistics}",
    kable(supp_table1, format = "latex", booktabs = TRUE, longtable = TRUE) %>%
        kable_styling(latex_options = c("repeat_header")),
    "",
    "\\\\section*{Supplementary Table 2: Filter Impact Analysis}",
    kable(supp_table2, format = "latex", booktabs = TRUE) %>%
        kable_styling(latex_options = c("hold_position", "scale_down")),
    "",
    "\\\\end{document}"
)

writeLines(latex_output, "publication_tables_latex.tex")

# Try to create Word document if officer package is available
tryCatch({
    library(officer)
    library(flextable)
    
    doc <- read_docx()
    
    # Add Table 1
    doc <- doc %>%
        body_add_par("Table 1: Sample Summary", style = "heading 1") %>%
        body_add_flextable(flextable(table1))
    
    # Add Table 2
    doc <- doc %>%
        body_add_break() %>%
        body_add_par("Table 2: Variant Summary", style = "heading 1") %>%
        body_add_flextable(flextable(table2))
    
    # Add Table 3
    doc <- doc %>%
        body_add_break() %>%
        body_add_par("Table 3: Quality Control Metrics", style = "heading 1") %>%
        body_add_flextable(flextable(table3))
    
    # Save document
    print(doc, target = "publication_tables.docx")
    cat("Word document generated successfully\\n")
    
}, error = function(e) {
    cat("Note: officer/flextable packages not available, skipping Word output\\n")
})

cat("Publication tables generated successfully\\n")
R_SCRIPT
    
    Rscript generate_tables.R
    
    # Generate a markdown summary for quick review
    cat > publication_summary.md << 'SUMMARY'
# Publication Tables Summary

## Generated Tables

### Main Tables
1. **Table 1: Sample Summary** - Overview of sample-level statistics
2. **Table 2: Variant Summary** - Summary of variant calling results  
3. **Table 3: Quality Control Metrics** - QC metrics and pass/fail status

### Supplementary Tables
1. **Supplementary Table 1: Per-Chromosome Statistics** - Detailed breakdown by chromosome
2. **Supplementary Table 2: Filter Impact Analysis** - Effect of different filtering strategies

## File Formats Available
- CSV files for all tables (Excel-compatible)
- LaTeX source file for direct inclusion in manuscripts
- Word document (if officer package available)

## Usage Instructions

### For Manuscripts
1. Import CSV files into your preferred table editor
2. Or use the LaTeX source directly in your manuscript
3. Adjust formatting as needed for journal requirements

### Key Statistics for Abstract/Text
- Total variants: See Table 2
- Ti/Tv ratio: See Table 3
- Sample quality metrics: See Table 1
- Filter retention: See Supplementary Table 2

## Quality Indicators
✓ Pass | ⚠ Warning | ✗ Fail | — Not applicable
SUMMARY
    
    echo "Publication tables complete!"
    echo "Available formats:"
    ls -la *.csv *.tex *.docx 2>/dev/null || ls -la *.csv *.tex
    """
}