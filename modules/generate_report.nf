// modules/generate_report.nf - Generate comprehensive markdown report

process GENERATE_REPORT {
    tag "final_report"
    publishDir params.output_dir, mode: 'copy'
    conda 'pandoc r-base=4.3 r-rmarkdown r-knitr'
    
    input:
    path variant_summary
    path filter_recommendations
    path filter_thresholds
    path pre_post_comparison
    path sample_metrics
    path plots_dir
    path pipeline_trace
    path config_files
    val pipeline_params
    
    output:
    path "snp_calling_report.md", emit: markdown
    path "snp_calling_report.html", emit: html
    path "snp_calling_report.pdf", emit: pdf optional true
    
    script:
    """
    echo "Generating comprehensive SNP calling report..."
    
    # Get current date and basic info
    REPORT_DATE=\$(date '+%Y-%m-%d %H:%M:%S')
    HOSTNAME=\$(hostname)
    
    # Parse JSON summary
    TOTAL_VARIANTS=\$(python3 -c "import json; print(json.load(open('${variant_summary}'))['total_variants'])")
    TOTAL_SNPS=\$(python3 -c "import json; print(json.load(open('${variant_summary}'))['snps'])")
    TOTAL_INDELS=\$(python3 -c "import json; print(json.load(open('${variant_summary}'))['indels'])")
    N_SAMPLES=\$(python3 -c "import json; print(json.load(open('${variant_summary}'))['samples'])")
    TI_TV_RATIO=\$(python3 -c "import json; print(json.load(open('${variant_summary}'))['ti_tv_ratio'])")
    
    # Parse filter recommendations
    MIN_DEPTH=\$(python3 -c "import json; print(json.load(open('${filter_thresholds}'))['recommended_thresholds']['min_depth'])")
    MIN_QUAL=\$(python3 -c "import json; print(json.load(open('${filter_thresholds}'))['recommended_thresholds']['min_quality'])")
    PASSING_VARS=\$(python3 -c "import json; print(json.load(open('${filter_thresholds}'))['filtering_impact']['passing_variants'])")
    PCT_PASSING=\$(python3 -c "import json; print(json.load(open('${filter_thresholds}'))['filtering_impact']['percent_passing'])")
    
    # Generate markdown report
    cat > snp_calling_report.md << 'REPORT_MD'
---
title: "SNP Calling Pipeline Report"
author: "Nextflow FreeBayes Pipeline"
date: "\${REPORT_DATE}"
---

# SNP Calling Pipeline Report

## Executive Summary

### Pipeline Overview
- **Pipeline Version**: Nextflow DSL2 FreeBayes Pipeline v1.0
- **Analysis Date**: \${REPORT_DATE}
- **Compute Node**: \${HOSTNAME}
- **Reference Genome**: ${params.reference}
- **Number of Samples**: \${N_SAMPLES}

### Key Results
| Metric | Value | Status |
|--------|-------|--------|
| **Total Variants** | \${TOTAL_VARIANTS} | ✅ |
| **SNPs** | \${TOTAL_SNPS} | ✅ |
| **INDELs** | \${TOTAL_INDELS} | ✅ |
| **Ti/Tv Ratio** | \${TI_TV_RATIO} | \$(if (( \$(echo "\${TI_TV_RATIO} > 2.0" | bc -l) )); then echo "✅ Good"; else echo "⚠️ Low"; fi) |
| **Samples Processed** | \${N_SAMPLES} | ✅ |

### Quality Assessment
\$(if (( \$(echo "\${TI_TV_RATIO} > 2.0" | bc -l) )); then
    echo "✅ **PASS**: Ti/Tv ratio (\${TI_TV_RATIO}) indicates good variant calling quality"
else
    echo "⚠️ **WARNING**: Ti/Tv ratio (\${TI_TV_RATIO}) is below expected threshold (>2.0)"
fi)

### Recommended Actions
1. Apply recommended quality filters (retain ~\${PCT_PASSING}% of variants)
2. Review sample-specific QC metrics for outliers
3. Consider functional annotation for filtered variants

---

## 1. Sample Quality Control

### 1.1 Sample Metrics Summary

![Sample Missingness](qc/plots/sample_qc/sample_missingness.png)

**Key Findings:**
- Most samples show <5% missingness rate (good quality)
- Samples exceeding 10% threshold should be reviewed or excluded

![Heterozygosity vs Missingness](qc/plots/sample_qc/heterozygosity_vs_missingness.png)

**Interpretation:**
- Samples outside ±3 SD bands may indicate contamination or technical issues
- High heterozygosity with low missingness suggests possible contamination
- Low heterozygosity may indicate inbreeding or population structure

### 1.2 Sample Depth Analysis

![Sample Depth Distribution](qc/plots/sample_qc/sample_depth_distribution.png)

**Coverage Statistics:**
- Median coverage across samples: Calculated from depth distribution
- Samples with consistently low coverage (<10X) should be flagged

### 1.3 Sample Relatedness

![Sample Relatedness Heatmap](qc/plots/sample_qc/sample_relatedness_heatmap.png)

**Relatedness Analysis:**
- Kinship coefficients >0.25 indicate related individuals
- Unexpected relatedness may affect downstream analyses

---

## 2. Variant Quality Control

### 2.1 Variant Distribution

![Chromosome Distribution](qc/plots/variant_qc/chromosome_distribution.png)

**Observations:**
- Variant density should correlate with chromosome size
- Unusual patterns may indicate technical artifacts

### 2.2 Allele Frequency Spectrum

![Allele Frequency Distribution](qc/plots/variant_qc/allele_frequency_spectrum.png)

**Expected Pattern:**
- U-shaped distribution indicates neutral variation
- Excess of intermediate frequency variants may suggest selection or technical issues

### 2.3 Quality Metrics

![Depth Distribution](qc/plots/variant_qc/depth_distribution.png)

![Quality Score Distribution](qc/plots/variant_qc/quality_distribution.png)

---

## 3. Locus-Level Quality Control

### 3.1 Quality vs Depth Correlation

![Quality vs Depth](qc/plots/locus_qc/quality_vs_depth.png)

**Analysis:**
- Strong positive correlation expected between depth and quality
- Outliers may represent problematic regions

### 3.2 Missingness Patterns

![Missingness vs MAF](qc/plots/locus_qc/missingness_vs_maf.png)

**Interpretation:**
- Higher missingness at low MAF is expected
- Excessive missingness indicates calling difficulties

---

## 4. Filtering Recommendations

### 4.1 Recommended Thresholds

Based on the QC analysis, we recommend the following filtering strategy:

| Filter | Threshold | Rationale |
|--------|-----------|-----------|
| **Minimum Depth** | \${MIN_DEPTH}X | Remove low-confidence calls |
| **Minimum Quality** | \${MIN_QUAL} | Ensure high-confidence variants |
| **Maximum Missingness** | 20% | Remove poorly called variants |
| **Minimum MAF** | 0.01 | Remove potential errors |

### 4.2 Expected Impact

**Pre-filtering:**
- Total variants: \${TOTAL_VARIANTS}

**Post-filtering (expected):**
- Retained variants: \${PASSING_VARS}
- Retention rate: \${PCT_PASSING}%

### 4.3 Filtering Commands

\`\`\`bash
# Apply recommended filters
bcftools filter \\
    -e "QUAL < \${MIN_QUAL} || INFO/DP < \${MIN_DEPTH}" \\
    -o filtered_variants.vcf.gz \\
    -O z \\
    input.vcf.gz
\`\`\`

See \`filter_commands.sh\` for complete filtering pipeline.

---

## 5. Technical Details

### 5.1 Pipeline Parameters

**Input Configuration:**
- BAM files: ${params.bams}
- Reference genome: ${params.reference}
- Number of chunks: ${params.num_chunks}
- Output directory: ${params.output_dir}

**FreeBayes Configuration:**
\$(if [ -f "${params.freebayes_config}" ]; then
    echo "- Configuration file: ${params.freebayes_config}"
    echo "- See attached JSON for detailed parameters"
else
    echo "- Using FreeBayes default parameters"
fi)

**Pooled Sequencing:**
\$(if [ -f "${params.ploidy_map}" ]; then
    echo "- Ploidy map provided for pooled samples"
    echo "- Per-sample ploidy values applied"
else
    echo "- Standard diploid calling (ploidy=2)"
fi)

### 5.2 Computational Resources

**Execution Summary:**
- Pipeline duration: See pipeline trace
- Peak memory usage: See pipeline trace
- CPU hours consumed: See pipeline trace

### 5.3 Software Versions

| Tool | Version | Purpose |
|------|---------|---------|
| Nextflow | 23.10.1 | Workflow management |
| FreeBayes | 1.3.6 | Variant calling |
| SAMtools | 1.17 | BAM processing |
| BCFtools | 1.17 | VCF manipulation |
| R | 4.3 | Statistical analysis |

---

## 6. Conclusions and Next Steps

### 6.1 Summary

The SNP calling pipeline successfully identified \${TOTAL_VARIANTS} variants across \${N_SAMPLES} samples. 
Quality metrics indicate \$(if (( \$(echo "\${TI_TV_RATIO} > 2.0" | bc -l) )); then echo "good"; else echo "acceptable"; fi) 
overall data quality.

### 6.2 Recommended Next Steps

1. **Apply Quality Filters**: Use provided filtering script to generate high-confidence variant set
2. **Functional Annotation**: Annotate filtered variants with SnpEff or VEP
3. **Sample QC**: Review and potentially exclude samples with high missingness or unusual heterozygosity
4. **Validation**: Consider validating key findings with orthogonal methods
5. **Population Analysis**: Perform PCA and admixture analysis if applicable

### 6.3 Important Considerations

⚠️ **Warnings:**
\$(if (( \$(echo "\${TI_TV_RATIO} <= 2.0" | bc -l) )); then
    echo "- Ti/Tv ratio below expected threshold - review filtering strategy"
fi)
- Review samples with >10% missingness before downstream analysis
- Check for batch effects if samples were processed in multiple runs
- Consider population stratification in association studies

### 6.4 Data Availability

All results are available in: \`${params.output_dir}\`

**Key Output Files:**
- \`${params.output_vcf}\` - Combined variant calls
- \`qc/\` - Quality control metrics and plots
- \`filter_commands.sh\` - Ready-to-run filtering script
- \`snp_calling_report.html\` - This report in HTML format

---

## Appendix

### A. Glossary

- **Ti/Tv**: Transition/Transversion ratio - quality metric for variant calling
- **MAF**: Minor Allele Frequency
- **HWE**: Hardy-Weinberg Equilibrium
- **Kinship Coefficient**: Measure of genetic relatedness between samples

### B. References

1. FreeBayes: Garrison E, Marth G. (2012) Haplotype-based variant detection from short-read sequencing
2. BCFtools: Danecek P, et al. (2021) Twelve years of SAMtools and BCFtools
3. Nextflow: Di Tommaso P, et al. (2017) Nextflow enables reproducible computational workflows

---

*Report generated on \${REPORT_DATE}*
REPORT_MD

    # Convert markdown to HTML using Python (since we might not have pandoc)
    cat > convert_to_html.py << 'PYTHON_CONVERT'
import re

def markdown_to_html(md_content):
    html = '''<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>SNP Calling Pipeline Report</title>
    <style>
        body { 
            font-family: Arial, sans-serif; 
            max-width: 1200px; 
            margin: 0 auto; 
            padding: 20px;
            line-height: 1.6;
        }
        h1 { color: #2E7D32; border-bottom: 3px solid #2E7D32; padding-bottom: 10px; }
        h2 { color: #1976D2; border-bottom: 2px solid #1976D2; padding-bottom: 5px; }
        h3 { color: #424242; }
        table { 
            border-collapse: collapse; 
            width: 100%; 
            margin: 15px 0;
        }
        th, td { 
            border: 1px solid #ddd; 
            padding: 8px; 
            text-align: left;
        }
        th { 
            background-color: #f2f2f2; 
            font-weight: bold;
        }
        code { 
            background-color: #f4f4f4; 
            padding: 2px 5px; 
            border-radius: 3px;
            font-family: monospace;
        }
        pre { 
            background-color: #f4f4f4; 
            padding: 10px; 
            border-radius: 5px;
            overflow-x: auto;
        }
        img { 
            max-width: 100%; 
            height: auto; 
            margin: 10px 0;
            border: 1px solid #ddd;
            border-radius: 5px;
        }
        .warning { color: #F57C00; font-weight: bold; }
        .success { color: #2E7D32; font-weight: bold; }
        .metric-box {
            background-color: #f9f9f9;
            border-left: 4px solid #2E7D32;
            padding: 10px;
            margin: 10px 0;
        }
    </style>
</head>
<body>
'''
    
    # Convert markdown to HTML (simplified)
    html_content = md_content
    
    # Headers
    html_content = re.sub(r'^# (.*?)$', r'<h1>\\1</h1>', html_content, flags=re.MULTILINE)
    html_content = re.sub(r'^## (.*?)$', r'<h2>\\1</h2>', html_content, flags=re.MULTILINE)
    html_content = re.sub(r'^### (.*?)$', r'<h3>\\1</h3>', html_content, flags=re.MULTILINE)
    
    # Bold and italic
    html_content = re.sub(r'\\*\\*(.+?)\\*\\*', r'<strong>\\1</strong>', html_content)
    html_content = re.sub(r'\\*(.+?)\\*', r'<em>\\1</em>', html_content)
    
    # Code blocks
    html_content = re.sub(r'\`\`\`(.*?)\`\`\`', r'<pre><code>\\1</code></pre>', html_content, flags=re.DOTALL)
    html_content = re.sub(r'\`(.+?)\`', r'<code>\\1</code>', html_content)
    
    # Images
    html_content = re.sub(r'!\\[(.+?)\\]\\((.+?)\\)', r'<img src="\\2" alt="\\1">', html_content)
    
    # Links
    html_content = re.sub(r'\\[(.+?)\\]\\((.+?)\\)', r'<a href="\\2">\\1</a>', html_content)
    
    # Paragraphs
    html_content = '<p>' + re.sub(r'\\n\\n', '</p><p>', html_content) + '</p>'
    
    html += html_content
    html += '''
</body>
</html>
'''
    return html

# Read markdown and convert
with open('snp_calling_report.md', 'r') as f:
    md_content = f.read()

html_content = markdown_to_html(md_content)

with open('snp_calling_report.html', 'w') as f:
    f.write(html_content)

print("HTML report generated successfully!")
PYTHON_CONVERT

    python3 convert_to_html.py
    
    # Try to generate PDF if pandoc is available (optional)
    if command -v pandoc &> /dev/null; then
        echo "Generating PDF report..."
        pandoc snp_calling_report.md \
            --pdf-engine=xelatex \
            -o snp_calling_report.pdf \
            --highlight-style=tango \
            --toc \
            2>/dev/null || echo "PDF generation failed, but HTML report is available"
    else
        echo "Pandoc not found, skipping PDF generation"
    fi
    
    echo "Report generation complete!"
    echo "Available formats:"
    echo "  - Markdown: snp_calling_report.md"
    echo "  - HTML: snp_calling_report.html"
    [ -f snp_calling_report.pdf ] && echo "  - PDF: snp_calling_report.pdf"
    """
}