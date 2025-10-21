// modules/locus_qc.nf - Calculate locus/variant-level QC metrics

process LOCUS_QC {
    tag "locus_qc"
    publishDir "${params.output_dir}/qc", mode: 'copy'
    
    input:
    path vcf
    path vcf_index
    
    output:
    path "locus_qc_metrics.tsv", emit: metrics
    path "problematic_regions.bed", emit: problem_regions
    path "high_missingness_loci.tsv", emit: high_miss
    path "depth_outlier_loci.tsv", emit: depth_outliers
    path "quality_outlier_loci.tsv", emit: qual_outliers
    path "allele_balance_issues.tsv", emit: ab_issues
    path "hardy_weinberg_failures.tsv", emit: hwe_failures
    
    script:
    """
    echo "Calculating locus-level QC metrics"
    
    # Calculate comprehensive locus metrics
    echo "chrom\tpos\tref\talt\ttype\tdepth\tquality\tmissingness_rate\tallele_freq\theterozygosity\thomozygosity\tallele_balance_mean\thwe_pvalue" > locus_qc_metrics.tsv
    
    # Extract variant information with calculated metrics
    bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%TYPE\\t%INFO/DP\\t%QUAL\\n' ${vcf} > variant_basic_info.txt
    
    # Get sample count for calculations
    n_samples=\$(bcftools query -l ${vcf} | wc -l)
    
    # Process each variant
    line_num=0
    bcftools view -H ${vcf} | while read -r line; do
        line_num=\$((line_num + 1))
        
        # Get basic info
        chrom=\$(echo "\$line" | cut -f1)
        pos=\$(echo "\$line" | cut -f2)
        ref=\$(echo "\$line" | cut -f4)
        alt=\$(echo "\$line" | cut -f5)
        qual=\$(echo "\$line" | cut -f6)
        info=\$(echo "\$line" | cut -f8)
        
        # Extract depth from INFO
        depth=\$(echo "\$info" | grep -oP 'DP=\\K[0-9]+' || echo "0")
        
        # Determine variant type
        if [ \${#ref} -eq 1 ] && [ \${#alt} -eq 1 ]; then
            vtype="SNP"
        elif [ \${#ref} -ne \${#alt} ]; then
            vtype="INDEL"
        else
            vtype="OTHER"
        fi
        
        # Get genotype data for this position
        bcftools view -H -r \${chrom}:\${pos}-\${pos} ${vcf} | cut -f10- | tr '\\t' '\\n' > temp_genotypes.txt
        
        # Calculate missingness
        missing=\$(grep -c "^\\./\\.\\|^\\.|\\." temp_genotypes.txt || true)
        miss_rate=\$(echo "scale=4; \$missing / \$n_samples" | bc)
        
        # Calculate allele frequency and genotype counts
        n_hom_ref=\$(grep -c "^0/0\\|^0|0" temp_genotypes.txt || true)
        n_het=\$(grep -c "^0/1\\|^1/0\\|^0|1\\|^1|0" temp_genotypes.txt || true)
        n_hom_alt=\$(grep -c "^1/1\\|^1|1" temp_genotypes.txt || true)
        
        n_called=\$((n_hom_ref + n_het + n_hom_alt))
        
        if [ \$n_called -gt 0 ]; then
            # Calculate allele frequency
            n_ref_alleles=\$((2 * n_hom_ref + n_het))
            n_alt_alleles=\$((2 * n_hom_alt + n_het))
            total_alleles=\$((n_ref_alleles + n_alt_alleles))
            
            if [ \$total_alleles -gt 0 ]; then
                af=\$(echo "scale=4; \$n_alt_alleles / \$total_alleles" | bc)
            else
                af=0
            fi
            
            # Calculate heterozygosity and homozygosity rates
            het_rate=\$(echo "scale=4; \$n_het / \$n_called" | bc)
            hom_rate=\$(echo "scale=4; (\$n_hom_ref + \$n_hom_alt) / \$n_called" | bc)
            
            # Hardy-Weinberg test (chi-square approximation)
            if [ \$n_called -ge 20 ]; then  # Only test if sufficient samples
                p=\$(echo "scale=6; \$n_ref_alleles / \$total_alleles" | bc)
                q=\$(echo "scale=6; 1 - \$p" | bc)
                
                # Expected counts under HWE
                exp_hom_ref=\$(echo "scale=2; \$n_called * \$p * \$p" | bc)
                exp_het=\$(echo "scale=2; \$n_called * 2 * \$p * \$q" | bc)
                exp_hom_alt=\$(echo "scale=2; \$n_called * \$q * \$q" | bc)
                
                # Chi-square statistic (simplified p-value calculation)
                if (( \$(echo "\$exp_hom_ref > 0 && \$exp_het > 0 && \$exp_hom_alt > 0" | bc -l) )); then
                    chi_sq=\$(echo "scale=4; \\
                        ((\$n_hom_ref - \$exp_hom_ref)^2 / \$exp_hom_ref) + \\
                        ((\$n_het - \$exp_het)^2 / \$exp_het) + \\
                        ((\$n_hom_alt - \$exp_hom_alt)^2 / \$exp_hom_alt)" | bc)
                    
                    # Very rough p-value approximation (chi-square with 1 df)
                    if (( \$(echo "\$chi_sq > 10.83" | bc -l) )); then
                        hwe_pval="<0.001"
                    elif (( \$(echo "\$chi_sq > 6.64" | bc -l) )); then
                        hwe_pval="<0.01"
                    elif (( \$(echo "\$chi_sq > 3.84" | bc -l) )); then
                        hwe_pval="<0.05"
                    else
                        hwe_pval=">0.05"
                    fi
                else
                    hwe_pval="NA"
                fi
            else
                hwe_pval="NA"
            fi
        else
            af=0
            het_rate=0
            hom_rate=0
            hwe_pval="NA"
        fi
        
        # Calculate allele balance for heterozygotes (simplified - would need AD field)
        ab_mean="NA"  # Would require allele depth information
        
        echo -e "\$chrom\\t\$pos\\t\$ref\\t\$alt\\t\$vtype\\t\$depth\\t\$qual\\t\$miss_rate\\t\$af\\t\$het_rate\\t\$hom_rate\\t\$ab_mean\\t\$hwe_pval" >> locus_qc_metrics.tsv
        
        # Process only first 10000 variants for speed (remove this limit for production)
        if [ \$line_num -ge 10000 ]; then
            echo "Processed first 10000 variants for QC metrics..."
            break
        fi
    done
    
    rm -f temp_genotypes.txt
    
    # Identify problematic loci based on various criteria
    
    # High missingness loci (>20% missing)
    echo "chrom\tpos\tref\talt\tmissingness_rate" > high_missingness_loci.tsv
    awk '\$8 > 0.2 && NR > 1 {print \$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$8}' locus_qc_metrics.tsv >> high_missingness_loci.tsv
    
    # Depth outliers (too low <5 or too high >99th percentile)
    echo "chrom\tpos\tref\talt\tdepth\tissue" > depth_outlier_loci.tsv
    
    # Calculate depth percentiles
    depth_p99=\$(awk 'NR > 1 && \$6 != "." {print \$6}' locus_qc_metrics.tsv | \
                 sort -n | awk '{all[NR] = \$0} END {print all[int(NR*0.99)]}')
    
    if [ -z "\$depth_p99" ]; then depth_p99=1000; fi
    
    awk -v p99="\$depth_p99" 'NR > 1 {
        if(\$6 < 5) print \$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$6"\\tlow_depth"
        else if(\$6 > p99) print \$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$6"\\thigh_depth"
    }' locus_qc_metrics.tsv >> depth_outlier_loci.tsv
    
    # Quality outliers (QUAL < 20)
    echo "chrom\tpos\tref\talt\tquality" > quality_outlier_loci.tsv
    awk '\$7 < 20 && NR > 1 && \$7 != "." {print \$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$7}' locus_qc_metrics.tsv >> quality_outlier_loci.tsv
    
    # Allele balance issues (for future implementation with AD field)
    echo "chrom\tpos\tref\talt\tallele_balance_mean\tissue" > allele_balance_issues.tsv
    # This would require AD (allele depth) field processing
    
    # Hardy-Weinberg failures (p < 0.001)
    echo "chrom\tpos\tref\talt\thwe_pvalue\thet_obs\thet_exp" > hardy_weinberg_failures.tsv
    awk '\$13 == "<0.001" && NR > 1 {print \$1"\\t"\$2"\\t"\$3"\\t"\$4"\\t"\$13"\\t"\$10"\\tNA"}' locus_qc_metrics.tsv >> hardy_weinberg_failures.tsv
    
    # Create BED file of problematic regions (combine all issues)
    echo "Creating BED file of problematic regions..."
    
    # Start with empty BED file
    > problematic_regions.bed
    
    # Add high missingness regions
    awk 'NR > 1 {print \$1"\\t"(\$2-1)"\\t"\$2"\\thigh_missingness"}' high_missingness_loci.tsv >> problematic_regions.bed
    
    # Add depth outlier regions
    awk 'NR > 1 {print \$1"\\t"(\$2-1)"\\t"\$2"\\t"\$6}' depth_outlier_loci.tsv >> problematic_regions.bed
    
    # Add quality outlier regions
    awk 'NR > 1 {print \$1"\\t"(\$2-1)"\\t"\$2"\\tlow_quality"}' quality_outlier_loci.tsv >> problematic_regions.bed
    
    # Add HWE failure regions
    awk 'NR > 1 {print \$1"\\t"(\$2-1)"\\t"\$2"\\thwe_failure"}' hardy_weinberg_failures.tsv >> problematic_regions.bed
    
    # Sort and merge overlapping regions
    if [ -s problematic_regions.bed ]; then
        sort -k1,1 -k2,2n problematic_regions.bed | \
        bedtools merge -i stdin -c 4 -o distinct > problematic_regions_merged.bed
        mv problematic_regions_merged.bed problematic_regions.bed
    fi
    
    # Summary statistics
    echo ""
    echo "Locus QC Summary:"
    echo "=================="
    total_loci=\$(tail -n +2 locus_qc_metrics.tsv | wc -l)
    echo "Total loci analyzed: \$total_loci"
    
    high_miss=\$(tail -n +2 high_missingness_loci.tsv | wc -l)
    echo "High missingness loci (>20%): \$high_miss"
    
    depth_issues=\$(tail -n +2 depth_outlier_loci.tsv | wc -l)
    echo "Depth outlier loci: \$depth_issues"
    
    qual_issues=\$(tail -n +2 quality_outlier_loci.tsv | wc -l)
    echo "Low quality loci (QUAL<20): \$qual_issues"
    
    hwe_issues=\$(tail -n +2 hardy_weinberg_failures.tsv | wc -l)
    echo "Hardy-Weinberg failures (p<0.001): \$hwe_issues"
    
    if [ -s problematic_regions.bed ]; then
        problem_regions=\$(wc -l < problematic_regions.bed)
        echo "Total problematic regions: \$problem_regions"
    fi
    
    echo "=================="
    echo "Locus QC complete"
    """
}