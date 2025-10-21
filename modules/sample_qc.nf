// modules/sample_qc.nf - Calculate sample-level QC metrics

process SAMPLE_QC {
    tag "sample_qc"
    publishDir "${params.output_dir}/qc", mode: 'copy'
    
    input:
    path vcf
    path vcf_index
    path ploidy_map  // Optional ploidy map for pooled samples
    
    output:
    path "sample_qc_metrics.tsv", emit: metrics
    path "sample_relatedness.tsv", emit: relatedness
    path "sample_heterozygosity.tsv", emit: heterozygosity
    path "sample_depth_stats.tsv", emit: depth_stats
    path "pool_specific_metrics.tsv", emit: pool_metrics
    path "sample_missingness.tsv", emit: missingness
    
    script:
    def has_ploidy_map = ploidy_map.name != 'NO_FILE'
    """
    echo "Calculating sample-level QC metrics"
    
    # Get sample list
    bcftools query -l ${vcf} > samples.txt
    
    # Calculate per-sample missingness and depth statistics
    echo "sample\ttotal_sites\tmissing_sites\tmissing_rate\tmean_depth\tmedian_depth\tmin_depth\tmax_depth\tsd_depth" > sample_depth_stats.tsv
    
    while read sample; do
        echo "Processing sample: \$sample"
        
        # Extract depth values for this sample
        bcftools query -s \$sample -f '[%DP\\n]' ${vcf} > \${sample}_depths.txt
        
        # Calculate statistics using awk
        awk -v sample="\$sample" '
        BEGIN {
            count = 0
            missing = 0
            sum = 0
            sum_sq = 0
            min = 999999
            max = 0
        }
        {
            if(\$1 == "." || \$1 == "") {
                missing++
            } else {
                depths[count++] = \$1
                sum += \$1
                sum_sq += \$1 * \$1
                if(\$1 < min) min = \$1
                if(\$1 > max) max = \$1
            }
        }
        END {
            total = count + missing
            miss_rate = missing / total
            
            if(count > 0) {
                mean = sum / count
                variance = (sum_sq / count) - (mean * mean)
                sd = sqrt(variance)
                
                # Sort for median
                asort(depths)
                if(count % 2 == 0) {
                    median = (depths[count/2] + depths[count/2 + 1]) / 2
                } else {
                    median = depths[(count+1)/2]
                }
            } else {
                mean = 0
                median = 0
                sd = 0
            }
            
            printf "%s\\t%d\\t%d\\t%.4f\\t%.2f\\t%.1f\\t%d\\t%d\\t%.2f\\n", 
                   sample, total, missing, miss_rate, mean, median, min, max, sd
        }' \${sample}_depths.txt >> sample_depth_stats.tsv
        
        rm \${sample}_depths.txt
    done < samples.txt
    
    # Calculate heterozygosity rates
    echo "sample\ttotal_genotypes\thomozygous_ref\thomozygous_alt\theterozygous\theterozygosity_rate\tinbreeding_coefficient" > sample_heterozygosity.tsv
    
    while read sample; do
        # Extract genotypes
        bcftools query -s \$sample -f '[%GT\\n]' ${vcf} > \${sample}_genotypes.txt
        
        # Count genotype classes
        total=\$(wc -l < \${sample}_genotypes.txt)
        hom_ref=\$(grep -c "^0/0\$\\|^0|0\$" \${sample}_genotypes.txt || true)
        het=\$(grep -c "^0/1\$\\|^1/0\$\\|^0|1\$\\|^1|0\$" \${sample}_genotypes.txt || true)
        hom_alt=\$(grep -c "^1/1\$\\|^1|1\$" \${sample}_genotypes.txt || true)
        
        # Calculate heterozygosity rate and F (inbreeding coefficient)
        if [ \$((hom_ref + het + hom_alt)) -gt 0 ]; then
            het_rate=\$(echo "scale=6; \$het / (\$hom_ref + \$het + \$hom_alt)" | bc)
            
            # Expected heterozygosity under HWE (rough estimate)
            p=\$(echo "scale=6; (2*\$hom_ref + \$het) / (2*(\$hom_ref + \$het + \$hom_alt))" | bc)
            q=\$(echo "scale=6; 1 - \$p" | bc)
            exp_het=\$(echo "scale=6; 2 * \$p * \$q" | bc)
            
            # Inbreeding coefficient
            if (( \$(echo "\$exp_het > 0" | bc -l) )); then
                f_coeff=\$(echo "scale=6; 1 - (\$het_rate / \$exp_het)" | bc)
            else
                f_coeff=0
            fi
        else
            het_rate=0
            f_coeff=0
        fi
        
        echo -e "\$sample\\t\$total\\t\$hom_ref\\t\$hom_alt\\t\$het\\t\$het_rate\\t\$f_coeff" >> sample_heterozygosity.tsv
        
        rm \${sample}_genotypes.txt
    done < samples.txt
    
    # Calculate sample missingness patterns
    echo "sample\ttotal_variants\tcalled_variants\tmissing_variants\tmissingness_rate\tsingleton_count\tdoubleton_count" > sample_missingness.tsv
    
    while read sample; do
        total=\$(bcftools view -H ${vcf} | wc -l)
        
        # Count missing genotypes
        missing=\$(bcftools query -s \$sample -f '[%GT\\n]' ${vcf} | grep -c "^\\./\\.\$\\|^\\.|\\.\$" || true)
        called=\$((total - missing))
        miss_rate=\$(echo "scale=6; \$missing / \$total" | bc)
        
        # Count singletons (private variants) for this sample
        singletons=\$(bcftools query -s \$sample -f '[%GT]\\t%AC\\n' ${vcf} | \
                     awk '\$1 != "./." && \$1 != "0/0" && \$2 == "1" {count++} END {print count+0}')
        
        # Count doubletons
        doubletons=\$(bcftools query -s \$sample -f '[%GT]\\t%AC\\n' ${vcf} | \
                     awk '\$1 != "./." && \$1 != "0/0" && \$2 == "2" {count++} END {print count+0}')
        
        echo -e "\$sample\\t\$total\\t\$called\\t\$missing\\t\$miss_rate\\t\$singletons\\t\$doubletons" >> sample_missingness.tsv
    done < samples.txt
    
    # Calculate pairwise relatedness (using method of moments for simplicity)
    echo "sample1\tsample2\tshared_variants\tibs0\tibs1\tibs2\tkinship_coefficient" > sample_relatedness.tsv
    
    # Get all pairwise combinations
    samples=(\$(cat samples.txt))
    n_samples=\${#samples[@]}
    
    if [ \$n_samples -gt 1 ]; then
        for ((i=0; i<\$n_samples-1; i++)); do
            for ((j=i+1; j<\$n_samples; j++)); do
                sample1=\${samples[\$i]}
                sample2=\${samples[\$j]}
                
                # Calculate IBS statistics
                bcftools query -s \$sample1,\$sample2 -f '[%GT\\t]\\n' ${vcf} | \
                awk -v s1="\$sample1" -v s2="\$sample2" '
                BEGIN {
                    ibs0 = 0; ibs1 = 0; ibs2 = 0; shared = 0
                }
                {
                    # Parse genotypes
                    gt1 = \$1; gt2 = \$2
                    
                    # Skip if either is missing
                    if(gt1 == "./." || gt2 == "./.") next
                    
                    shared++
                    
                    # Count alleles
                    if(gt1 == "0/0" || gt1 == "0|0") a1 = 0
                    else if(gt1 == "0/1" || gt1 == "0|1" || gt1 == "1/0" || gt1 == "1|0") a1 = 1
                    else if(gt1 == "1/1" || gt1 == "1|1") a1 = 2
                    else next
                    
                    if(gt2 == "0/0" || gt2 == "0|0") a2 = 0
                    else if(gt2 == "0/1" || gt2 == "0|1" || gt2 == "1/0" || gt2 == "1|0") a2 = 1
                    else if(gt2 == "1/1" || gt2 == "1|1") a2 = 2
                    else next
                    
                    # Calculate IBS
                    diff = a1 - a2
                    if(diff < 0) diff = -diff
                    
                    if(diff == 0) ibs2++
                    else if(diff == 1) ibs1++
                    else if(diff == 2) ibs0++
                }
                END {
                    if(shared > 0) {
                        # Estimate kinship coefficient (simplified)
                        kinship = (ibs1 * 0.25 + ibs2 * 0.5) / shared
                        printf "%s\\t%s\\t%d\\t%d\\t%d\\t%d\\t%.4f\\n", 
                               s1, s2, shared, ibs0, ibs1, ibs2, kinship
                    }
                }' >> sample_relatedness.tsv
            done
        done
    fi
    
    # Pool-specific metrics if ploidy map is provided
    if [ "${has_ploidy_map}" = "true" ]; then
        echo "Calculating pool-specific metrics..."
        echo "sample\tpool_size\testimated_individuals\teffective_alleles\tminor_allele_freq_mean\tminor_allele_freq_median\tpolymorphic_sites" > pool_specific_metrics.tsv
        
        # Parse ploidy map
        while read line; do
            [[ "\$line" =~ ^#.*\$ ]] && continue
            [[ -z "\$line" ]] && continue
            
            bam_name=\$(echo "\$line" | awk '{print \$1}')
            ploidy=\$(echo "\$line" | awk '{print \$2}')
            pool_size=\$((ploidy / 2))  # Assuming diploid
            
            # Get sample name from BAM
            sample_name="\${bam_name%.bam}"
            sample_name="\${sample_name%.filtered}"
            
            # Check if this sample exists in VCF
            if bcftools query -l ${vcf} | grep -q "^\$sample_name\$"; then
                # Calculate pool-specific metrics
                bcftools query -s \$sample_name -f '%AF\\n' -i 'GT!="./." && GT!="0/0"' ${vcf} > \${sample_name}_afs.txt
                
                if [ -s \${sample_name}_afs.txt ]; then
                    polymorphic=\$(wc -l < \${sample_name}_afs.txt)
                    
                    # Calculate AF statistics
                    af_stats=\$(awk '
                    BEGIN {sum=0; count=0}
                    {
                        af = \$1
                        if(af > 0.5) af = 1 - af  # Convert to minor allele freq
                        afs[count++] = af
                        sum += af
                    }
                    END {
                        if(count > 0) {
                            mean = sum / count
                            asort(afs)
                            if(count % 2 == 0)
                                median = (afs[count/2] + afs[count/2 + 1]) / 2
                            else
                                median = afs[(count+1)/2]
                        } else {
                            mean = 0
                            median = 0
                        }
                        printf "%.4f\\t%.4f", mean, median
                    }' \${sample_name}_afs.txt)
                    
                    maf_mean=\$(echo "\$af_stats" | cut -f1)
                    maf_median=\$(echo "\$af_stats" | cut -f2)
                    
                    # Estimate effective number of alleles
                    effective_alleles=\$(echo "scale=2; \$pool_size * 2 * \$maf_mean" | bc)
                else
                    polymorphic=0
                    maf_mean=0
                    maf_median=0
                    effective_alleles=0
                fi
                
                echo -e "\$sample_name\\t\$pool_size\\t\$pool_size\\t\$effective_alleles\\t\$maf_mean\\t\$maf_median\\t\$polymorphic" >> pool_specific_metrics.tsv
                
                rm -f \${sample_name}_afs.txt
            fi
        done < ${ploidy_map}
    else
        # Create empty file if no ploidy map
        echo "sample\tpool_size\testimated_individuals\teffective_alleles\tminor_allele_freq_mean\tminor_allele_freq_median\tpolymorphic_sites" > pool_specific_metrics.tsv
    fi
    
    echo "Sample QC metrics calculation complete"
    """
}