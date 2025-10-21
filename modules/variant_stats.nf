// modules/variant_stats.nf - Extract comprehensive variant statistics

process VARIANT_STATS {
    tag "variant_statistics"
    publishDir "${params.output_dir}/qc", mode: 'copy'
    
    input:
    path vcf
    path vcf_index
    
    output:
    path "variant_stats.txt", emit: stats
    path "variant_summary.json", emit: summary
    path "per_sample_stats.tsv", emit: sample_stats
    path "per_chromosome_stats.tsv", emit: chr_stats
    path "allele_frequency_dist.tsv", emit: af_dist
    path "depth_distribution.tsv", emit: depth_dist
    path "quality_distribution.tsv", emit: qual_dist
    
    script:
    """
    echo "Extracting variant statistics from ${vcf}"
    
    # Basic statistics using bcftools
    bcftools stats ${vcf} > variant_stats.txt
    
    # Extract per-sample statistics
    echo "sample\ttotal_variants\thomozygous_ref\thomozygous_alt\theterozygous\tmissing\tmissing_rate" > per_sample_stats.tsv
    bcftools query -l ${vcf} | while read sample; do
        echo "Processing sample: \$sample"
        
        # Count genotype types for this sample
        total=\$(bcftools view -H ${vcf} | wc -l)
        hom_ref=\$(bcftools query -s \$sample -f '[%GT\\n]' ${vcf} | grep -c "^0/0\$" || true)
        hom_alt=\$(bcftools query -s \$sample -f '[%GT\\n]' ${vcf} | grep -c "^1/1\$" || true)
        het=\$(bcftools query -s \$sample -f '[%GT\\n]' ${vcf} | grep -E -c "^0/1\$|^1/0\$" || true)
        missing=\$(bcftools query -s \$sample -f '[%GT\\n]' ${vcf} | grep -c "^\\./\\.\$" || true)
        missing_rate=\$(echo "scale=4; \$missing / \$total" | bc)
        
        echo -e "\$sample\\t\$total\\t\$hom_ref\\t\$hom_alt\\t\$het\\t\$missing\\t\$missing_rate" >> per_sample_stats.tsv
    done
    
    # Extract per-chromosome statistics
    echo "chromosome\ttotal_variants\tsnps\tindels\tother\tavg_depth\tavg_qual" > per_chromosome_stats.tsv
    bcftools index -s ${vcf} | cut -f1 | while read chr; do
        echo "Processing chromosome: \$chr"
        
        total=\$(bcftools view -H -r \$chr ${vcf} | wc -l)
        snps=\$(bcftools view -H -r \$chr -v snps ${vcf} | wc -l)
        indels=\$(bcftools view -H -r \$chr -v indels ${vcf} | wc -l)
        other=\$((total - snps - indels))
        
        # Calculate average depth and quality
        avg_depth=\$(bcftools query -r \$chr -f '%INFO/DP\\n' ${vcf} | \
                    awk '{sum+=\$1; count++} END {if(count>0) printf "%.2f", sum/count; else print "0"}')
        avg_qual=\$(bcftools query -r \$chr -f '%QUAL\\n' ${vcf} | \
                   awk '{sum+=\$1; count++} END {if(count>0) printf "%.2f", sum/count; else print "0"}')
        
        echo -e "\$chr\\t\$total\\t\$snps\\t\$indels\\t\$other\\t\$avg_depth\\t\$avg_qual" >> per_chromosome_stats.tsv
    done
    
    # Extract allele frequency distribution (for all variants)
    echo "af_bin\tcount\tvariant_type" > allele_frequency_dist.tsv
    
    # For biallelic SNPs
    bcftools query -f '%AF\\n' -i 'TYPE="snp" && N_ALT=1' ${vcf} | \
    awk '{
        if(\$1 < 0.01) bin="0-0.01"
        else if(\$1 < 0.05) bin="0.01-0.05"
        else if(\$1 < 0.1) bin="0.05-0.1"
        else if(\$1 < 0.25) bin="0.1-0.25"
        else if(\$1 < 0.5) bin="0.25-0.5"
        else if(\$1 < 0.75) bin="0.5-0.75"
        else if(\$1 < 0.95) bin="0.75-0.95"
        else bin="0.95-1.0"
        count[bin]++
    }
    END {
        for(b in count) print b"\t"count[b]"\tSNP"
    }' | sort >> allele_frequency_dist.tsv
    
    # For INDELs
    bcftools query -f '%AF\\n' -i 'TYPE="indel"' ${vcf} | \
    awk '{
        if(\$1 < 0.01) bin="0-0.01"
        else if(\$1 < 0.05) bin="0.01-0.05"
        else if(\$1 < 0.1) bin="0.05-0.1"
        else if(\$1 < 0.25) bin="0.1-0.25"
        else if(\$1 < 0.5) bin="0.25-0.5"
        else if(\$1 < 0.75) bin="0.5-0.75"
        else if(\$1 < 0.95) bin="0.75-0.95"
        else bin="0.95-1.0"
        count[bin]++
    }
    END {
        for(b in count) print b"\t"count[b]"\tINDEL"
    }' | sort >> allele_frequency_dist.tsv
    
    # Extract depth distribution
    echo "depth_bin\tcount" > depth_distribution.tsv
    bcftools query -f '%INFO/DP\\n' ${vcf} | \
    awk '{
        if(\$1 < 5) bin="0-5"
        else if(\$1 < 10) bin="5-10"
        else if(\$1 < 20) bin="10-20"
        else if(\$1 < 30) bin="20-30"
        else if(\$1 < 50) bin="30-50"
        else if(\$1 < 100) bin="50-100"
        else if(\$1 < 200) bin="100-200"
        else bin="200+"
        count[bin]++
    }
    END {
        for(b in count) print b"\t"count[b]
    }' | sort >> depth_distribution.tsv
    
    # Extract quality distribution
    echo "quality_bin\tcount" > quality_distribution.tsv
    bcftools query -f '%QUAL\\n' ${vcf} | \
    awk '{
        if(\$1 < 10) bin="0-10"
        else if(\$1 < 20) bin="10-20"
        else if(\$1 < 30) bin="20-30"
        else if(\$1 < 50) bin="30-50"
        else if(\$1 < 100) bin="50-100"
        else if(\$1 < 200) bin="100-200"
        else if(\$1 < 500) bin="200-500"
        else bin="500+"
        count[bin]++
    }
    END {
        for(b in count) print b"\t"count[b]
    }' | sort >> quality_distribution.tsv
    
    # Create JSON summary for easy parsing
    cat > variant_summary.json << EOF
{
    "total_variants": \$(bcftools view -H ${vcf} | wc -l),
    "snps": \$(bcftools view -H -v snps ${vcf} | wc -l),
    "indels": \$(bcftools view -H -v indels ${vcf} | wc -l),
    "multiallelic": \$(bcftools view -H -m3 ${vcf} | wc -l),
    "samples": \$(bcftools query -l ${vcf} | wc -l),
    "chromosomes": \$(bcftools index -s ${vcf} | wc -l),
    "ti_tv_ratio": \$(bcftools stats ${vcf} | grep "TSTV" | tail -1 | cut -f5)
}
EOF
    
    echo "Variant statistics extraction complete"
    """
}
