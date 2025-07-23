// modules/combine_vcfs.nf

process COMBINE_VCFS {
    tag "combining_vcfs"
    publishDir params.output_dir, mode: 'copy'
    
    input:
    path vcf_files
    val output_name
    
    output:
    path output_name
    path "${output_name}.tbi"
    
    script:
    """
    echo "Combining ${vcf_files.size()} VCF files..."
    
    # Create list of VCF files sorted by genomic position
    # Extract chromosome and start position for sorting
    for vcf in ${vcf_files.join(' ')}; do
        # Get the chunk info from filename (format: chunk_chr_start_end.vcf.gz)
        base=\$(basename \$vcf .vcf.gz)
        # Extract chr, start, end from filename
        chrom=\$(echo \$base | cut -d'_' -f2)
        start=\$(echo \$base | cut -d'_' -f3)
        end=\$(echo \$base | cut -d'_' -f4)
        echo -e "\$chrom\t\$start\t\$end\t\$vcf"
    done | sort -k1,1V -k2,2n > vcf_sort_list.txt
    
    # Extract just the VCF filenames in sorted order
    cut -f4 vcf_sort_list.txt > vcf_list.txt
    
    echo "VCF files will be combined in this order:"
    cat vcf_list.txt
    
    # Check if we have any non-empty VCFs
    non_empty_vcfs=()
    while read vcf; do
        # Check if VCF has variants (more than just header)
        if [ \$(zcat \$vcf | grep -v '^#' | wc -l) -gt 0 ]; then
            non_empty_vcfs+=(\$vcf)
        fi
    done < vcf_list.txt
    
    if [ \${#non_empty_vcfs[@]} -eq 0 ]; then
        echo "No variants found in any chunks. Creating empty VCF..."
        # Create minimal VCF with proper header
        echo "##fileformat=VCFv4.2" > empty.vcf
        echo "##reference=${params.reference}" >> empty.vcf
        echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> empty.vcf
        bgzip -c empty.vcf > ${output_name}
        tabix -p vcf ${output_name}
    else
        echo "Found \${#non_empty_vcfs[@]} non-empty VCF files"
        
        # Create list of non-empty VCFs
        printf '%s\n' "\${non_empty_vcfs[@]}" > non_empty_vcf_list.txt
        
        # Combine VCFs using bcftools
        bcftools concat \\
            --file-list non_empty_vcf_list.txt \\
            --output-type z \\
            --output ${output_name}
        
        # Index the final VCF
        tabix -p vcf ${output_name}
    fi
    
    echo "Final VCF created: ${output_name}"
    
    # Print summary statistics
    echo "=== SUMMARY ==="
    echo "Total variants: \$(zcat ${output_name} | grep -v '^#' | wc -l)"
    echo "File size: \$(ls -lh ${output_name} | awk '{print \$5}')"
    """
}