// modules/combine_vcfs.nf - Updated for new chunk naming

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
    
    # Create list of VCF files sorted by chunk number
    # New filename format: chunk_0.vcf.gz, chunk_1.vcf.gz, etc.
    for vcf in ${vcf_files.join(' ')}; do
        # Get the chunk number from filename (format: chunk_N.vcf.gz)
        base=\$(basename \$vcf .vcf.gz)
        chunk_num=\$(echo \$base | sed 's/chunk_//')
        
        # Pad with zeros for proper sorting
        printf "%05d\\t%s\\n" "\$chunk_num" "\$vcf"
    done | sort -n > vcf_sort_list.txt
    
    # Extract just the VCF filenames in sorted order
    cut -f2 vcf_sort_list.txt > vcf_list.txt
    
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
        cat > empty.vcf << 'VCF_HEADER'
##fileformat=VCFv4.2
##reference=${params.reference}
##source=freebayes-parallel
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
VCF_HEADER
        bgzip -c empty.vcf > ${output_name}
        tabix -p vcf ${output_name}
    else
        echo "Found \${#non_empty_vcfs[@]} non-empty VCF files out of ${vcf_files.size()} total"
        
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
    echo "Chunks processed: ${vcf_files.size()}"
    
    # Print first few variants for verification
    echo "=== FIRST FEW VARIANTS ==="
    zcat ${output_name} | grep -v '^#' | head -3
    """
}