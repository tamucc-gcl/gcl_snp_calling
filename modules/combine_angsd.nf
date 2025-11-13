// modules/combine_angsd.nf - Combine ANGSD chunk outputs

process COMBINE_ANGSD {
    tag "combining_angsd_outputs"
    publishDir params.output_dir, mode: 'copy'
    
    input:
    path chunk_files  // All chunk output files
    val output_prefix
    val output_format
    
    output:
    path "${output_prefix}.beagle.gz", emit: beagle, optional: true
    path "${output_prefix}.mafs.gz", emit: mafs, optional: true
    path "${output_prefix}.vcf.gz", emit: vcf, optional: true
    path "${output_prefix}.vcf.gz.tbi", emit: vcf_index, optional: true
    path "${output_prefix}.sites", emit: sites
    path "${output_prefix}_summary.txt", emit: summary
    
    script:
    """
    echo "Combining ANGSD outputs"
    echo "Output format: ${output_format}"
    echo "Output prefix: ${output_prefix}"
    
    # Sort chunk files by number
    for file in ${chunk_files.join(' ')}; do
        if [[ \$file == *".beagle.gz" ]]; then
            # Extract chunk number from filename
            chunk_num=\$(echo \$file | sed 's/.*chunk_\\([0-9]*\\).*/\\1/')
            printf "%05d\\t%s\\n" "\$chunk_num" "\$file"
        fi
    done | sort -n > beagle_files.txt
    
    # Check if we have any beagle files
    if [ -s beagle_files.txt ]; then
        echo "Found beagle files to combine:"
        cat beagle_files.txt
        
        # Extract just filenames
        cut -f2 beagle_files.txt > beagle_list.txt
        
        # Combine beagle files (genotype likelihoods)
        echo "Combining genotype likelihood files..."
        
        # Get header from first file
        first_file=\$(head -n1 beagle_list.txt)
        zcat \$first_file | head -n1 > ${output_prefix}.beagle
        
        # Combine all data (skip headers)
        while read file; do
            zcat \$file | tail -n +2 >> ${output_prefix}.beagle
        done < beagle_list.txt
        
        # Compress final beagle file
        bgzip ${output_prefix}.beagle
        
        # Count total sites
        TOTAL_SITES=\$(zcat ${output_prefix}.beagle.gz | tail -n +2 | wc -l)
        echo "Total sites with genotype likelihoods: \$TOTAL_SITES"
    else
        echo "No beagle files found to combine"
        touch ${output_prefix}.beagle.gz
    fi
    
    # Combine MAF files if they exist
    echo "Checking for MAF files..."
    ls chunk_*.mafs.gz 2>/dev/null | sort -V > maf_files.txt || true
    
    if [ -s maf_files.txt ]; then
        echo "Combining MAF (allele frequency) files..."
        
        # Get header
        first_maf=\$(head -n1 maf_files.txt)
        zcat \$first_maf | head -n1 > ${output_prefix}.mafs
        
        # Combine data
        while read file; do
            zcat \$file | tail -n +2 >> ${output_prefix}.mafs
        done < maf_files.txt
        
        bgzip ${output_prefix}.mafs
        
        TOTAL_MAFS=\$(zcat ${output_prefix}.mafs.gz | tail -n +2 | wc -l)
        echo "Total sites in MAF file: \$TOTAL_MAFS"
    else
        echo "No MAF files found"
        touch ${output_prefix}.mafs.gz
    fi
    
    # Handle VCF output if requested
    if [ "${output_format}" = "vcf" ] || [ "${output_format}" = "all" ]; then
        echo "Processing VCF files..."
        
        # Sort VCF files by chunk number
        for file in chunk_*.vcf.gz; do
            if [ -f "\$file" ]; then
                chunk_num=\$(echo \$file | sed 's/.*chunk_\\([0-9]*\\).*/\\1/')
                printf "%05d\\t%s\\n" "\$chunk_num" "\$file"
            fi
        done | sort -n > vcf_files.txt
        
        if [ -s vcf_files.txt ]; then
            cut -f2 vcf_files.txt > vcf_list.txt
            
            # Combine VCFs using bcftools
            bcftools concat \\
                --file-list vcf_list.txt \\
                --output-type z \\
                --output ${output_prefix}.vcf.gz
            
            # Index VCF
            tabix -p vcf ${output_prefix}.vcf.gz
            
            TOTAL_VARIANTS=\$(zcat ${output_prefix}.vcf.gz | grep -v '^#' | wc -l)
            echo "Total variants in VCF: \$TOTAL_VARIANTS"
        else
            echo "No VCF files found"
            # Create empty VCF with header
            echo '##fileformat=VCFv4.2' > empty.vcf
            echo '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' >> empty.vcf
            bgzip empty.vcf
            mv empty.vcf.gz ${output_prefix}.vcf.gz
            tabix -p vcf ${output_prefix}.vcf.gz
        fi
    fi
    
    # Create sites file for future ANGSD runs
    echo "Creating sites file for reference..."
    
    if [ -s ${output_prefix}.beagle.gz ]; then
        # Extract sites from beagle file
        # Format: marker allele1 allele2 Ind0 Ind0 Ind0 Ind1 Ind1 Ind1...
        zcat ${output_prefix}.beagle.gz | tail -n +2 | \\
        awk '{
            split(\$1, a, "_")
            print a[1], a[2]
        }' | sort -k1,1 -k2,2n > ${output_prefix}.sites
        
        TOTAL_SITES_FILE=\$(wc -l < ${output_prefix}.sites)
        echo "Sites file created with \$TOTAL_SITES_FILE positions"
    elif [ -s ${output_prefix}.mafs.gz ]; then
        # Extract from MAF file if beagle not available
        zcat ${output_prefix}.mafs.gz | tail -n +2 | \\
        awk '{print \$1, \$2}' | sort -k1,1 -k2,2n > ${output_prefix}.sites
    else
        touch ${output_prefix}.sites
    fi
    
    # Create summary file
    cat > ${output_prefix}_summary.txt << EOF
    ANGSD Pipeline Summary
    ======================
    Output prefix: ${output_prefix}
    Timestamp: \$(date)
    
    Files Generated:
    ----------------
    EOF
    
    if [ -s ${output_prefix}.beagle.gz ]; then
        echo "- Genotype likelihoods: ${output_prefix}.beagle.gz (\$(zcat ${output_prefix}.beagle.gz | tail -n +2 | wc -l) sites)" >> ${output_prefix}_summary.txt
    fi
    
    if [ -s ${output_prefix}.mafs.gz ]; then
        echo "- Allele frequencies: ${output_prefix}.mafs.gz (\$(zcat ${output_prefix}.mafs.gz | tail -n +2 | wc -l) sites)" >> ${output_prefix}_summary.txt
    fi
    
    if [ -f ${output_prefix}.vcf.gz ]; then
        echo "- VCF file: ${output_prefix}.vcf.gz (\$(zcat ${output_prefix}.vcf.gz | grep -v '^#' | wc -l) variants)" >> ${output_prefix}_summary.txt
    fi
    
    if [ -s ${output_prefix}.sites ]; then
        echo "- Sites file: ${output_prefix}.sites (\$(wc -l < ${output_prefix}.sites) positions)" >> ${output_prefix}_summary.txt
    fi
    
    echo "" >> ${output_prefix}_summary.txt
    echo "Processing complete!" >> ${output_prefix}_summary.txt
    
    cat ${output_prefix}_summary.txt
    """
}