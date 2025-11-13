// modules/combine_angsd.nf - Combine ANGSD chunk outputs

process COMBINE_ANGSD {
    tag "combining_angsd_outputs"
    publishDir params.output_dir, mode: 'copy'
    
    input:
    path chunk_files  // All chunk output files
    val output_prefix
    
    output:
    path "${output_prefix}.beagle.gz", emit: beagle
    path "${output_prefix}.mafs.gz", emit: mafs
    path "${output_prefix}.vcf.gz", emit: vcf
    path "${output_prefix}.vcf.gz.tbi", emit: vcf_index
    path "${output_prefix}.sites", emit: sites
    path "${output_prefix}_summary.txt", emit: summary
    
    script:
    """
    echo "Combining ANGSD outputs"
    echo "Output prefix: ${output_prefix}"
    echo "Processing both Beagle and BCF/VCF formats"
    
    # Sort chunk files by number for Beagle files
    for file in ${chunk_files.join(' ')}; do
        if [[ \$file == *".beagle.gz" ]]; then
            # Extract chunk number from filename
            chunk_num=\$(echo \$file | sed 's/.*chunk_\\([0-9]*\\).*/\\1/')
            printf "%05d\\t%s\\n" "\$chunk_num" "\$file"
        fi
    done | sort -n > beagle_files.txt
    
    # Combine beagle files (genotype likelihoods)
    if [ -s beagle_files.txt ]; then
        echo "Found beagle files to combine:"
        cat beagle_files.txt
        
        # Extract just filenames
        cut -f2 beagle_files.txt > beagle_list.txt
        
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
        
        TOTAL_SITES=\$(zcat ${output_prefix}.beagle.gz | tail -n +2 | wc -l)
        echo "Total sites with genotype likelihoods: \$TOTAL_SITES"
    else
        echo "ERROR: No beagle files found to combine"
        exit 1
    fi
    
    # Combine MAF files
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
        echo "No MAF files found - creating empty file"
        echo "chromo\tposition\tmajor\tminor\tref\tanc\tknownEM\tnInd" | bgzip > ${output_prefix}.mafs.gz
    fi
    
    # Handle BCF files and convert to VCF
    echo "Processing BCF files..."
    
    # Sort BCF files by chunk number
    for file in chunk_*.bcf; do
        if [ -f "\$file" ]; then
            chunk_num=\$(echo \$file | sed 's/.*chunk_\\([0-9]*\\).*/\\1/')
            printf "%05d\\t%s\\n" "\$chunk_num" "\$file"
        fi
    done | sort -n > bcf_files.txt
    
    if [ -s bcf_files.txt ]; then
        echo "Found BCF files to combine:"
        cat bcf_files.txt
        
        cut -f2 bcf_files.txt > bcf_list.txt
        
        # Combine BCFs using bcftools
        echo "Concatenating BCF files..."
        bcftools concat \\
            --file-list bcf_list.txt \\
            --output-type b \\
            --output combined.bcf
        
        # Convert BCF to VCF
        echo "Converting BCF to compressed VCF..."
        bcftools view \\
            --output-type z \\
            --output ${output_prefix}.vcf.gz \\
            combined.bcf
        
        # Index VCF
        tabix -p vcf ${output_prefix}.vcf.gz
        
        TOTAL_VARIANTS=\$(zcat ${output_prefix}.vcf.gz | grep -v '^#' | wc -l)
        echo "Total variants in VCF: \$TOTAL_VARIANTS"
        
        # Clean up intermediate BCF
        rm -f combined.bcf
    else
        echo "No BCF files found - checking for VCF files as fallback"
        
        # Fallback to VCF files if BCF not found (shouldn't happen with new config)
        ls chunk_*.vcf.gz 2>/dev/null | sort -V > vcf_files.txt || true
        
        if [ -s vcf_files.txt ]; then
            bcftools concat \\
                --file-list vcf_files.txt \\
                --output-type z \\
                --output ${output_prefix}.vcf.gz
            
            tabix -p vcf ${output_prefix}.vcf.gz
        else
            echo "No BCF or VCF files found - creating minimal VCF"
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
    - Genotype likelihoods: ${output_prefix}.beagle.gz (\$(zcat ${output_prefix}.beagle.gz | tail -n +2 | wc -l) sites)
    - Allele frequencies: ${output_prefix}.mafs.gz (\$(zcat ${output_prefix}.mafs.gz | tail -n +2 | wc -l) sites)
    - VCF file: ${output_prefix}.vcf.gz (\$(zcat ${output_prefix}.vcf.gz | grep -v '^#' | wc -l) variants)
    - Sites file: ${output_prefix}.sites (\$(wc -l < ${output_prefix}.sites) positions)
    
    Processing complete!
    EOF
    
    cat ${output_prefix}_summary.txt
    """
}