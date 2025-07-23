// modules/freebayes_chunk.nf - Fixed BED file format

process FREEBAYES_CHUNK {
    tag "chunk_${chunk_id}"
    
    input:
    path reference
    path bams
    path bam_indices
    tuple val(chunk_id), val(regions_string)
    
    output:
    tuple val(chunk_id), path("chunk_${chunk_id}.vcf.gz")
    
    script:
    def bam_args = bams.collect{ "--bam $it" }.join(' ')
    """
    # Ensure BAM indices exist
    for bam in ${bams.join(' ')}; do
        if [ ! -f \${bam}.bai ]; then
            echo "Creating index for \${bam}"
            samtools index \${bam}
        fi
    done
    
    echo "Processing chunk ${chunk_id}"
    echo "Regions: ${regions_string}"
    
    # Parse regions string and create BED file for this chunk
    echo "${regions_string}" | tr ',' '\\n' > regions_list.txt
    
    # Convert regions to BED format for freebayes --targets
    # Create BED file WITHOUT comment lines (freebayes doesn't like them)
    touch chunk_targets.bed
    
    while read region; do
        if [ ! -z "\$region" ]; then
            # Parse region format: chr:start-end
            chrom=\$(echo \$region | cut -d':' -f1)
            positions=\$(echo \$region | cut -d':' -f2)
            start=\$(echo \$positions | cut -d'-' -f1)
            end=\$(echo \$positions | cut -d'-' -f2)
            
            # Convert to 0-based BED format
            bed_start=\$((start - 1))
            
            echo -e "\$chrom\\t\$bed_start\\t\$end" >> chunk_targets.bed
            echo "  Added region: \$chrom:\$start-\$end (BED: \$chrom:\$bed_start-\$end)"
        fi
    done < regions_list.txt
    
    echo "Created BED file for chunk ${chunk_id}:"
    cat chunk_targets.bed
    
    # Verify BED file is not empty
    if [ ! -s chunk_targets.bed ]; then
        echo "ERROR: BED file is empty for chunk ${chunk_id}"
        exit 1
    fi
    
    # Run freebayes with targets file
    echo "Running freebayes on chunk ${chunk_id}..."
    
    freebayes \\
        --fasta-reference ${reference} \\
        --targets chunk_targets.bed \\
        ${bam_args} \\
        --vcf chunk_${chunk_id}.vcf
    
    # Check if freebayes produced output
    if [ -s chunk_${chunk_id}.vcf ]; then
        echo "Freebayes completed successfully for chunk ${chunk_id}"
        
        # Count variants found
        variant_count=\$(grep -v '^#' chunk_${chunk_id}.vcf | wc -l)
        echo "Found \$variant_count variants in chunk ${chunk_id}"
        
        # Compress and index the VCF
        bgzip -c chunk_${chunk_id}.vcf > chunk_${chunk_id}.vcf.gz
        tabix -p vcf chunk_${chunk_id}.vcf.gz
        
        echo "Successfully compressed and indexed chunk ${chunk_id}"
    else
        echo "No variants found in chunk ${chunk_id}, creating empty VCF"
        
        # Create minimal VCF header for empty chunks
        cat > chunk_${chunk_id}.vcf << 'VCF_HEADER'
##fileformat=VCFv4.2
##reference=${reference}
##source=freebayes
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
VCF_HEADER
        
        bgzip -c chunk_${chunk_id}.vcf > chunk_${chunk_id}.vcf.gz
        tabix -p vcf chunk_${chunk_id}.vcf.gz
    fi
    
    # Clean up intermediate files
    rm -f chunk_${chunk_id}.vcf regions_list.txt chunk_targets.bed
    
    echo "Chunk ${chunk_id} processing complete"
    """
}