// modules/freebayes_chunk.nf

process FREEBAYES_CHUNK {
    tag "chunk_${chunk_id}"
    
    input:
    path reference
    path bams
    path bam_indices
    tuple val(chunk_id), val(region)
    
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
    
    echo "Running freebayes on region: ${region}"
    echo "Using BAM files: ${bams.join(', ')}"
    
    # Run freebayes on this chunk
    freebayes \\
        --fasta-reference ${reference} \\
        --region ${region} \\
        ${bam_args} \\
        --vcf chunk_${chunk_id}.vcf
    
    # Compress and index the VCF
    if [ -s chunk_${chunk_id}.vcf ]; then
        bgzip -c chunk_${chunk_id}.vcf > chunk_${chunk_id}.vcf.gz
        tabix -p vcf chunk_${chunk_id}.vcf.gz
        echo "Successfully processed chunk ${chunk_id}"
    else
        echo "No variants found in chunk ${chunk_id}, creating empty VCF"
        # Create minimal VCF header for empty chunks
        echo "##fileformat=VCFv4.2" > chunk_${chunk_id}.vcf
        echo "##reference=${reference}" >> chunk_${chunk_id}.vcf
        echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> chunk_${chunk_id}.vcf
        bgzip -c chunk_${chunk_id}.vcf > chunk_${chunk_id}.vcf.gz
        tabix -p vcf chunk_${chunk_id}.vcf.gz
    fi
    
    # Clean up uncompressed VCF
    rm -f chunk_${chunk_id}.vcf
    """
}