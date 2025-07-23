// modules/create_chunks.nf

process CREATE_CHUNKS {
    tag "chunking_genome"
    
    input:
    path reference
    val num_chunks
    
    output:
    path "chunks.bed"
    
    script:
    """
    # Create fasta index if it doesn't exist
    if [ ! -f ${reference}.fai ]; then
        samtools faidx ${reference}
    fi
    
    # Create chunks using a python script
    python3 << 'EOF'
    import sys
    import math

    num_chunks = ${num_chunks}
    chromosomes = []
    total_length = 0

    # Read the fasta index file to get chromosome info
    with open('${reference}.fai', 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            chrom = parts[0]
            length = int(parts[1])
            chromosomes.append((chrom, length))
            total_length += length

    print(f"Total genome length: {total_length:,} bp")
    print(f"Target number of chunks: {num_chunks}")

    # Calculate target chunk size
    target_chunk_size = total_length // num_chunks
    print(f"Target chunk size: {target_chunk_size:,} bp")

    chunks = []

    # Create chunks for each chromosome
    for chrom, chrom_length in chromosomes:
        # Calculate how many chunks this chromosome should contribute
        chrom_chunks = max(1, round(chrom_length / target_chunk_size))
        actual_chunk_size = chrom_length // chrom_chunks
        
        print(f"Chromosome {chrom}: {chrom_length:,} bp -> {chrom_chunks} chunks of ~{actual_chunk_size:,} bp each")
        
        # Create chunks for this chromosome
        start = 1
        for i in range(chrom_chunks):
            if i == chrom_chunks - 1:  # Last chunk gets remainder
                end = chrom_length
            else:
                end = start + actual_chunk_size - 1
            
            chunks.append(f"{chrom}\t{start-1}\t{end}")
            start = end + 1

    # Write chunks to BED file
    with open('chunks.bed', 'w') as f:
        for chunk in chunks:
            f.write(f"{chunk}\n")

    print(f"Created {len(chunks)} total chunks (target was {num_chunks})")
    EOF
    """
}