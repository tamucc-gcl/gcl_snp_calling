// modules/create_chunks.nf - Genome-based chunking (exact chunk count)

process CREATE_CHUNKS {
    tag "chunking_genome"
    
    input:
    path reference
    val num_chunks
    
    output:
    path "chunks.bed"
    path "chunk_regions.txt"  // Maps chunk_id to list of regions for freebayes
    
    script:
    """
    # Create fasta index if it doesn't exist
    if [ ! -f ${reference}.fai ]; then
        samtools faidx ${reference}
    fi
    
    # Genome-based chunking - divides genome into exact number of chunks
    cat > create_chunks.py << 'PYTHON_SCRIPT'
import sys

num_chunks = ${num_chunks}
chromosomes = []
total_length = 0

# Read all contigs and build genome map
with open('${reference}.fai', 'r') as f:
    for line in f:
        parts = line.strip().split('\\t')
        chrom = parts[0]
        length = int(parts[1])
        chromosomes.append((chrom, length))
        total_length += length

print(f"Total genome: {total_length:,} bp across {len(chromosomes)} contigs")
print(f"Creating exactly {num_chunks} chunks")

# Calculate chunk size
chunk_size = total_length // num_chunks
remainder = total_length % num_chunks
print(f"Base chunk size: {chunk_size:,} bp")
print(f"Remainder: {remainder:,} bp (will be distributed across first {remainder} chunks)")

# Create cumulative position map
cumulative_pos = []
current_pos = 0
for chrom, length in chromosomes:
    cumulative_pos.append((chrom, current_pos, current_pos + length - 1))
    current_pos += length

print(f"Genome position map created")

# Function to find contig and position given global position
def find_contig_pos(global_pos):
    for chrom, start, end in cumulative_pos:
        if start <= global_pos <= end:
            local_pos = global_pos - start + 1  # 1-based position within contig
            return chrom, local_pos
    return None, None

# Create chunks
chunks = []
chunk_regions = {}  # chunk_id -> list of regions

for chunk_id in range(num_chunks):
    # Calculate start and end positions for this chunk
    chunk_start = chunk_id * chunk_size + min(chunk_id, remainder)
    
    if chunk_id < remainder:
        current_chunk_size = chunk_size + 1
    else:
        current_chunk_size = chunk_size
    
    chunk_end = chunk_start + current_chunk_size - 1
    
    print(f"\\nChunk {chunk_id}: global positions {chunk_start:,} - {chunk_end:,} ({current_chunk_size:,} bp)")
    
    # Find all contigs that overlap with this chunk
    chunk_regions_list = []
    
    current_pos = chunk_start
    while current_pos <= chunk_end:
        # Find which contig this position is in
        contig_found = False
        for chrom, contig_start, contig_end in cumulative_pos:
            if contig_start <= current_pos <= contig_end:
                # This contig overlaps with our chunk
                
                # Calculate the overlap region
                region_start_global = max(current_pos, contig_start)
                region_end_global = min(chunk_end, contig_end)
                
                # Convert to local coordinates (1-based)
                region_start_local = region_start_global - contig_start + 1
                region_end_local = region_end_global - contig_start + 1
                
                # Add BED entry (0-based start, 1-based end)
                bed_start = region_start_local - 1
                bed_end = region_end_local
                chunks.append(f"{chrom}\\t{bed_start}\\t{bed_end}")
                
                # Add region for freebayes (1-based)
                region = f"{chrom}:{region_start_local}-{region_end_local}"
                chunk_regions_list.append(region)
                
                region_length = region_end_local - region_start_local + 1
                print(f"  {region} ({region_length:,} bp)")
                
                # Move to next position
                current_pos = contig_end + 1
                contig_found = True
                break
        
        if not contig_found:
            print(f"ERROR: Could not find contig for position {current_pos}")
            break
    
    chunk_regions[f"chunk_{chunk_id}"] = chunk_regions_list
    
    total_chunk_bp = sum(
        int(region.split(':')[1].split('-')[1]) - int(region.split(':')[1].split('-')[0]) + 1 
        for region in chunk_regions_list
    )
    print(f"  Total: {len(chunk_regions_list)} regions, {total_chunk_bp:,} bp")

# Write BED file
with open('chunks.bed', 'w') as f:
    for chunk in chunks:
        f.write(f"{chunk}\\n")

# Write chunk regions file for freebayes
with open('chunk_regions.txt', 'w') as f:
    for chunk_id, regions in chunk_regions.items():
        f.write(f"{chunk_id}\\t{','.join(regions)}\\n")

print(f"\\nSummary:")
print(f"- Created exactly {num_chunks} logical chunks")
print(f"- Total BED regions: {len(chunks)}")
print(f"- Chunk regions file created for freebayes processing")

# Verify coverage
total_covered = 0
for chunk in chunks:
    parts = chunk.split('\\t')
    start, end = int(parts[1]), int(parts[2])
    total_covered += (end - start)

print(f"- Total coverage: {total_covered:,} bp ({100*total_covered/total_length:.1f}% of genome)")

if total_covered == total_length:
    print("✅ SUCCESS: Complete genome coverage with exact chunk count!")
else:
    print(f"⚠️  Coverage difference: {total_length - total_covered:,} bp")
PYTHON_SCRIPT

    python3 create_chunks.py
    """
}