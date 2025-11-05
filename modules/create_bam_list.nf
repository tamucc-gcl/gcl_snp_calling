// modules/create_bam_list.nf
process CREATE_BAM_LIST {
    tag "creating_bam_list"
    
    input:
    val(bam_pattern)
    
    output:
    path("bam_list.txt")
    
    script:
    """
    # Find all BAM files matching the pattern and sort them deterministically
    # This creates a stable, cacheable file list
    find \$(dirname "${bam_pattern}") -maxdepth 1 -name "\$(basename "${bam_pattern}")" -type f | sort > bam_list.txt
    
    # Verify we found some BAMs
    if [ ! -s bam_list.txt ]; then
        echo "ERROR: No BAM files found matching pattern: ${bam_pattern}" >&2
        exit 1
    fi
    
    # Show what we found
    echo "Found \$(wc -l < bam_list.txt) BAM files:"
    cat bam_list.txt
    
    # Verify indices exist
    while read bam; do
        if [ ! -f "\${bam}.bai" ]; then
            echo "ERROR: Index not found for \${bam}" >&2
            exit 1
        fi
    done < bam_list.txt
    
    echo "All BAM files have indices"
    """
}