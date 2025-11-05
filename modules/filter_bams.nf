// modules/filter_bams.nf - Filter BAM files based on configurable parameters

process FILTER_BAMS {
    tag "${bam.baseName}"
    
    input:
    path bam
    path bai
    path filter_config
    
    output:
    tuple path("${bam.baseName}.filtered.bam"), path("${bam.baseName}.filtered.bam.bai")
    
    script:
    def has_config = filter_config.name != 'NO_FILE'
    """
    echo "Filtering BAM file: ${bam}"
    
    # Initialize filtering parameters with defaults
    MAPPING_MIN_QUALITY=20
    REMOVE_UNMAPPED="no"
    REMOVE_UNPAIRED_MAPPED="yes"
    REMOVE_SECONDARY="yes"
    REMOVE_VENDOR_FAILED="no"
    REMOVE_DUPLICATES="no"
    REMOVE_SUPPLEMENTARY="no"
    KEEP_ONLY_PROPER_PAIRS="yes"
    CUSTOM_F_FLAG=0
    CUSTOM_f_FLAG=0
    MAX_SOFT_CLIP=20
    MIN_ALIGN_SCORE=50
    REMOVE_ORPHANS="yes"
    
    # Load parameters from config if provided
    if [ "${has_config}" = "true" ]; then
        echo "Loading filter configuration from ${filter_config}"
        
        # Parse JSON config using Python
        python3 << 'PYTHON_SCRIPT' > filter_params.sh
import json
import sys

try:
    with open('${filter_config}', 'r') as f:
        config = json.load(f)
    
    # Extract filter parameters with defaults
    filters = config.get('filter_parameters', {})
    
    print(f"MAPPING_MIN_QUALITY={filters.get('mapping_min_quality', 20)}")
    print(f"REMOVE_UNMAPPED='{filters.get('remove_unmapped_reads', 'no')}'")
    print(f"REMOVE_UNPAIRED_MAPPED='{filters.get('remove_read_pair_if_one_is_unmapped', 'yes')}'")
    print(f"REMOVE_SECONDARY='{filters.get('remove_secondary_alignments', 'yes')}'")
    print(f"REMOVE_VENDOR_FAILED='{filters.get('remove_reads_not_passing_platform_vendor_filters', 'no')}'")
    print(f"REMOVE_DUPLICATES='{filters.get('remove_pcr_or_optical_duplicates', 'no')}'")
    print(f"REMOVE_SUPPLEMENTARY='{filters.get('remove_supplementary_alignments', 'no')}'")
    print(f"KEEP_ONLY_PROPER_PAIRS='{filters.get('keep_only_properly_aligned_read_pairs', 'yes')}'")
    print(f"CUSTOM_F_FLAG={filters.get('custom_samtools_view_F_bit_value', 0)}")
    print(f"CUSTOM_f_FLAG={filters.get('custom_samtools_view_f_bit_value', 0)}")
    print(f"MAX_SOFT_CLIP={filters.get('remove_reads_with_excessive_soft_clipping', 20)}")
    print(f"MIN_ALIGN_SCORE={filters.get('remove_reads_with_alignment_score_below', 50)}")
    print(f"REMOVE_ORPHANS='{filters.get('remove_reads_orphaned_by_filters', 'yes')}'")
    
except Exception as e:
    print(f"# Error reading config: {e}", file=sys.stderr)
    sys.exit(1)
PYTHON_SCRIPT
        
        # Source the parameters
        source filter_params.sh
    fi
    
    echo "Filter parameters:"
    echo "  Mapping quality threshold: \$MAPPING_MIN_QUALITY"
    echo "  Remove unmapped: \$REMOVE_UNMAPPED"
    echo "  Remove unpaired mapped: \$REMOVE_UNPAIRED_MAPPED"
    echo "  Remove secondary alignments: \$REMOVE_SECONDARY"
    echo "  Remove vendor failed: \$REMOVE_VENDOR_FAILED"
    echo "  Remove duplicates: \$REMOVE_DUPLICATES"
    echo "  Remove supplementary: \$REMOVE_SUPPLEMENTARY"
    echo "  Keep only proper pairs: \$KEEP_ONLY_PROPER_PAIRS"
    echo "  Custom -F flag: \$CUSTOM_F_FLAG"
    echo "  Custom -f flag: \$CUSTOM_f_FLAG"
    echo "  Max soft clipping: \$MAX_SOFT_CLIP"
    echo "  Min alignment score: \$MIN_ALIGN_SCORE"
    echo "  Remove orphans: \$REMOVE_ORPHANS"
    
    # Build samtools filter flags
    F_FLAGS=0
    f_FLAGS=0
    
    # Build -F flags (exclude reads with these flags)
    [ "\$REMOVE_UNMAPPED" = "yes" ] && F_FLAGS=\$((F_FLAGS | 4))
    [ "\$REMOVE_UNPAIRED_MAPPED" = "yes" ] && F_FLAGS=\$((F_FLAGS | 8))
    [ "\$REMOVE_SECONDARY" = "yes" ] && F_FLAGS=\$((F_FLAGS | 256))
    [ "\$REMOVE_VENDOR_FAILED" = "yes" ] && F_FLAGS=\$((F_FLAGS | 512))
    [ "\$REMOVE_DUPLICATES" = "yes" ] && F_FLAGS=\$((F_FLAGS | 1024))
    [ "\$REMOVE_SUPPLEMENTARY" = "yes" ] && F_FLAGS=\$((F_FLAGS | 2048))
    
    # Build -f flags (require reads with these flags)
    [ "\$KEEP_ONLY_PROPER_PAIRS" = "yes" ] && f_FLAGS=\$((f_FLAGS | 2))
    
    echo "Computed samtools flags: -F \$F_FLAGS -f \$f_FLAGS"
    
    # First pass: basic samtools filtering
    SAMTOOLS_CMD="samtools view -@ ${task.cpus} -b"
    SAMTOOLS_CMD="\$SAMTOOLS_CMD -q \$MAPPING_MIN_QUALITY"
    [ \$F_FLAGS -ne 0 ] && SAMTOOLS_CMD="\$SAMTOOLS_CMD -F \$F_FLAGS"
    [ \$f_FLAGS -ne 0 ] && SAMTOOLS_CMD="\$SAMTOOLS_CMD -f \$f_FLAGS"
    SAMTOOLS_CMD="\$SAMTOOLS_CMD ${bam}"
    
    echo "Samtools command: \$SAMTOOLS_CMD"
    
    # Apply basic filters and pipe to temporary file
    \$SAMTOOLS_CMD > temp_filtered.bam
    
    # Apply custom F and f flags if specified
    if [ \$CUSTOM_F_FLAG -ne 0 ] || [ \$CUSTOM_f_FLAG -ne 0 ]; then
        echo "Applying custom flags: -F \$CUSTOM_F_FLAG -f \$CUSTOM_f_FLAG"
        CUSTOM_CMD="samtools view -@ ${task.cpus} -b"
        [ \$CUSTOM_F_FLAG -ne 0 ] && CUSTOM_CMD="\$CUSTOM_CMD -F \$CUSTOM_F_FLAG"
        [ \$CUSTOM_f_FLAG -ne 0 ] && CUSTOM_CMD="\$CUSTOM_CMD -f \$CUSTOM_f_FLAG"
        CUSTOM_CMD="\$CUSTOM_CMD temp_filtered.bam"
        
        \$CUSTOM_CMD > temp_filtered2.bam
        mv temp_filtered2.bam temp_filtered.bam
    fi
    
    # Apply soft clipping and alignment score filters
    if [ "\$MAX_SOFT_CLIP" != "no" ] && [ "\$MAX_SOFT_CLIP" -gt 0 ] || [ \$MIN_ALIGN_SCORE -gt 0 ]; then
        echo "Applying soft clipping (max: \$MAX_SOFT_CLIP) and alignment score (min: \$MIN_ALIGN_SCORE) filters"
        
        samtools view -h temp_filtered.bam | \\
        awk -v max_clip="\$MAX_SOFT_CLIP" -v min_score="\$MIN_ALIGN_SCORE" '
        BEGIN {OFS="\\t"}
        /^@/ {print; next}  # Print header lines
        {
            # Check alignment score (AS tag)
            as_score = 0
            for(i=12; i<=NF; i++) {
                if(\$i ~ /^AS:i:/) {
                    split(\$i, a, ":")
                    as_score = a[3]
                    break
                }
            }
            if(min_score > 0 && as_score < min_score) next
            
            # Check soft clipping in CIGAR
            if(max_clip != "no" && max_clip > 0) {
                cigar = \$6
                # Extract all soft clip values
                total_soft_clip = 0
                while(match(cigar, /[0-9]+S/)) {
                    clip_val = substr(cigar, RSTART, RLENGTH-1)
                    total_soft_clip += clip_val
                    cigar = substr(cigar, RSTART+RLENGTH)
                }
                if(total_soft_clip > max_clip) next
            }
            
            print
        }' | samtools view -b > temp_filtered3.bam
        
        mv temp_filtered3.bam temp_filtered.bam
    fi
    
    # Remove orphaned reads if requested
    if [ "\$REMOVE_ORPHANS" = "yes" ]; then
        echo "Removing orphaned reads..."
        
        # Sort by name for fixmate
        samtools sort -@ ${task.cpus} -n temp_filtered.bam -o temp_namesorted.bam
        
        # Fix mate information and remove orphans
        samtools fixmate -@ ${task.cpus} -r temp_namesorted.bam temp_fixmate.bam
        
        # Sort back by coordinate
        samtools sort -@ ${task.cpus} temp_fixmate.bam -o ${bam.baseName}.filtered.bam
        
        rm temp_namesorted.bam temp_fixmate.bam
    else
        # Just rename the file
        mv temp_filtered.bam ${bam.baseName}.filtered.bam
    fi
    
    # 2) Name-sort (required before fixmate)
    samtools sort -@ ${task.cpus ?: 8} -n -o ${bam.baseName}.nsrt.bam ${bam.baseName}.filtered.bam

    # 3) Fixmate (adds required ms/MC tags for paired-end)
    samtools fixmate -@ ${task.cpus ?: 8} -m ${bam.baseName}.nsrt.bam ${bam.baseName}.fxmt.bam

    # 4) Coordinate-sort for markdup
    samtools sort -@ ${task.cpus ?: 8} -o ${bam.baseName}.csrt.bam ${bam.baseName}.fxmt.bam

    # 5) Mark duplicates (marks; use -r to remove)
    samtools markdup -@ ${task.cpus ?: 8} ${bam.baseName}.csrt.bam ${bam.baseName}.bam

    # Optional: clean up intermediates to save space
    rm -f ${bam.baseName}.filtered.bam ${bam.baseName}.nsrt.bam ${bam.baseName}.fxmt.bam ${bam.baseName}.csrt.bam temp_filtered.bam

    # 6) Index final BAM
    mv ${bam.baseName}.bam ${bam.baseName}.filtered.bam
    samtools index -@ ${task.cpus ?: 8} ${bam.baseName}.filtered.bam
    
    # REMOVED: Statistics that reference input BAM
    # This was breaking cache because it reads ${bam} again
    echo "Filtering complete for ${bam.baseName}"
    """
}