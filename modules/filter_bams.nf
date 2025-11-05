// modules/filter_bams.nf - Step 7 FIXED: Properly escape bash variable

process FILTER_BAMS {
    tag "${bam.simpleName}"
    
    input:
    path bam
    path bai
    path filter_config
    
    output:
    tuple path("${bam.simpleName}.filtered.bam"), path("${bam.simpleName}.filtered.bam.bai")
    
    script:
    def has_config = filter_config.name != 'NO_FILE'
    """
    if [ "${has_config}" = "true" ]; then
        cp ${filter_config} filter_config.json
    fi
    
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
    
    if [ "${has_config}" = "true" ]; then
        python3 << 'PYTHON_SCRIPT' > filter_params.sh
import json
import sys

try:
    with open('filter_config.json', 'r') as f:
        config = json.load(f)
    
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
    print(f"# Error: {e}", file=sys.stderr)
    sys.exit(1)
PYTHON_SCRIPT
        
        source filter_params.sh
    fi
    
    F_FLAGS=0
    f_FLAGS=0
    
    [ "\$REMOVE_UNMAPPED" = "yes" ] && F_FLAGS=\$((F_FLAGS | 4))
    [ "\$REMOVE_UNPAIRED_MAPPED" = "yes" ] && F_FLAGS=\$((F_FLAGS | 8))
    [ "\$REMOVE_SECONDARY" = "yes" ] && F_FLAGS=\$((F_FLAGS | 256))
    [ "\$REMOVE_VENDOR_FAILED" = "yes" ] && F_FLAGS=\$((F_FLAGS | 512))
    [ "\$REMOVE_DUPLICATES" = "yes" ] && F_FLAGS=\$((F_FLAGS | 1024))
    [ "\$REMOVE_SUPPLEMENTARY" = "yes" ] && F_FLAGS=\$((F_FLAGS | 2048))
    [ "\$KEEP_ONLY_PROPER_PAIRS" = "yes" ] && f_FLAGS=\$((f_FLAGS | 2))
    
    samtools view -@ ${task.cpus ?: 8} -b -q \$MAPPING_MIN_QUALITY -F \$F_FLAGS -f \$f_FLAGS ${bam} > temp_filtered.bam
    
    if [ \$CUSTOM_F_FLAG -ne 0 ] || [ \$CUSTOM_f_FLAG -ne 0 ]; then
        CUSTOM_CMD="samtools view -@ ${task.cpus ?: 8} -b"
        [ \$CUSTOM_F_FLAG -ne 0 ] && CUSTOM_CMD="\$CUSTOM_CMD -F \$CUSTOM_F_FLAG"
        [ \$CUSTOM_f_FLAG -ne 0 ] && CUSTOM_CMD="\$CUSTOM_CMD -f \$CUSTOM_f_FLAG"
        CUSTOM_CMD="\$CUSTOM_CMD temp_filtered.bam"
        \$CUSTOM_CMD > temp_filtered2.bam
        mv temp_filtered2.bam temp_filtered.bam
    fi
    
    if [ "\$MAX_SOFT_CLIP" != "no" ] && [ "\$MAX_SOFT_CLIP" -gt 0 ] || [ \$MIN_ALIGN_SCORE -gt 0 ]; then
        samtools view -h temp_filtered.bam | \\
        awk -v max_clip="\$MAX_SOFT_CLIP" -v min_score="\$MIN_ALIGN_SCORE" '
        BEGIN {OFS="\\t"}
        /^@/ {print; next}
        {
            as_score = 0
            for(i=12; i<=NF; i++) {
                if(\$i ~ /^AS:i:/) {
                    split(\$i, a, ":")
                    as_score = a[3]
                    break
                }
            }
            if(min_score > 0 && as_score < min_score) next
            
            if(max_clip != "no" && max_clip > 0) {
                cigar = \$6
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
    
    # Remove orphans - output to unique temp name
    if [ "\$REMOVE_ORPHANS" = "yes" ]; then
        samtools sort -@ ${task.cpus ?: 8} -n temp_filtered.bam -o temp_namesorted.bam
        samtools fixmate -@ ${task.cpus ?: 8} -r temp_namesorted.bam temp_fixmate.bam
        samtools sort -@ ${task.cpus ?: 8} temp_fixmate.bam -o temp_orphan_removed.bam
        rm temp_namesorted.bam temp_fixmate.bam temp_filtered.bam
    else
        # Keep temp_filtered.bam as is
        mv temp_filtered.bam temp_orphan_removed.bam
    fi
    
    # Final fixmate/markdup - always use temp_orphan_removed.bam
    samtools sort -@ ${task.cpus ?: 8} -n -o temp.nsrt.bam temp_orphan_removed.bam
    samtools fixmate -@ ${task.cpus ?: 8} -m temp.nsrt.bam temp.fxmt.bam
    samtools sort -@ ${task.cpus ?: 8} -o temp.csrt.bam temp.fxmt.bam
    samtools markdup -@ ${task.cpus ?: 8} temp.csrt.bam ${bam.simpleName}.filtered.bam
    
    rm -f temp_orphan_removed.bam temp.nsrt.bam temp.fxmt.bam temp.csrt.bam
    
    samtools index -@ ${task.cpus ?: 8} ${bam.simpleName}.filtered.bam
    """
}