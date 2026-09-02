// modules/decompose_vcf.nf - Decompose the combined FreeBayes callset
//
// Runs once on the combined VCF, after COMBINE_VCFS. Measured on this project:
// 5,746,640 records x 874 samples decomposed in 8 minutes, so there is no case
// for fanning this out per chunk. Doing it after the concat also removes an
// entire class of bug: norm left-aligns, and 880,464 records (15% of the
// callset) were realigned in that run. Any of those sitting near a chunk
// boundary could shift upstream of its own chunk's start coordinate, which
// per-chunk decomposition would then have to repair with `concat -a` plus a
// post-concat deduplication pass. Post-concat decomposition never creates the
// problem.
//
// Pipeline:
//   norm -m -any   split multiallelic records into one record per ALT
//   norm -a        atomize MNPs and complex alleles into individual SNPs
//   norm -d exact  atomization can emit duplicates; deduplicate
//   +fill-tags     recompute AC/AN/AF/NS/MAF/F_MISSING, which are all invalid
//                  after splitting and atomization
//
// Deliberately does NOT filter by variant type. Indels are retained because
// they are what lets downstream filtering remove SNPs near indels, which
// matters given the chimeric-junction rate in this dataset.
//
// --atom-overlaps '.'  writes a missing genotype where atomized alleles
//                      overlap, rather than inserting a star allele. Avoids
//                      fake multiallelic records at the cost of a small
//                      call-rate reduction. Set params.atom_overlaps = '*'
//                      to keep them instead.
//
// --check-ref w  warns on REF mismatch rather than exiting, and the warning
//                count is reported. A hard exit here would discard an
//                otherwise complete callset over a handful of records.

process DECOMPOSE_VCF {
    tag "decompose_${vcf.simpleName}"
    publishDir params.output_dir, mode: 'copy', pattern: "*.vcf.gz*"
    publishDir "${params.output_dir}/pipeline", mode: 'copy', pattern: "*.log"
    publishDir "${params.output_dir}/snp_qc", mode: 'copy', pattern: "decompose_summary.tsv"

    input:
    path vcf
    path vcf_index
    path reference
    path reference_fai
    val output_name

    output:
    path output_name,              emit: vcf
    path "${output_name}.tbi",     emit: index
    path "decompose_summary.tsv",  emit: summary
    path "decompose.log",           emit: log

    script:
    def atom_overlaps = params.atom_overlaps != null ? params.atom_overlaps : '.'
    """
    set -euo pipefail

    LOG=decompose.log

    {
        echo "=== DECOMPOSE_VCF ==="
        echo "Input:     ${vcf}"
        echo "Output:    ${output_name}"
        echo "Reference: ${reference}"
        echo "--atom-overlaps '${atom_overlaps}'"
        bcftools --version | head -1
        echo
    } > \${LOG}

    N_IN=\$(bcftools index -n "${vcf}" 2>/dev/null || bcftools view -H "${vcf}" | wc -l)
    N_SAMP_IN=\$(bcftools query -l "${vcf}" | wc -l)
    echo "Records in: \${N_IN}"  >> \${LOG}
    echo "Samples in: \${N_SAMP_IN}" >> \${LOG}
    echo >> \${LOG}

    if [ "\${N_SAMP_IN}" -eq 0 ]; then
        echo "ERROR: no samples in ${vcf}" >&2
        exit 1
    fi

    # --- decompose -----------------------------------------------------------
    # Each norm stage gets its own stderr file so its counters can be
    # attributed. Piped together they would interleave into one stream.
    #
    # Note that atomization is NOT counted in norm's 'split' column. The
    # evidence it did anything is the record count rising between the
    # stage-2 and stage-3 totals.
    bcftools norm -m -any -f "${reference}" -c w \
             --threads ${task.cpus} -Ou "${vcf}" 2> norm1.err \
      | bcftools norm -a --atom-overlaps '${atom_overlaps}' -f "${reference}" -c w \
               --threads ${task.cpus} -Ou 2> norm2.err \
      | bcftools norm -d exact --threads ${task.cpus} -Ou 2> norm3.err \
      | bcftools +fill-tags --threads ${task.cpus} -Oz -o "${output_name}" \
               -- -t AC,AN,AF,NS,MAF,F_MISSING 2> filltags.err

    tabix -f -p vcf "${output_name}"

    {
        echo "--- bcftools norm counters ---"
        echo "[stage 1: split multiallelics]"
        cat norm1.err
        echo "[stage 2: atomize]"
        cat norm2.err
        echo "[stage 3: deduplicate]"
        cat norm3.err
        echo "[fill-tags]"
        cat filltags.err
        echo
    } >> \${LOG}

    # --- parse counters into a machine-readable summary ----------------------
    # The 'Lines' line format is:
    #   total/split/joined/realigned/mismatch_removed/dup_removed/skipped
    parse_norm() {
        local f="\$1" field="\$2"
        grep -h '^Lines' "\$f" 2>/dev/null \
          | tail -n1 \
          | sed 's/.*:[[:space:]]*//' \
          | cut -d'/' -f"\$field" \
          || echo 0
    }

    N_SPLIT=\$(parse_norm norm1.err 2)
    N_REALIGNED=\$(parse_norm norm1.err 4)
    N_STAGE2_IN=\$(parse_norm norm2.err 1)
    N_STAGE3_IN=\$(parse_norm norm3.err 1)
    N_DUP_REMOVED=\$(parse_norm norm3.err 6)

    # REF mismatches are warned about rather than removed under -c w, so they
    # appear as free-text stderr lines instead of in the counter row.
    N_REF_MISMATCH=\$(grep -hic 'mismatch' norm1.err norm2.err 2>/dev/null | paste -sd+ | bc || echo 0)

    N_OUT=\$(bcftools index -n "${output_name}")
    N_SAMP_OUT=\$(bcftools query -l "${output_name}" | wc -l)

    {
        printf "raw_records\\t%s\\n"           "\${N_IN}"
        printf "decomposed_records\\t%s\\n"   "\${N_OUT}"
        printf "n_samples\\t%s\\n"            "\${N_SAMP_OUT}"
        printf "multiallelic_split\\t%s\\n"   "\${N_SPLIT}"
        printf "realigned\\t%s\\n"            "\${N_REALIGNED}"
        printf "atomize_input\\t%s\\n"        "\${N_STAGE2_IN}"
        printf "dedup_input\\t%s\\n"          "\${N_STAGE3_IN}"
        printf "exact_duplicates_removed\\t%s\\n" "\${N_DUP_REMOVED}"
        printf "ref_mismatch_warnings\\t%s\\n"    "\${N_REF_MISMATCH}"
    } > decompose_summary.tsv

    # --- guards --------------------------------------------------------------
    # Sample count must survive decomposition. A change means the output would
    # silently mislabel every genotype column.
    if [ "\${N_SAMP_IN}" -ne "\${N_SAMP_OUT}" ]; then
        echo "ERROR: sample count changed \${N_SAMP_IN} -> \${N_SAMP_OUT}" >&2
        exit 80
    fi

    if [ "\${N_OUT}" -eq 0 ]; then
        echo "ERROR: decomposition produced an empty callset from \${N_IN} records." >&2
        exit 81
    fi

    {
        echo "--- SUMMARY ---"
        cat decompose_summary.tsv
        echo
        echo "Sample count check OK: \${N_SAMP_OUT}"
        echo "File size: \$(ls -lh ${output_name} | awk '{print \$5}')"
    } >> \${LOG}

    if [ "\${N_REF_MISMATCH}" -gt 0 ]; then
        echo "NOTE: \${N_REF_MISMATCH} REF-mismatch warnings. Records were kept" >> \${LOG}
        echo "      (--check-ref w). A large count suggests the reference passed" >> \${LOG}
        echo "      here is not the one the BAMs were mapped against." >> \${LOG}
    fi

    cat \${LOG}
    """
}
