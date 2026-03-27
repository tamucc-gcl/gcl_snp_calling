process REPORT_SNP_CALLING_SUMMARY {
    tag "${prefix}"
    publishDir params.output_dir, mode: 'copy'

    input:
    val caller
    val raw_vcf_link
    tuple val(prefix),
          path(stats_txt),
          path(site_qc_tsv),
          path(sample_qc_tsv),
          path(missing_site),
          path(missing_indv),
          path(freq_tsv),
          path(standardized_summary),
          path(summary_plot),
          path(extra_plot),
          path(pca_plot),
          path(sample_qc_derived),
          path(locus_qc_derived),
          path(worst_samples),
          path(worst_loci)

    output:
    path("${prefix}_snp_calling_report.md")

    script:
    """
    set -euo pipefail

    python3 "${projectDir}/py_scripts/report_snp_calling.py" \
        "${prefix}" \
        "${caller}" \
        "${raw_vcf_link}" \
        "${stats_txt}" \
        "${standardized_summary}" \
        "${summary_plot}" \
        "${extra_plot}" \
        "${pca_plot}" \
        "${worst_samples}" \
        "${worst_loci}" \
        "${sample_qc_derived}" \
        "${locus_qc_derived}" \
        "${missing_indv}" \
        "${missing_site}" \
        "${freq_tsv}" \
        "${site_qc_tsv}" \
        "${sample_qc_tsv}"
    """
}
