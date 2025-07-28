process SUMMARIZE_VCFS {
    tag "${vcf.baseName}"
    publishDir params.output_dir, mode: 'copy'

    input:
    path vcf

    output:
    path("${vcf.simpleName}.stats.txt")

    script:
    """
    bcftools stats ${vcf} > ${vcf.simpleName}.stats.txt
    """
}