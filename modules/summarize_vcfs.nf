process SUMMARIZE_VCFS {
    tag "${vcf.baseName}"

    input:
    path vcf

    output:
    path("${vcf.simpleName}.stats.txt")

    script:
    """
    bcftools stats ${vcf} > ${vcf.simpleName}.stats.txt
    """
}