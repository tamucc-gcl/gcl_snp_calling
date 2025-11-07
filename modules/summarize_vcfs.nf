process SUMMARIZE_VCFS {
    tag "${vcf.baseName}"
    publishDir params.output_dir, mode: 'copy'

    input:
    path vcf
    path ploidy_map 

    output:
    path("${vcf.simpleName}.stats.txt")
    path("${vcf.simpleName}_summary_plots.png")
    path("${vcf.simpleName}_pca.png")

    script:
    def has_ploidy_map = ploidy_map.name != 'NO_FILE'
    def ploidy_arg = has_ploidy_map ? "${ploidy_map}" : "NO_PLOIDY_MAP"
    """
    # Generate bcftools stats
    bcftools stats ${vcf} > ${vcf.simpleName}.stats.txt

    # Run R script for QC plots
    # Copy the R script from the module directory
    cp ${projectDir}/modules/vcf_qc_plots.R .
    
    # Check if VCF has enough variants for analysis
    VARIANT_COUNT=\$(zcat ${vcf} | grep -v '^#' | wc -l)
    echo "Total variants in VCF: \$VARIANT_COUNT"
    
    # Debug: Show what we're passing to R
    echo "Has ploidy map: ${has_ploidy_map}"
    echo "Ploidy argument to R: ${ploidy_arg}"

    if [ \$VARIANT_COUNT -gt 10 ]; then
        echo "Running QC plot generation..."
        # Pass the ploidy argument directly without bash variable expansion issues
        Rscript vcf_qc_plots.R ${vcf} ${vcf.simpleName} ${ploidy_arg}
    else
        echo "Too few variants (\$VARIANT_COUNT) for meaningful QC plots. Creating placeholder images..."
        
        # Create placeholder images using R
        Rscript -e "
        library(ggplot2)
        placeholder <- ggplot() + 
            annotate('text', x = 0.5, y = 0.5, 
                    label = 'Too few variants (n=\$VARIANT_COUNT)\\nfor QC analysis', 
                    size = 8) +
            theme_void() +
            theme(panel.border = element_rect(fill=NA))
        ggsave('${vcf.simpleName}_summary_plots.png', placeholder, width=10, height=10)
        ggsave('${vcf.simpleName}_pca.png', placeholder, width=5, height=5)
        "
    fi
    
    echo "Summarization complete for ${vcf}"
    """
}