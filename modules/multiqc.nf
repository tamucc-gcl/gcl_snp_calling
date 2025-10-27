// modules/multiqc.nf
process multiqc {
    label 'multiqc'
    tag "${step_name}"
    
    // Only publish the HTML report to publishDir
    publishDir "${params.outdir}/multiqc_reports", mode: 'copy', pattern: "*.html"

    input:
        path(input_files)
        val(step_name)

    output:
        path("multiqc_${step_name}.html")
        path("multiqc_${step_name}_general_stats.txt")

    script:
    """
    # Run MultiQC
    multiqc \\
        --title "MultiQC Report - ${step_name}" \\
        --filename multiqc_${step_name}.html \\
        --export \\
        --data-dir \\
        .
    
    # The data directory is named after the output filename
    # So multiqc_${step_name}.html creates multiqc_${step_name}_data directory
    DATA_DIR="multiqc_${step_name}_data"
    
    echo "Looking for data directory: \$DATA_DIR"
    ls -la \$DATA_DIR/ 2>/dev/null || echo "Data directory not found"
    
    # Copy the general stats file from the correct location
    if [ -f "\$DATA_DIR/multiqc_general_stats.txt" ]; then
        echo "Found general stats file in \$DATA_DIR"
        cp "\$DATA_DIR/multiqc_general_stats.txt" "multiqc_${step_name}_general_stats.txt"
        echo "Copied general stats file with size: \$(wc -l multiqc_${step_name}_general_stats.txt)"
    else
        echo "General stats file not found in \$DATA_DIR"
        echo "Creating empty general stats file"
        touch "multiqc_${step_name}_general_stats.txt"
    fi
    
    # Verify the final output
    echo "Final output files:"
    ls -la multiqc_${step_name}*
    """
}