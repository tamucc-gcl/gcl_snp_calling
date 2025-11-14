// modules/pcangsd.nf - PCAngsd analysis on ANGSD genotype likelihoods
// 
// This module runs PCAngsd on the combined Beagle file output from ANGSD
// to perform population genetic analyses including:
// - Principal Component Analysis (PCA) from genotype likelihoods
// - Admixture analysis with automatic K selection
// - Neighbor-joining tree construction
// - MAF estimation and nucleotide diversity calculation

process PCANGSD {
    tag "pcangsd_${output_prefix}"
    publishDir params.output_dir, mode: 'copy'
    
    input:
    path beagle_file
    val output_prefix
    
    output:
    path "${output_prefix}.pcangsd.cov", emit: covariance
    path "${output_prefix}.pcangsd.maf.npy", emit: maf, optional: true
    path "${output_prefix}.pcangsd.pi.npy", emit: pi, optional: true
    path "${output_prefix}.pcangsd.admix.Q", emit: admix_q, optional: true
    path "${output_prefix}.pcangsd.admix.P", emit: admix_p, optional: true
    path "${output_prefix}.pcangsd.tree.nwk", emit: tree, optional: true
    path "${output_prefix}.pcangsd.log", emit: log
    path "${output_prefix}_pcangsd_summary.txt", emit: summary
    
    script:
    """
    echo "Running PCAngsd on: ${beagle_file}"
    echo "Output prefix: ${output_prefix}"
    echo "Available threads: ${task.cpus}"
    
    # Verify beagle file exists and has data
    if [ ! -f "${beagle_file}" ]; then
        echo "ERROR: Beagle file not found!"
        exit 1
    fi
    
    SITES=\$(zcat ${beagle_file} | tail -n +2 | wc -l)
    echo "Sites in Beagle file: \$SITES"
    
    if [ \$SITES -eq 0 ]; then
        echo "ERROR: Beagle file is empty!"
        exit 1
    fi
    
    # Run PCAngsd with all requested analyses
    echo "Starting PCAngsd analysis..."
    pcangsd \\
        -b ${beagle_file} \\
        -t ${task.cpus} \\
        -o ${output_prefix}.pcangsd \\
        --tree \\
        --maf_save \\
        --pi_save \\
        --admix \\
        --admix_auto 10 \\
        2>&1 | tee ${output_prefix}.pcangsd.log
    
    PCANGSD_EXIT=\$?
    
    echo "PCAngsd exit status: \$PCANGSD_EXIT"
    
    if [ \$PCANGSD_EXIT -ne 0 ]; then
        echo "ERROR: PCAngsd failed!"
        exit 1
    fi
    
    # Create summary file
    cat > ${output_prefix}_pcangsd_summary.txt << EOF
PCAngsd Analysis Summary
========================
Timestamp: \$(date)
Input file: ${beagle_file}
Output prefix: ${output_prefix}
Threads used: ${task.cpus}

Input Data:
-----------
Sites processed: \$SITES

Output Files:
-------------
EOF
    
    # Check which outputs were created
    if [ -f "${output_prefix}.pcangsd.cov" ]; then
        N_SAMPLES=\$(head -n1 ${output_prefix}.pcangsd.cov | wc -w)
        echo "- Covariance matrix: ${output_prefix}.pcangsd.cov (\$N_SAMPLES samples)" >> ${output_prefix}_pcangsd_summary.txt
    fi
    
    if [ -f "${output_prefix}.pcangsd.maf.npy" ]; then
        echo "- MAF estimates: ${output_prefix}.pcangsd.maf.npy" >> ${output_prefix}_pcangsd_summary.txt
    fi
    
    if [ -f "${output_prefix}.pcangsd.pi.npy" ]; then
        echo "- Nucleotide diversity: ${output_prefix}.pcangsd.pi.npy" >> ${output_prefix}_pcangsd_summary.txt
    fi
    
    if [ -f "${output_prefix}.pcangsd.admix.Q" ]; then
        K=\$(head -n1 ${output_prefix}.pcangsd.admix.Q | wc -w)
        echo "- Admixture proportions: ${output_prefix}.pcangsd.admix.Q (K=\$K)" >> ${output_prefix}_pcangsd_summary.txt
    fi
    
    if [ -f "${output_prefix}.pcangsd.admix.P" ]; then
        echo "- Admixture frequencies: ${output_prefix}.pcangsd.admix.P" >> ${output_prefix}_pcangsd_summary.txt
    fi
    
    if [ -f "${output_prefix}.pcangsd.tree.nwk" ]; then
        echo "- Neighbor-joining tree: ${output_prefix}.pcangsd.tree.nwk" >> ${output_prefix}_pcangsd_summary.txt
    fi
    
    echo "" >> ${output_prefix}_pcangsd_summary.txt
    echo "Analysis completed successfully!" >> ${output_prefix}_pcangsd_summary.txt
    
    # Display summary
    echo ""
    echo "========================================"
    cat ${output_prefix}_pcangsd_summary.txt
    echo "========================================"
    
    echo "PCAngsd analysis complete!"
    """
}