// modules/generate_plots.nf - Generate comprehensive QC plots using R and tidyverse

process GENERATE_PLOTS {
    tag "qc_plots"
    publishDir "${params.output_dir}/qc/plots", mode: 'copy'
    conda 'r-base=4.3 r-tidyverse=2.0 r-ggplot2=3.4 r-scales r-viridis r-gridextra r-cowplot r-ggrepel r-pheatmap'
    
    input:
    path variant_stats
    path sample_stats
    path chr_stats
    path af_dist
    path depth_dist
    path qual_dist
    path sample_metrics
    path sample_heterozygosity
    path sample_depth_stats
    path sample_missingness
    path sample_relatedness
    path locus_metrics
    path pool_metrics
    
    output:
    path "*.pdf", emit: plots
    path "*.png", emit: pngs
    path "qc_plot_data.RData", emit: rdata
    
    script:
    """
    #!/usr/bin/env Rscript
    
    # Load required libraries
    suppressPackageStartupMessages({
        library(tidyverse)
        library(ggplot2)
        library(scales)
        library(viridis)
        library(gridExtra)
        library(cowplot)
        library(ggrepel)
        library(pheatmap)
    })
    
    # Set theme
    theme_set(theme_bw(base_size = 12))
    
    # Color palette
    qc_colors <- c("#2E7D32", "#1976D2", "#F57C00", "#C62828", "#7B1FA2", "#00838F")
    
    # Read all data files
    cat("Reading QC data files...\\n")
    
    # Sample-level data
    sample_stats_df <- read_tsv("${sample_stats}", show_col_types = FALSE)
    sample_het <- read_tsv("${sample_heterozygosity}", show_col_types = FALSE)
    sample_depth <- read_tsv("${sample_depth_stats}", show_col_types = FALSE)
    sample_miss <- read_tsv("${sample_missingness}", show_col_types = FALSE)
    sample_rel <- read_tsv("${sample_relatedness}", show_col_types = FALSE)
    
    # Variant-level data
    chr_stats_df <- read_tsv("${chr_stats}", show_col_types = FALSE)
    af_dist_df <- read_tsv("${af_dist}", show_col_types = FALSE)
    depth_dist_df <- read_tsv("${depth_dist}", show_col_types = FALSE)
    qual_dist_df <- read_tsv("${qual_dist}", show_col_types = FALSE)
    
    # Locus-level data (sample first 10000 for plotting)
    locus_qc <- read_tsv("${locus_metrics}", show_col_types = FALSE, n_max = 10000)
    
    # Pool metrics if available
    if (file.size("${pool_metrics}") > 50) {
        pool_df <- read_tsv("${pool_metrics}", show_col_types = FALSE)
        has_pools <- TRUE
    } else {
        has_pools <- FALSE
    }
    
    # Create plots directory structure
    dir.create("sample_qc", showWarnings = FALSE)
    dir.create("variant_qc", showWarnings = FALSE)
    dir.create("locus_qc", showWarnings = FALSE)
    
    # ============= SAMPLE-LEVEL QC PLOTS =============
    
    cat("Generating sample-level QC plots...\\n")
    
    # 1. Sample missingness barplot
    p_miss <- ggplot(sample_miss, aes(x = reorder(sample, missingness_rate), 
                                       y = missingness_rate)) +
        geom_bar(stat = "identity", fill = qc_colors[1]) +
        geom_hline(yintercept = 0.1, linetype = "dashed", color = "red") +
        coord_flip() +
        scale_y_continuous(labels = percent_format()) +
        labs(title = "Sample Missingness Rate",
             subtitle = "Red line indicates 10% threshold",
             x = "Sample",
             y = "Missingness Rate") +
        theme(axis.text.y = element_text(size = 8))
    
    ggsave("sample_qc/sample_missingness.pdf", p_miss, width = 10, height = 8)
    ggsave("sample_qc/sample_missingness.png", p_miss, width = 10, height = 8, dpi = 300)
    
    # 2. Sample heterozygosity vs missingness
    sample_combined <- sample_het %>%
        left_join(sample_miss, by = "sample")
    
    p_het_miss <- ggplot(sample_combined, aes(x = missingness_rate, 
                                               y = heterozygosity_rate)) +
        geom_point(size = 3, alpha = 0.7, color = qc_colors[2]) +
        geom_text_repel(aes(label = sample), size = 3, max.overlaps = 15) +
        geom_vline(xintercept = 0.1, linetype = "dashed", color = "red", alpha = 0.5) +
        geom_hline(yintercept = mean(sample_combined\$heterozygosity_rate, na.rm = TRUE) + 
                              3 * sd(sample_combined\$heterozygosity_rate, na.rm = TRUE),
                   linetype = "dashed", color = "orange", alpha = 0.5) +
        geom_hline(yintercept = mean(sample_combined\$heterozygosity_rate, na.rm = TRUE) - 
                              3 * sd(sample_combined\$heterozygosity_rate, na.rm = TRUE),
                   linetype = "dashed", color = "orange", alpha = 0.5) +
        scale_x_continuous(labels = percent_format()) +
        labs(title = "Sample Heterozygosity vs Missingness",
             subtitle = "Outliers may indicate sample quality issues or contamination",
             x = "Missingness Rate",
             y = "Heterozygosity Rate")
    
    ggsave("sample_qc/heterozygosity_vs_missingness.pdf", p_het_miss, width = 10, height = 8)
    ggsave("sample_qc/heterozygosity_vs_missingness.png", p_het_miss, width = 10, height = 8, dpi = 300)
    
    # 3. Sample depth distribution
    p_depth_box <- ggplot(sample_depth, aes(x = reorder(sample, mean_depth), 
                                            y = mean_depth)) +
        geom_boxplot(aes(ymin = mean_depth - sd_depth, 
                        ymax = mean_depth + sd_depth,
                        lower = mean_depth - sd_depth/2,
                        upper = mean_depth + sd_depth/2,
                        middle = median_depth),
                    stat = "identity", fill = qc_colors[3], alpha = 0.7) +
        geom_point(aes(y = median_depth), color = "red", size = 2) +
        coord_flip() +
        labs(title = "Sample Read Depth Distribution",
             subtitle = "Box: mean ± SD, Red dot: median",
             x = "Sample",
             y = "Read Depth") +
        theme(axis.text.y = element_text(size = 8))
    
    ggsave("sample_qc/sample_depth_distribution.pdf", p_depth_box, width = 10, height = 10)
    ggsave("sample_qc/sample_depth_distribution.png", p_depth_box, width = 10, height = 10, dpi = 300)
    
    # 4. Sample relatedness heatmap (if multiple samples)
    if (nrow(sample_rel) > 0) {
        # Create matrix from pairwise data
        samples_unique <- unique(c(sample_rel\$sample1, sample_rel\$sample2))
        rel_matrix <- matrix(0, nrow = length(samples_unique), ncol = length(samples_unique))
        rownames(rel_matrix) <- samples_unique
        colnames(rel_matrix) <- samples_unique
        
        for (i in 1:nrow(sample_rel)) {
            s1 <- sample_rel\$sample1[i]
            s2 <- sample_rel\$sample2[i]
            k <- sample_rel\$kinship_coefficient[i]
            rel_matrix[s1, s2] <- k
            rel_matrix[s2, s1] <- k
        }
        diag(rel_matrix) <- 0.5  # Self-relatedness
        
        # Plot heatmap
        pdf("sample_qc/sample_relatedness_heatmap.pdf", width = 12, height = 10)
        pheatmap(rel_matrix,
                color = colorRampPalette(c("white", "yellow", "red"))(50),
                main = "Sample Relatedness (Kinship Coefficient)",
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                show_rownames = TRUE,
                show_colnames = TRUE,
                annotation_legend = TRUE)
        dev.off()
        
        # Also save as PNG
        png("sample_qc/sample_relatedness_heatmap.png", width = 1200, height = 1000, res = 100)
        pheatmap(rel_matrix,
                color = colorRampPalette(c("white", "yellow", "red"))(50),
                main = "Sample Relatedness (Kinship Coefficient)",
                cluster_rows = TRUE,
                cluster_cols = TRUE)
        dev.off()
    }
    
    # 5. Inbreeding coefficient distribution
    p_inbreed <- ggplot(sample_het, aes(x = inbreeding_coefficient)) +
        geom_histogram(bins = 30, fill = qc_colors[4], alpha = 0.7, color = "black") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
        labs(title = "Distribution of Inbreeding Coefficients",
             subtitle = "F > 0 indicates excess homozygosity",
             x = "Inbreeding Coefficient (F)",
             y = "Number of Samples")
    
    ggsave("sample_qc/inbreeding_distribution.pdf", p_inbreed, width = 8, height = 6)
    ggsave("sample_qc/inbreeding_distribution.png", p_inbreed, width = 8, height = 6, dpi = 300)
    
    # ============= VARIANT-LEVEL QC PLOTS =============
    
    cat("Generating variant-level QC plots...\\n")
    
    # 6. Chromosome variant distribution
    p_chr <- ggplot(chr_stats_df, aes(x = reorder(chromosome, -total_variants), 
                                      y = total_variants)) +
        geom_bar(stat = "identity", aes(fill = "Total"), alpha = 0.7) +
        geom_bar(stat = "identity", aes(y = snps, fill = "SNPs"), alpha = 0.7) +
        geom_bar(stat = "identity", aes(y = indels, fill = "INDELs"), alpha = 0.7) +
        scale_fill_manual(values = qc_colors[1:3]) +
        labs(title = "Variant Distribution by Chromosome",
             x = "Chromosome",
             y = "Number of Variants",
             fill = "Type") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave("variant_qc/chromosome_distribution.pdf", p_chr, width = 12, height = 6)
    ggsave("variant_qc/chromosome_distribution.png", p_chr, width = 12, height = 6, dpi = 300)
    
    # 7. Allele frequency spectrum
    # Order AF bins properly
    af_order <- c("0-0.01", "0.01-0.05", "0.05-0.1", "0.1-0.25", 
                  "0.25-0.5", "0.5-0.75", "0.75-0.95", "0.95-1.0")
    af_dist_df\$af_bin <- factor(af_dist_df\$af_bin, levels = af_order)
    
    p_af <- ggplot(af_dist_df, aes(x = af_bin, y = count, fill = variant_type)) +
        geom_bar(stat = "identity", position = "dodge") +
        scale_y_log10(labels = comma_format()) +
        scale_fill_manual(values = qc_colors[1:2]) +
        labs(title = "Allele Frequency Distribution",
             subtitle = "Log scale; Expected U-shape for neutral variation",
             x = "Allele Frequency Bin",
             y = "Count (log10)",
             fill = "Variant Type") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave("variant_qc/allele_frequency_spectrum.pdf", p_af, width = 10, height = 6)
    ggsave("variant_qc/allele_frequency_spectrum.png", p_af, width = 10, height = 6, dpi = 300)
    
    # 8. Depth distribution
    # Order depth bins properly
    depth_order <- c("0-5", "5-10", "10-20", "20-30", "30-50", "50-100", "100-200", "200+")
    depth_dist_df\$depth_bin <- factor(depth_dist_df\$depth_bin, levels = depth_order)
    
    p_depth_dist <- ggplot(depth_dist_df, aes(x = depth_bin, y = count)) +
        geom_bar(stat = "identity", fill = qc_colors[3], alpha = 0.7) +
        geom_vline(xintercept = which(depth_order == "5-10"), 
                  linetype = "dashed", color = "red", alpha = 0.5) +
        scale_y_continuous(labels = comma_format()) +
        labs(title = "Variant Depth Distribution",
             subtitle = "Red line indicates low depth threshold",
             x = "Depth Bin",
             y = "Number of Variants") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave("variant_qc/depth_distribution.pdf", p_depth_dist, width = 10, height = 6)
    ggsave("variant_qc/depth_distribution.png", p_depth_dist, width = 10, height = 6, dpi = 300)
    
    # 9. Quality score distribution
    qual_order <- c("0-10", "10-20", "20-30", "30-50", "50-100", "100-200", "200-500", "500+")
    qual_dist_df\$quality_bin <- factor(qual_dist_df\$quality_bin, levels = qual_order)
    
    p_qual_dist <- ggplot(qual_dist_df, aes(x = quality_bin, y = count)) +
        geom_bar(stat = "identity", fill = qc_colors[5], alpha = 0.7) +
        geom_vline(xintercept = which(qual_order == "20-30"), 
                  linetype = "dashed", color = "red", alpha = 0.5) +
        scale_y_log10(labels = comma_format()) +
        labs(title = "Variant Quality Score Distribution",
             subtitle = "Log scale; Red line indicates typical quality threshold",
             x = "Quality Score Bin",
             y = "Number of Variants (log10)") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave("variant_qc/quality_distribution.pdf", p_qual_dist, width = 10, height = 6)
    ggsave("variant_qc/quality_distribution.png", p_qual_dist, width = 10, height = 6, dpi = 300)
    
    # ============= LOCUS-LEVEL QC PLOTS =============
    
    cat("Generating locus-level QC plots...\\n")
    
    # 10. Quality vs Depth scatter plot
    locus_subset <- locus_qc %>%
        filter(!is.na(depth) & !is.na(quality)) %>%
        sample_n(min(5000, n()))  # Subsample for plotting
    
    p_qual_depth <- ggplot(locus_subset, aes(x = depth, y = quality)) +
        geom_point(alpha = 0.3, size = 1, color = qc_colors[1]) +
        geom_smooth(method = "loess", color = "red", se = TRUE) +
        scale_x_log10(labels = comma_format()) +
        scale_y_log10(labels = comma_format()) +
        geom_vline(xintercept = 10, linetype = "dashed", alpha = 0.5) +
        geom_hline(yintercept = 30, linetype = "dashed", alpha = 0.5) +
        labs(title = "Variant Quality vs Read Depth",
             subtitle = "Dashed lines show typical filtering thresholds",
             x = "Read Depth (log10)",
             y = "Quality Score (log10)")
    
    ggsave("locus_qc/quality_vs_depth.pdf", p_qual_depth, width = 10, height = 8)
    ggsave("locus_qc/quality_vs_depth.png", p_qual_depth, width = 10, height = 8, dpi = 300)
    
    # 11. Missingness vs MAF
    p_miss_maf <- ggplot(locus_qc, aes(x = allele_freq, y = missingness_rate)) +
        geom_point(alpha = 0.3, size = 1, color = qc_colors[2]) +
        geom_smooth(method = "loess", color = "red", se = TRUE) +
        geom_hline(yintercept = 0.2, linetype = "dashed", color = "red", alpha = 0.5) +
        labs(title = "Variant Missingness vs Minor Allele Frequency",
             subtitle = "High missingness at low MAF may indicate calling issues",
             x = "Allele Frequency",
             y = "Missingness Rate")
    
    ggsave("locus_qc/missingness_vs_maf.pdf", p_miss_maf, width = 10, height = 6)
    ggsave("locus_qc/missingness_vs_maf.png", p_miss_maf, width = 10, height = 6, dpi = 300)
    
    # 12. Heterozygosity by variant type
    p_het_type <- ggplot(locus_qc, aes(x = type, y = heterozygosity, fill = type)) +
        geom_boxplot(alpha = 0.7) +
        scale_fill_manual(values = qc_colors[1:3]) +
        labs(title = "Heterozygosity Distribution by Variant Type",
             x = "Variant Type",
             y = "Heterozygosity Rate") +
        theme(legend.position = "none")
    
    ggsave("locus_qc/heterozygosity_by_type.pdf", p_het_type, width = 8, height = 6)
    ggsave("locus_qc/heterozygosity_by_type.png", p_het_type, width = 8, height = 6, dpi = 300)
    
    # ============= POOL-SPECIFIC PLOTS (if applicable) =============
    
    if (has_pools) {
        cat("Generating pool-specific QC plots...\\n")
        
        p_pool_maf <- ggplot(pool_df, aes(x = reorder(sample, -minor_allele_freq_mean), 
                                          y = minor_allele_freq_mean)) +
            geom_bar(stat = "identity", fill = qc_colors[6], alpha = 0.7) +
            geom_errorbar(aes(ymin = minor_allele_freq_median, 
                            ymax = minor_allele_freq_mean),
                         width = 0.2) +
            coord_flip() +
            labs(title = "Pool-specific Minor Allele Frequencies",
                 subtitle = "Bar: mean MAF, Error bar to median",
                 x = "Pool",
                 y = "Mean Minor Allele Frequency")
        
        ggsave("sample_qc/pool_maf_distribution.pdf", p_pool_maf, width = 10, height = 8)
        ggsave("sample_qc/pool_maf_distribution.png", p_pool_maf, width = 10, height = 8, dpi = 300)
        
        p_pool_poly <- ggplot(pool_df, aes(x = pool_size, y = polymorphic_sites)) +
            geom_point(size = 3, color = qc_colors[1]) +
            geom_smooth(method = "lm", se = TRUE, color = "red") +
            geom_text_repel(aes(label = sample), size = 3) +
            labs(title = "Polymorphic Sites vs Pool Size",
                 subtitle = "Expected: larger pools have more polymorphic sites",
                 x = "Pool Size (number of individuals)",
                 y = "Number of Polymorphic Sites")
        
        ggsave("sample_qc/pool_polymorphic_sites.pdf", p_pool_poly, width = 10, height = 8)
        ggsave("sample_qc/pool_polymorphic_sites.png", p_pool_poly, width = 10, height = 8, dpi = 300)
    }
    
    # ============= COMBINED QC DASHBOARD =============
    
    cat("Creating QC dashboard...\\n")
    
    # Create a multi-panel dashboard
    dashboard_plots <- list(
        p_miss + theme(legend.position = "none"),
        p_het_miss + theme(legend.position = "none"),
        p_depth_box + theme(legend.position = "none", axis.text.y = element_blank()),
        p_inbreed + theme(legend.position = "none"),
        p_af + theme(legend.position = "bottom"),
        p_qual_depth + theme(legend.position = "none")
    )
    
    dashboard <- plot_grid(plotlist = dashboard_plots, 
                          ncol = 2, 
                          labels = c("A", "B", "C", "D", "E", "F"))
    
    ggsave("qc_dashboard.pdf", dashboard, width = 16, height = 18)
    ggsave("qc_dashboard.png", dashboard, width = 16, height = 18, dpi = 300)
    
    # Save R data for potential reuse
    save.image("qc_plot_data.RData")
    
    cat("QC plots generation complete!\\n")
    """
}