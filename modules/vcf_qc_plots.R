#!/usr/bin/env Rscript

#### Parse command line arguments ####
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript vcf_qc_plots.R <input.vcf.gz> <output_prefix>")
}

vcf_file <- args[1]
output_prefix <- args[2]

cat("Processing VCF:", vcf_file, "\n")
cat("Output prefix:", output_prefix, "\n")

#### Libraries ####
library(adegenet)
library(vcfR)
library(patchwork)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)

#### Data ####
raw_vcf <- read.vcfR(vcf_file)

#### QC Plots ####
depth_data <- extract.gt(raw_vcf, 
                         element = "DP", 
                         as.numeric = TRUE)

locus_missingness_plot <- depth_data %>%
  is.na() %>%
  rowMeans(na.rm = TRUE) %>%
  enframe(name = 'locus',
          value = 'pct_samples_missing') %>%
  ggplot(aes(x = pct_samples_missing)) +
  geom_histogram(bins = 50) +
  scale_x_continuous(labels = scales::percent_format()) +
  scale_y_continuous(labels = scales::comma_format()) +
  coord_cartesian(xlim = c(0, 1), clip = 'on') +
  labs(x = '% of samples missing locus',
       y = 'Number of Loci') 

sample_missingness_plot <- depth_data %>%
  is.na() %>%
  colMeans(na.rm = TRUE) %>%
  enframe(name = 'sample',
          value = 'pct_samples_missing') %>%
  ggplot(aes(x = pct_samples_missing)) +
  geom_histogram(bins = 50) +
  scale_x_continuous(labels = scales::percent_format()) +
  scale_y_continuous(labels = scales::comma_format()) +
  coord_cartesian(xlim = c(0, 1), clip = 'on') +
  labs(x = '% of loci missing in sample',
       y = 'Number of Samples')


depth_plot <- depth_data %>%
  as.numeric() %>%
  tibble(depth = .) %>%
  mutate(depth = replace_na(depth, 0L)) %>%
  filter(depth > 0) %>%
  ggplot(aes(x = depth)) +
  geom_histogram(bins = 50) +
  scale_x_continuous(transform = scales::log10_trans(),
                     labels = scales::comma_format(),
                     guide = guide_axis_logticks()) +
  scale_y_continuous(transform = scales::log10_trans(),
                     labels = scales::comma_format(), 
                     guide = guide_axis_logticks()) +
  labs(x = 'Sequencing Depth',
       y = 'Number of Sites')


meanDepth_plot <- depth_data %>%
  rowMeans(na.rm = TRUE) %>%
  enframe(name = 'locus',
          value = 'mean_depth') %>%
  ggplot(aes(x = mean_depth)) +
  geom_histogram(bins = 50) +
  scale_x_continuous(transform = scales::log10_trans(),
                     labels = scales::comma_format(),
                     guide = guide_axis_logticks()) +
  scale_y_continuous(transform = scales::log10_trans(),
                     labels = scales::comma_format(), 
                     guide = guide_axis_logticks()) +
  labs(x = 'Mean Sequencing Depth',
       y = 'Number of Loci')

snp_summary_stats_plot  <- locus_missingness_plot + 
  sample_missingness_plot +
  depth_plot +
  meanDepth_plot + 
  plot_annotation(title = 'Raw SNP Summary Plots') &
  theme_classic(base_size = 16) +
  theme(panel.background = element_rect(colour = 'black')) 
ggsave(paste0(output_prefix, '_summary_plots.png'),
       plot = snp_summary_stats,
       height = 10,
       width = 10)

#### PCA ####
raw_genlight <- vcfR2genlight(raw_vcf)
pca_out <- glPca(raw_genlight,
                 nf = 2, 
                 loadings = FALSE)

pct_var <- scales::percent(pca_out$eig / sum(pca_out$eig), accuracy = 0.1)

#Identify outliers (Mahalanobis or Robust Mahalanobis distance (MCD estimator))
flag_outliers <- function(data, pca, alpha = 0.975){
  
  md <- mahalanobis(pca$scores, colMeans(pca$scores), cov(pca$scores))
  
  mcd <- MASS::cov.mcd(pca$scores)
  md_robust <- mahalanobis(pca$scores, mcd$center, mcd$cov)
  
  
  mutate(data,
         outlier_score = md,
         outlier = md > qchisq(alpha, 
                               df = ncol(pca$scores)),
         robust_outlier_score = md_robust,
         robust_outlier = md_robust > qchisq(alpha, 
                                      df = ncol(pca$scores)))
}

raw_pca_scores <- pca_out$scores %>%
  as_tibble(rownames = 'sample_id') %>%
  flag_outliers(pca_out)

raw_snp_pca <- raw_pca_scores %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_point() +
  geom_text_repel(data = . %>% filter(outlier),
            aes(label = sample_id),
            hjust = 'inward',
            vjust = 'inward') +
  labs(x = str_c('PC1 (', pct_var[1], ")"),
       y = str_c('PC1 (', pct_var[2], ")")) +
  theme_classic(base_size = 16) +
  theme(panel.background = element_rect(colour = 'black')) 

ggsave(paste0(output_prefix, '_pca.png'),
       plot = snp_pca,
       height = 5,
       width = 5)

cat("QC plots successfully generated!\n")