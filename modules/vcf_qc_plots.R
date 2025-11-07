#!/usr/bin/env Rscript

#### Parse command line arguments ####
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript vcf_qc_plots.R <input.vcf.gz> <output_prefix>")
}

vcf_file <- args[1]
output_prefix <- args[2]
ploidy_map_file <- if (length(args) >= 3 && args[3] != "NO_PLOIDY_MAP") args[3] else NULL


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
library(forcats)

#### Data ####
raw_vcf <- read.vcfR(vcf_file)

if (!is.null(ploidy_map_file)) {
  cat("Ploidy map file:", ploidy_map_file, "\n")
  
  # Try to read the ploidy map with error handling
  tryCatch({
    ploidy_map <- read.table(ploidy_map_file, 
                             header = FALSE, 
                             col.names = c("sample", "ploidy"),
                             stringsAsFactors = FALSE,
                             comment.char = "#")
    
    cat("Loaded ploidy information for", nrow(ploidy_map), "samples\n")
    cat("Ploidy map class:", class(ploidy_map), "\n")
    cat("Ploidy map dimensions:", dim(ploidy_map), "\n")
  }, error = function(e) {
    cat("Error reading ploidy map file:", e$message, "\n")
    ploidy_map <<- NULL
  })
  
} else {
  cat("No ploidy map provided - setting to NULL\n")
  ploidy_map <- NULL
  cat("Ploidy map is NULL:", is.null(ploidy_map), "\n")
}

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
       y = 'Number of Genotypes')


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

distribution_shared_loci <- depth_data %>%
  is.na %>%
  `!` %>%
  rowSums() %>%
  tibble(number_individuals = .) %>%
  count(number_individuals) %>%
  ggplot(aes(x = number_individuals, y = n)) +
  geom_col() +
  scale_x_continuous(labels = scales::comma_format(),
                     limits = c(0, ncol(depth_data))) +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(x = 'Number of Samples with locus',
       y = 'Number of Loci')

sample_n_loci_plot <- depth_data %>%
  is.na %>%
  `!` %>%
  colSums() %>%
  enframe(name = 'sample_id',
          value = 'number_loci') %>%
  mutate(sample_id = fct_reorder(sample_id, number_loci)) %>%
  ggplot(aes(x = number_loci, y = sample_id)) +
  geom_col() +
  scale_x_continuous(labels = scales::comma_format()) +
  labs(x = 'Number of Loci',
       y = NULL)

snp_summary_stats_plot  <- locus_missingness_plot + 
  sample_missingness_plot +
  depth_plot +
  meanDepth_plot + 
  distribution_shared_loci +
  sample_n_loci_plot +
  plot_layout(ncol = 2) +
  plot_annotation(title = 'Raw SNP Summary Plots') &
  theme_classic(base_size = 16) +
  theme(panel.background = element_rect(colour = 'black')) 
ggsave(paste0(output_prefix, '_summary_plots.png'),
       plot = snp_summary_stats_plot,
       height = 15,
       width = 10)

#### PCA ####

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

genind2genlight <- function(gi){
  locna <- gi@loc.n.all
  ccc <- 1
  for (i in 2:length(locna)) {
    if (locna[i - 1] == 1) {
      ccc[i] <- ccc[i - 1] + 1
    }
    else {
      ccc[i] <- ccc[i - 1] + 2
    }
  }
  
  new("genlight", gi@tab[, ccc], pop = pop(gi), other = gi@other, 
      ploidy = ploidy(gi), loc.names = locNames(gi), ind.names = indNames(gi))
}


if(is.null(ploidy_map)){
  raw_genlight <- vcfR2genlight(raw_vcf)
  
} else {
  raw_genlight <- ploidy_map %>%
    arrange(match(sample, colnames(raw_vcf@gt)[-1])) %>%
    pull(ploidy) %>%
    vcfR2genind(raw_vcf,
                ploidy = .) %>%
    genind2genlight()
}

pca_out <- glPca(raw_genlight,
                 nf = 2, 
                 loadings = FALSE)

pct_var <- scales::percent(pca_out$eig / sum(pca_out$eig), accuracy = 0.1)


raw_pca_scores <- pca_out$scores %>%
  as_tibble(rownames = 'sample_id') %>%
  flag_outliers(pca_out)


raw_snp_pca <- raw_pca_scores %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed')

# More robust check for ploidy_map
# Check if ploidy_map exists, is a data frame, and has the expected structure
has_valid_ploidy_map <- FALSE
if(!is.null(ploidy_map)) {
  if(is.data.frame(ploidy_map)) {
    if(nrow(ploidy_map) > 0) {
      has_valid_ploidy_map <- TRUE
      cat("Valid ploidy map detected with", nrow(ploidy_map), "samples\n")
    }
  }
}

# Use the validated flag for the conditional
if(has_valid_ploidy_map && nrow(ploidy_map) < 10){
  raw_snp_pca <- raw_snp_pca +
    geom_text(aes(label = sample_id)) 
} else {
  raw_snp_pca <- raw_snp_pca +
    geom_point() +
    geom_text_repel(data = . %>% filter(outlier),
                    aes(label = sample_id),
                    hjust = 'inward',
                    vjust = 'inward')
}

raw_snp_pca <- raw_snp_pca +
  labs(x = str_c('PC1 (', pct_var[1], ")"),
       y = str_c('PC2 (', pct_var[2], ")")) +  # Fixed typo: PC2 instead of PC1
  theme_classic(base_size = 16) +
  theme(panel.background = element_rect(colour = 'black')) 

ggsave(paste0(output_prefix, '_pca.png'),
       plot = raw_snp_pca,
       height = 5,
       width = 5)

cat("QC plots successfully generated!\n")