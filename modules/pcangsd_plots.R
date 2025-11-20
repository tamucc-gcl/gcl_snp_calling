#!/usr/bin/env Rscript

# PCAngsd Results Visualization Script
# Usage: Rscript pcangsd_plots.R <output_prefix> <sample_names_file>

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggrepel)
suppressMessages(suppressWarnings(library(ggtree)))
suppressMessages(suppressWarnings(library(treeio)))
# suppressMessages(suppressWarnings(library(tidytree)))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript pcangsd_plots.R <output_prefix> <sample_names_file>")
}

output_prefix <- args[1]
sample_names_file <- args[2]

cat("Visualizing PCAngsd results for:", output_prefix, "\n")

# Read sample names
if (!file.exists(sample_names_file)) {
  stop("Sample names file not found: ", sample_names_file)
}

cat("Reading sample names from:", sample_names_file, "\n")
sample_ids <- read.table(sample_names_file, header = FALSE, stringsAsFactors = FALSE)$V1
n_samples <- length(sample_ids)
cat("Loaded", n_samples, "sample names\n")

# Read covariance matrix and perform PCA
cov_file <- paste0(output_prefix, ".pcangsd.cov")
if (!file.exists(cov_file)) {
  stop("Covariance file not found: ", cov_file)
}

cat("Reading covariance matrix...\n")
cov_matrix <- as.matrix(read.table(cov_file, header = FALSE))

# Perform eigen decomposition
cat("Performing eigenvalue decomposition...\n")
eigen_result <- eigen(cov_matrix)

# Calculate variance explained
var_explained <- eigen_result$values / sum(eigen_result$values)
pct_var <- scales::percent(var_explained, accuracy = 0.1)

#Identify outliers (Mahalanobis or Robust Mahalanobis distance (MCD estimator))
flag_outliers <- function(data, alpha = 0.975){
  pca_mat <- select(data, -sample_id) %>%
    as.matrix
  md <- mahalanobis(pca_mat, colMeans(pca_mat), cov(pca_mat))
  
  mcd <- MASS::cov.mcd(pca_mat)
  md_robust <- mahalanobis(pca_mat, mcd$center, mcd$cov)
  
  
  mutate(data,
         outlier_score = md,
         outlier = md > qchisq(alpha, 
                               df = ncol(pca_mat)),
         robust_outlier_score = md_robust,
         robust_outlier = md_robust > qchisq(alpha, 
                                             df = ncol(pca_mat)))
}

# Create PCA data frame
pca_data <- tibble(sample_id = sample_ids,
       PC1 = eigen_result$vectors[, 1],
       PC2 = eigen_result$vectors[, 2],
       PC3 = eigen_result$vectors[, 3],
       PC4 = eigen_result$vectors[, 4]) %>%
  flag_outliers()



# PCA Plot: PC1 vs PC2
p1 <- pca_data %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_point() +
  geom_text_repel(data = . %>% filter(outlier),
                  aes(label = sample_id),
                  hjust = 'inward',
                  vjust = 'inward') +
  labs(x = paste0('PC1 (', pct_var[1], ")"),
       y = paste0('PC2 (', pct_var[2], ")")) +
  theme_classic(base_size = 16) +
  theme(panel.background = element_rect(colour = 'black')) 

# PCA Plot: PC3 vs PC4
p2 <- pca_data %>%
  ggplot(aes(x = PC3, y = PC4)) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_point() +
  geom_text_repel(data = . %>% filter(outlier),
                  aes(label = sample_id),
                  hjust = 'inward',
                  vjust = 'inward') +
  labs(x = paste0('PC3 (', pct_var[3], ")"),
       y = paste0('PC4 (', pct_var[4], ")")) +
  theme_classic(base_size = 16) +
  theme(panel.background = element_rect(colour = 'black')) 

# Scree plot
scree_data <- data.frame(
  PC = 1:min(20, length(var_explained)),
  variance = var_explained[1:min(20, length(var_explained))]
)

p3 <- ggplot(scree_data, aes(x = PC, y = variance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_line(group = 1) +
  geom_point() +
  labs(
    x = "Principal Component",
    y = "Proportion of Variance Explained",
    title = "Scree Plot"
  ) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_classic(base_size = 14) +
  theme(panel.border = element_rect(fill = NA, color = "black"))

# Combine plots
combined_pca <- (p1 | p2) +#/ p3 +
  plot_annotation(
    title = "PCAngsd Population Structure Analysis",
    subtitle = paste0("Based on genotype likelihoods from ", n_samples, " samples")
  )

# Save PCA plots
cat("Saving PCA plots...\n")
ggsave(
  paste0(output_prefix, "_pcangsd_pca.png"),
  combined_pca,
  width = 14,
  height = 10,
  dpi = 300
)

# Read and plot admixture results if available
admix_q_file <- paste0(output_prefix, ".pcangsd.admix.Q")
if (file.exists(admix_q_file)) {
  cat("Reading admixture results...\n")
  admix_q <- read.table(admix_q_file, header = FALSE)
  K <- ncol(admix_q)
  
  # Format admixture data
  admix_data <- admix_q %>%
    mutate(sample_id = sample_ids) %>%
    pivot_longer(
      cols = starts_with("V"),
      names_to = "cluster",
      values_to = "proportion"
    ) %>%
    mutate(cluster = gsub("V", "K", cluster))
  
  # Admixture barplot
  p_admix <- ggplot(admix_data, aes(x = sample_id, y = proportion, fill = cluster)) +
    geom_bar(stat = "identity", width = 1) +
    labs(
      x = "Sample",
      y = "Ancestry Proportion",
      title = paste0("Admixture Analysis (K = ", K, ")"),
      fill = "Cluster"
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_classic(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.border = element_rect(fill = NA, color = "black"),
      legend.position = "right"
    )
  
  # Save admixture plot
  cat("Saving admixture plot...\n")
  ggsave(
    paste0(output_prefix, "_pcangsd_admixture.png"),
    p_admix,
    width = max(10, n_samples * 0.2),
    height = 6,
    dpi = 300
  )
}

# Read and summarize tree if available
tree_file <- paste0(output_prefix, ".pcangsd.tree")
if (file.exists(tree_file)) {
  cat("Tree file found:", tree_file, "\n")
  
  
  tree_data <- read.newick(tree_file) %>%
    as.treedata() %>%
    as_tibble() %>%
    mutate(label = sample_ids[as.integer(label)]) %>%
    as.treedata() 
  
  tree_plot <- tree_data %>%
    ggtree() +
    # layout_dendrogram() +
    geom_tippoint() +
    geom_tiplab(size = 6) +
    geom_treescale(fontsize = 6, linesize = 1, offset = 1) +
    theme(legend.position = 'bottom') +
    scale_x_continuous(expand = c(0,0.1))
  
  # Save admixture plot
  cat("Saving tree plot...\n")
  ggsave(
    paste0(output_prefix, "_pcangsd_tree.png"),
    tree_plot,
    width = 6,
    height = 6,
    dpi = 300
  )
}

cat("\nVisualization complete!\n")
cat("Output files:\n")
cat("  -", paste0(output_prefix, "_pcangsd_pca.png"), "\n")
if (file.exists(admix_q_file)) {
  cat("  -", paste0(output_prefix, "_pcangsd_admixture.png"), "\n")
}