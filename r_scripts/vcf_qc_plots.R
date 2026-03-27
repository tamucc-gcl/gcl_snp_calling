#!/usr/bin/env Rscript

#### Parse command line arguments ####
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript vcf_qc_plots.R <input.vcf.gz> <output_prefix> [ploidy_map]")
}

vcf_file <- args[1]
output_prefix <- args[2]

# Initialize ploidy_map to NULL first to ensure it always exists
ploidy_map <- NULL

# Check if a ploidy map file was provided
if (length(args) >= 3 && args[3] != "NO_PLOIDY_MAP" && args[3] != "NO_FILE") {
  ploidy_map_file <- args[3]
} else {
  ploidy_map_file <- NULL
}

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
library(readr)
library(scales)
library(purrr)
library(MASS)
select <- dplyr::select

#### Helper functions ####

safe_read_tsv <- function(file) {
  if (is.null(file) || !file.exists(file)) {
    return(NULL)
  }
  
  cat("Reading companion file:", file, "\n")
  
  tryCatch({
    readr::read_tsv(
      file,
      show_col_types = FALSE,
      progress = FALSE,
      comment = "#"
    )
  }, error = function(e) {
    cat("Failed to read", file, ":", e$message, "\n")
    NULL
  })
}

find_companion_file <- function(prefix, suffixes) {
  candidate_files <- unlist(lapply(suffixes, function(x) {
    c(
      paste0(prefix, x),
      paste0(basename(prefix), x)
    )
  }))
  
  candidate_files <- unique(candidate_files)
  existing <- candidate_files[file.exists(candidate_files)]
  
  if (length(existing) == 0) {
    return(NULL)
  }
  
  existing[1]
}

calc_maf_from_gt <- function(gt_vec) {
  gt_vec <- gt_vec[!is.na(gt_vec)]
  gt_vec <- gt_vec[!gt_vec %in% c("./.", ".|.")]
  
  if (length(gt_vec) == 0) {
    return(NA_real_)
  }
  
  alleles <- unlist(strsplit(gt_vec, "[/|]"))
  alleles <- alleles[alleles != "."]
  
  if (length(alleles) == 0) {
    return(NA_real_)
  }
  
  alt_count <- sum(alleles != "0")
  af <- alt_count / length(alleles)
  pmin(af, 1 - af)
}

flag_outliers <- function(data, pca, alpha = 0.975) {
  if (nrow(pca$scores) < 3) {
    return(
      mutate(
        data,
        outlier_score = NA_real_,
        outlier = FALSE,
        robust_outlier_score = NA_real_,
        robust_outlier = FALSE
      )
    )
  }
  
  md <- mahalanobis(pca$scores, colMeans(pca$scores), cov(pca$scores))
  
  mcd <- tryCatch({
    MASS::cov.mcd(pca$scores)
  }, error = function(e) {
    NULL
  })
  
  if (is.null(mcd)) {
    md_robust <- rep(NA_real_, nrow(pca$scores))
    robust_flag <- rep(FALSE, nrow(pca$scores))
  } else {
    md_robust <- mahalanobis(pca$scores, mcd$center, mcd$cov)
    robust_flag <- md_robust > qchisq(alpha, df = ncol(pca$scores))
  }
  
  mutate(
    data,
    outlier_score = md,
    outlier = md > qchisq(alpha, df = ncol(pca$scores)),
    robust_outlier_score = md_robust,
    robust_outlier = robust_flag
  )
}

genind2genlight <- function(gi) {
  locna <- gi@loc.n.all
  ccc <- 1
  for (i in 2:length(locna)) {
    if (locna[i - 1] == 1) {
      ccc[i] <- ccc[i - 1] + 1
    } else {
      ccc[i] <- ccc[i - 1] + 2
    }
  }
  
  new(
    "genlight",
    gi@tab[, ccc],
    pop = pop(gi),
    other = gi@other,
    ploidy = ploidy(gi),
    loc.names = locNames(gi),
    ind.names = indNames(gi)
  )
}

extract_first_af <- function(x) {
  if (is.na(x) || x == ".") {
    return(NA_real_)
  }
  
  x <- str_split(as.character(x), ",", simplify = TRUE)[1]
  
  suppressWarnings(as.numeric(x))
}

make_threshold_label <- function(x, digits = 1) {
  paste0(percent(x, accuracy = 0.1), " cutoff")
}

#### Data ####
raw_vcf <- read.vcfR(vcf_file)

if (!is.null(ploidy_map_file) && ploidy_map_file != "NO_FILE" && ploidy_map_file != "NO_PLOIDY_MAP") {
  cat("Attempting to read ploidy map file:", ploidy_map_file, "\n")
  
  tryCatch({
    ploidy_map <- read.table(
      ploidy_map_file,
      header = FALSE,
      col.names = c("sample", "ploidy"),
      stringsAsFactors = FALSE,
      comment.char = "#"
    )
    
    cat("Successfully loaded ploidy information for", nrow(ploidy_map), "samples\n")
    cat("Ploidy map class:", class(ploidy_map), "\n")
    cat("Ploidy map dimensions:", dim(ploidy_map), "\n")
  }, error = function(e) {
    cat("Error reading ploidy map file:", e$message, "\n")
    cat("Setting ploidy_map to NULL\n")
    ploidy_map <<- NULL
  })
  
} else {
  cat("No ploidy map provided (file:", ploidy_map_file, ")\n")
  cat("Setting ploidy_map to NULL\n")
  ploidy_map <- NULL
}

# Verify the state of ploidy_map
cat("Final ploidy_map status:\n")
cat("  - is.null(ploidy_map):", is.null(ploidy_map), "\n")
if (!is.null(ploidy_map)) {
  cat("  - class(ploidy_map):", class(ploidy_map), "\n")
  cat("  - nrow(ploidy_map):", nrow(ploidy_map), "\n")
}

#### Optional companion files from summarize_vcfs.nf ####
site_qc_file <- find_companion_file(output_prefix, c(".site_qc.tsv.gz", ".site_qc.tsv"))
sample_qc_file <- find_companion_file(output_prefix, c(".sample_qc.tsv.gz", ".sample_qc.tsv"))
freq_file <- find_companion_file(output_prefix, c(".freq.tsv"))
missing_site_file <- find_companion_file(output_prefix, c(".missing_site.tsv"))
missing_indv_file <- find_companion_file(output_prefix, c(".missing_indv.tsv"))

site_qc_tbl <- safe_read_tsv(site_qc_file)
sample_qc_tbl <- safe_read_tsv(sample_qc_file)
freq_tbl <- safe_read_tsv(freq_file)
missing_site_tbl <- safe_read_tsv(missing_site_file)
missing_indv_tbl <- safe_read_tsv(missing_indv_file)

#### Derive core matrices from VCF ####
depth_data <- extract.gt(
  raw_vcf,
  element = "DP",
  as.numeric = TRUE
)

gt_data <- extract.gt(
  raw_vcf,
  element = "GT",
  as.numeric = FALSE
)

depth_data[depth_data == 0] <- NA

sample_ids <- colnames(depth_data)
locus_ids <- rownames(depth_data)

#### Derive fallback QC summaries directly from VCF ####
fallback_locus_qc <- tibble(
  locus = locus_ids,
  pct_samples_missing = is.na(depth_data) %>% rowMeans(),
  number_samples_with_locus = (!is.na(depth_data)) %>% rowSums(),
  mean_depth = rowMeans(depth_data, na.rm = TRUE)
)

fallback_sample_qc <- tibble(
  sample = sample_ids,
  pct_loci_missing = is.na(depth_data) %>% colMeans(),
  number_loci = (!is.na(depth_data)) %>% colSums(),
  mean_depth_called = apply(depth_data, 2, function(x) mean(x, na.rm = TRUE))
)

#### Build locus QC table ####
if (!is.null(site_qc_tbl)) {
  cat("Using site_qc companion table for locus metrics where possible\n")
  
  locus_qc <- site_qc_tbl %>%
    mutate(
      locus = paste(chromo, position, sep = ":"),
      af = if ("AF" %in% names(.)) extract_first_af(AF) else NA_real_
    ) %>%
    select(any_of(c("locus", "chromo", "position", "REF", "ALT", "QUAL", "FILTER", "NS", "DP", "af"))) %>%
    left_join(fallback_locus_qc, by = "locus")
  
} else {
  cat("No site_qc companion table found; deriving locus metrics from VCF only\n")
  
  locus_qc <- fallback_locus_qc %>%
    mutate(
      chromo = raw_vcf@fix[, "CHROM"],
      position = as.integer(raw_vcf@fix[, "POS"]),
      REF = raw_vcf@fix[, "REF"],
      ALT = raw_vcf@fix[, "ALT"],
      QUAL = raw_vcf@fix[, "QUAL"],
      FILTER = raw_vcf@fix[, "FILTER"],
      NS = number_samples_with_locus,
      DP = rowSums(depth_data, na.rm = TRUE),
      af = apply(gt_data, 1, calc_maf_from_gt)
    )
}

locus_qc <- locus_qc %>%
  mutate(
    pct_samples_missing = replace_na(pct_samples_missing, 1),
    number_samples_with_locus = replace_na(number_samples_with_locus, 0),
    mean_depth = replace_na(mean_depth, 0),
    af = suppressWarnings(as.numeric(af)),
    maf = pmin(af, 1 - af)
  )

#### Build sample QC table ####
if (!is.null(sample_qc_tbl)) {
  cat("Using sample_qc companion table for sample metrics where possible\n")
  
  sample_qc <- sample_qc_tbl %>%
    rename(sample = any_of(c("sample", "INDV"))) %>%
    mutate(
      pct_loci_missing = if ("f_missing" %in% names(.)) f_missing else pct_loci_missing,
      number_loci = if ("sites_called" %in% names(.)) sites_called else number_loci,
      mean_depth_called = if ("mean_dp_called" %in% names(.)) mean_dp_called else mean_depth_called
    ) %>%
    select(any_of(c(
      "sample", "sites_total", "sites_called", "sites_missing",
      "pct_loci_missing", "number_loci", "mean_depth_called",
      "het", "hom_ref", "hom_alt"
    )))
  
} else {
  cat("No sample_qc companion table found; deriving sample metrics from VCF only\n")
  
  sample_qc <- fallback_sample_qc
}

if (!is.null(missing_indv_tbl)) {
  sample_qc <- sample_qc %>%
    left_join(
      missing_indv_tbl %>%
        rename(sample = INDV, f_missing_vcftools = F_MISS),
      by = "sample"
    ) %>%
    mutate(
      pct_loci_missing = if_else(
        !is.na(f_missing_vcftools),
        f_missing_vcftools,
        pct_loci_missing
      )
    ) %>%
    select(-f_missing_vcftools)
}

sample_qc <- sample_qc %>%
  mutate(
    pct_loci_missing = replace_na(pct_loci_missing, 1),
    number_loci = replace_na(number_loci, 0),
    mean_depth_called = replace_na(mean_depth_called, 0)
  ) %>%
  arrange(desc(pct_loci_missing), sample) %>%
  mutate(
    sample_rank = row_number(),
    flag_high_missing = pct_loci_missing > 0.3,
    flag_low_loci = number_loci < quantile(number_loci, 0.05, na.rm = TRUE),
    flag_low_depth = mean_depth_called < quantile(mean_depth_called[mean_depth_called > 0], 0.05, na.rm = TRUE)
  )

#### MAF table ####
if (!is.null(freq_tbl)) {
  cat("Attempting to parse MAF from frequency table\n")
  
  if (all(c("CHROM", "POS") %in% names(freq_tbl))) {
    freq_long <- freq_tbl %>%
      mutate(
        locus = paste(CHROM, POS, sep = ":")
      )
    
    allele_cols <- setdiff(names(freq_long), c("CHROM", "POS", "N_ALLELES", "N_CHR", "locus"))
    
    if (length(allele_cols) > 0) {
      maf_tbl <- freq_long %>%
        pivot_longer(
          cols = all_of(allele_cols),
          names_to = "allele_col",
          values_to = "allele_freq_raw"
        ) %>%
        mutate(
          allele_freq = str_extract(allele_freq_raw, "[0-9.eE+-]+$"),
          allele_freq = suppressWarnings(as.numeric(allele_freq))
        ) %>%
        filter(!is.na(allele_freq)) %>%
        group_by(locus) %>%
        summarise(
          maf_from_freq = min(allele_freq, 1 - allele_freq, na.rm = TRUE),
          .groups = "drop"
        )
      
      locus_qc <- locus_qc %>%
        left_join(maf_tbl, by = "locus") %>%
        mutate(maf = coalesce(maf_from_freq, maf)) %>%
        select(-maf_from_freq)
    }
  }
}

#### Save derived tables ####
readr::write_tsv(sample_qc, paste0(output_prefix, "_sample_qc_derived.tsv"))
readr::write_tsv(locus_qc, paste0(output_prefix, "_locus_qc_derived.tsv"))

#### Thresholds used for interpretation ####
missing_thresholds <- c(0.1, 0.2, 0.3)
shared_threshold <- median(locus_qc$number_samples_with_locus, na.rm = TRUE)
sample_depth_threshold <- quantile(sample_qc$mean_depth_called[sample_qc$mean_depth_called > 0], 0.1, na.rm = TRUE)
locus_depth_low <- quantile(locus_qc$mean_depth[locus_qc$mean_depth > 0], 0.1, na.rm = TRUE)
locus_depth_high <- quantile(locus_qc$mean_depth[locus_qc$mean_depth > 0], 0.95, na.rm = TRUE)

#### QC Plots ####

locus_missingness_plot <- locus_qc %>%
  ggplot(aes(x = pct_samples_missing)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = missing_thresholds,
             linetype = "dashed") +
  scale_x_continuous(
    labels = scales::percent_format(),
    limits = c(0, 1)
  ) +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(
    x = "% of samples missing locus",
    y = "Number of Loci",
    subtitle = paste(
      "Median:",
      percent(median(locus_qc$pct_samples_missing, na.rm = TRUE), accuracy = 0.1),
      "| >30%:",
      comma(sum(locus_qc$pct_samples_missing > 0.3, na.rm = TRUE))
    )
  )

sample_missingness_plot <- sample_qc %>%
  ggplot(aes(x = pct_loci_missing)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = missing_thresholds,
             linetype = "dashed") +
  scale_x_continuous(
    labels = scales::percent_format(),
    limits = c(0, 1)
  ) +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(
    x = "% of loci missing in sample",
    y = "Number of Samples",
    subtitle = paste(
      "Median:",
      percent(median(sample_qc$pct_loci_missing, na.rm = TRUE), accuracy = 0.1),
      "| >30%:",
      comma(sum(sample_qc$pct_loci_missing > 0.3, na.rm = TRUE))
    )
  )

depth_plot <- depth_data %>%
  as.numeric() %>%
  tibble(depth = .) %>%
  mutate(depth = replace_na(depth, 0)) %>%
  filter(depth > 0) %>%
  ggplot(aes(x = depth)) +
  geom_histogram(bins = 60) +
  scale_x_continuous(
    transform = scales::log10_trans(),
    labels = scales::comma_format(),
    guide = guide_axis_logticks()
  ) +
  scale_y_continuous(
    transform = scales::log10_trans(),
    labels = scales::comma_format(),
    guide = guide_axis_logticks()
  ) +
  labs(
    x = "Sequencing Depth",
    y = "Number of Genotypes",
    subtitle = "Called genotypes only"
  )

meanDepth_plot <- locus_qc %>%
  filter(mean_depth > 0) %>%
  ggplot(aes(x = mean_depth)) +
  geom_histogram(bins = 60) +
  geom_vline(
    xintercept = c(locus_depth_low, locus_depth_high),
    linetype = "dashed"
  ) +
  scale_x_continuous(
    transform = scales::log10_trans(),
    labels = scales::comma_format(),
    guide = guide_axis_logticks()
  ) +
  scale_y_continuous(
    transform = scales::log10_trans(),
    labels = scales::comma_format(),
    guide = guide_axis_logticks()
  ) +
  labs(
    x = "Mean Sequencing Depth",
    y = "Number of Loci",
    subtitle = paste(
      "Low / high reference:",
      round(locus_depth_low, 1),
      "/",
      round(locus_depth_high, 1)
    )
  )

distribution_shared_loci <- locus_qc %>%
  ggplot(aes(x = number_samples_with_locus)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = shared_threshold,
             linetype = "dashed") +
  scale_x_continuous(labels = scales::comma_format()) +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(
    x = "Number of Samples with Locus",
    y = "Number of Loci",
    subtitle = paste(
      "Median shared samples:",
      comma(round(shared_threshold))
    )
  )

ranked_sample_qc_plot <- sample_qc %>%
  arrange(desc(number_loci), sample) %>%
  mutate(rank = row_number()) %>%
  ggplot(aes(x = rank, y = number_loci)) +
  geom_segment(aes(xend = rank, y = 0, yend = number_loci),
               linewidth = 0.3) +
  geom_point(aes(color = flag_high_missing), size = 1.5) +
  scale_color_manual(values = c(`TRUE` = "firebrick", `FALSE` = "black")) +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(
    x = "Samples ranked by called loci",
    y = "Number of Loci",
    color = "High missingness",
    subtitle = "Worst samples flagged by >30% missing loci"
  )

snp_summary_stats_plot <- locus_missingness_plot +
  sample_missingness_plot +
  depth_plot +
  meanDepth_plot +
  distribution_shared_loci +
  ranked_sample_qc_plot +
  plot_layout(ncol = 2) +
  plot_annotation(title = "Raw SNP Summary Plots") &
  theme_classic(base_size = 16) +
  theme(panel.background = element_rect(colour = "black"))

ggsave(
  paste0(output_prefix, "_summary_plots.png"),
  plot = snp_summary_stats_plot,
  height = 15,
  width = 10
)

#### Additional QC plots ####

sample_mean_depth_plot <- sample_qc %>%
  arrange(mean_depth_called) %>%
  mutate(rank = row_number()) %>%
  ggplot(aes(x = rank, y = mean_depth_called)) +
  geom_segment(aes(xend = rank, y = 0, yend = mean_depth_called),
               linewidth = 0.3) +
  geom_point(aes(color = flag_low_depth), size = 1.5) +
  scale_color_manual(values = c(`TRUE` = "firebrick", `FALSE` = "black")) +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(
    x = "Samples ranked by mean called depth",
    y = "Mean Called Depth",
    color = "Low depth",
    subtitle = paste("10th percentile reference:", round(sample_depth_threshold, 2))
  )

maf_plot <- locus_qc %>%
  filter(!is.na(maf), maf >= 0, maf <= 0.5) %>%
  ggplot(aes(x = maf)) +
  geom_histogram(bins = 60) +
  geom_vline(xintercept = c(0.01, 0.05, 0.1),
             linetype = "dashed") +
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01),
                     limits = c(0, 0.5)) +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(
    x = "Minor Allele Frequency",
    y = "Number of Loci",
    subtitle = "From site AF when available"
  )

locus_depth_missingness_plot <- locus_qc %>%
  filter(mean_depth > 0) %>%
  ggplot(aes(x = pct_samples_missing, y = mean_depth)) +
  geom_point(alpha = 0.15, size = 0.75) +
  geom_vline(xintercept = missing_thresholds,
             linetype = "dashed") +
  scale_x_continuous(labels = scales::percent_format(),
                     limits = c(0, 1)) +
  scale_y_continuous(
    transform = scales::log10_trans(),
    labels = scales::comma_format(),
    guide = guide_axis_logticks()
  ) +
  labs(
    x = "% of samples missing locus",
    y = "Mean Sequencing Depth",
    subtitle = "Locus depth versus locus missingness"
  )

advanced_summary_plot <- sample_mean_depth_plot +
  maf_plot +
  locus_depth_missingness_plot +
  plot_spacer() +
  plot_layout(ncol = 2) +
  plot_annotation(title = "Additional Raw SNP QC Plots") &
  theme_classic(base_size = 16) +
  theme(panel.background = element_rect(colour = "black"))

ggsave(
  paste0(output_prefix, "_summary_plots_extra.png"),
  plot = advanced_summary_plot,
  height = 10,
  width = 10
)

#### PCA ####

# Build genlight object
if (is.null(ploidy_map)) {
  raw_genlight <- vcfR2genlight(raw_vcf)
  ploidy(raw_genlight) <- max(ploidy(raw_genlight), na.rm = TRUE)
  
} else {
  raw_genlight <- ploidy_map %>%
    arrange(match(sample, colnames(raw_vcf@gt)[-1])) %>%
    pull(ploidy) %>%
    vcfR2genind(raw_vcf, ploidy = .) %>%
    genind2genlight()
}

# Remove loci that are entirely NA in the genlight representation
to_remove <- is.na(glMean(raw_genlight, alleleAsUnit = FALSE))
if (any(to_remove)) {
  raw_genlight <- raw_genlight[, !to_remove]
  cat("Removed", sum(to_remove), "loci with all-NA genlight means before PCA\n")
}

pca_success <- TRUE
pca_out <- tryCatch({
  glPca(raw_genlight, nf = 2, loadings = FALSE)
}, error = function(e) {
  cat("PCA failed:", e$message, "\n")
  pca_success <<- FALSE
  NULL
})

if (pca_success && !is.null(pca_out)) {
  pct_var <- scales::percent(pca_out$eig / sum(pca_out$eig), accuracy = 0.1)
  
  raw_pca_scores <- pca_out$scores %>%
    as_tibble(rownames = "sample") %>%
    flag_outliers(pca_out) %>%
    left_join(sample_qc, by = "sample")
  
  readr::write_tsv(raw_pca_scores, paste0(output_prefix, "_pca_scores.tsv"))
  
  use_text_labels <- FALSE
  
  # Extremely defensive check for ploidy_map / sample count
  if (!is.null(ploidy_map)) {
    cat("Checking ploidy_map for PCA label decision...\n")
    if (is.data.frame(ploidy_map)) {
      if (nrow(ploidy_map) > 0) {
        if (nrow(ploidy_map) < 10) {
          use_text_labels <- TRUE
          cat("Will use text labels (small dataset with", nrow(ploidy_map), "samples)\n")
        }
      }
    }
  }
  
  raw_snp_pca <- raw_pca_scores %>%
    ggplot(aes(
      x = PC1,
      y = PC2,
      color = pct_loci_missing,
      size = mean_depth_called
    )) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed")
  
  if (use_text_labels) {
    raw_snp_pca <- raw_snp_pca +
      geom_text(aes(label = sample))
  } else {
    raw_snp_pca <- raw_snp_pca +
      geom_point(alpha = 0.85) 
    
    pca_label_data <- raw_pca_scores %>%
      filter(outlier | robust_outlier | flag_high_missing | flag_low_depth)
    
    if (nrow(pca_label_data) > 0) {
      raw_snp_pca <- raw_snp_pca +
        geom_text_repel(
          data = pca_label_data,
          aes(label = sample),
          hjust = "inward",
          vjust = "inward",
          max.overlaps = 50
        )
    }
  }
  
  raw_snp_pca <- raw_snp_pca +
    scale_color_viridis_c(
      labels = scales::percent_format(),
      option = "D",
      na.value = "grey50"
    ) +
    scale_size_continuous(range = c(1.5, 5)) +
    labs(
      x = str_c("PC1 (", pct_var[1], ")"),
      y = str_c("PC2 (", pct_var[2], ")"),
      color = "% missing",
      size = "Mean\ncalled depth",
      title = "Raw SNP PCA",
      subtitle = "Colored by sample missingness; labels flag outliers / poor-QC samples"
    ) +
    theme_classic(base_size = 16) +
    theme(panel.background = element_rect(colour = "black"))
  
  ggsave(
    paste0(output_prefix, "_pca.png"),
    plot = raw_snp_pca,
    height = 5,
    width = 6
  )
  
} else {
  cat("Creating placeholder PCA plot because PCA did not succeed\n")
  
  placeholder <- ggplot() +
    annotate(
      "text",
      x = 0.5, y = 0.5,
      label = "PCA could not be computed",
      size = 8
    ) +
    theme_void() +
    theme(panel.border = element_rect(fill = NA))
  
  ggsave(
    paste0(output_prefix, "_pca.png"),
    plot = placeholder,
    height = 5,
    width = 6
  )
}

#### Flagged summary outputs ####

worst_samples_tbl <- sample_qc %>%
  arrange(desc(pct_loci_missing), mean_depth_called) %>%
  slice_head(n = 20)

worst_loci_tbl <- locus_qc %>%
  arrange(desc(pct_samples_missing), desc(mean_depth)) %>%
  slice_head(n = 20)

readr::write_tsv(worst_samples_tbl, paste0(output_prefix, "_worst_samples.tsv"))
readr::write_tsv(worst_loci_tbl, paste0(output_prefix, "_worst_loci.tsv"))

cat("QC plots successfully generated!\n")
cat("Wrote derived sample QC table:", paste0(output_prefix, "_sample_qc_derived.tsv"), "\n")
cat("Wrote derived locus QC table:", paste0(output_prefix, "_locus_qc_derived.tsv"), "\n")
cat("Wrote extra summary plot:", paste0(output_prefix, "_summary_plots_extra.png"), "\n")
cat("Wrote worst-sample table:", paste0(output_prefix, "_worst_samples.tsv"), "\n")
cat("Wrote worst-locus table:", paste0(output_prefix, "_worst_loci.tsv"), "\n")