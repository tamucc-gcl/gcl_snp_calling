#!/usr/bin/env Rscript

#### Parse command line arguments ####
# Args: <pca_subset.vcf.gz> <output_prefix> [ploidy_map] [dp_matrix.tsv.gz]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript vcf_qc_plots.R <pca_subset.vcf.gz|NO_PCA_VCF> <output_prefix> [ploidy_map] [dp_matrix.tsv.gz]")
}

pca_vcf_file  <- args[1]
output_prefix <- args[2]
ploidy_map    <- NULL

if (length(args) >= 3 && !args[3] %in% c("NO_PLOIDY_MAP", "NO_FILE")) {
  ploidy_map_file <- args[3]
} else {
  ploidy_map_file <- NULL
}

dp_matrix_file <- if (length(args) >= 4 && file.exists(args[4])) args[4] else NULL

have_pca_vcf <- pca_vcf_file != "NO_PCA_VCF" && file.exists(pca_vcf_file)

cat("Output prefix:", output_prefix, "\n")
cat("PCA VCF (small subset):", if (have_pca_vcf) pca_vcf_file else "not available", "\n")
cat("DP matrix:", dp_matrix_file %||% "not provided", "\n")

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

`%||%` <- function(a, b) if (!is.null(a)) a else b

#### Helper functions ####

safe_read_tsv <- function(file) {
  if (is.null(file) || !file.exists(file)) return(NULL)
  cat("Reading companion file:", file, "\n")
  tryCatch(
    readr::read_tsv(file, show_col_types = FALSE, progress = FALSE, comment = "#"),
    error = function(e) { cat("Failed:", e$message, "\n"); NULL }
  )
}

find_companion_file <- function(prefix, suffixes) {
  candidates <- unlist(lapply(suffixes, function(x)
    c(paste0(prefix, x), paste0(basename(prefix), x))))
  existing <- unique(candidates)[file.exists(unique(candidates))]
  if (length(existing) == 0) NULL else existing[1]
}

flag_outliers <- function(data, pca, alpha = 0.975) {
  if (nrow(pca$scores) < 3)
    return(mutate(data, outlier_score = NA_real_, outlier = FALSE,
                  robust_outlier_score = NA_real_, robust_outlier = FALSE))
  md  <- mahalanobis(pca$scores, colMeans(pca$scores), cov(pca$scores))
  mcd <- tryCatch(MASS::cov.mcd(pca$scores), error = function(e) NULL)
  if (is.null(mcd)) {
    md_robust   <- rep(NA_real_, nrow(pca$scores))
    robust_flag <- rep(FALSE,    nrow(pca$scores))
  } else {
    md_robust   <- mahalanobis(pca$scores, mcd$center, mcd$cov)
    robust_flag <- md_robust > qchisq(alpha, df = ncol(pca$scores))
  }
  mutate(data,
    outlier_score        = md,
    outlier              = md > qchisq(alpha, df = ncol(pca$scores)),
    robust_outlier_score = md_robust,
    robust_outlier       = robust_flag)
}

genind2genlight <- function(gi) {
  locna <- gi@loc.n.all
  ccc   <- 1L
  for (i in 2:length(locna))
    ccc[i] <- ccc[i-1] + if (locna[i-1] == 1L) 1L else 2L
  new("genlight", gi@tab[, ccc],
      pop = pop(gi), other = gi@other, ploidy = ploidy(gi),
      loc.names = locNames(gi), ind.names = indNames(gi))
}

extract_first_af <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == "."] <- NA_character_
  suppressWarnings(as.numeric(stringr::str_split_fixed(x, ",", 2)[, 1]))
}

# -----------------------------------------------------------------------
# Load DP matrix (pre-extracted by bcftools in the Nextflow process)
# Uses data.table for speed; falls back to readr if not available.
# -----------------------------------------------------------------------
load_dp_matrix <- function(path) {
  cat("Reading DP matrix:", path, "\n")
  if (requireNamespace("data.table", quietly = TRUE)) {
    tbl <- data.table::fread(path, sep = "\t", header = TRUE,
                              na.strings = c(".", "NA"), data.table = FALSE)
  } else {
    tbl <- readr::read_tsv(path, show_col_types = FALSE, progress = FALSE,
                            na = c(".", "NA"))
  }
  rnames <- paste(tbl[[1]], tbl[[2]], sep = "_")
  mat    <- as.matrix(tbl[, -(1:2), drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat)     <- rnames
  mat
}

# -----------------------------------------------------------------------
# Data loading
# -----------------------------------------------------------------------

# DP matrix — pre-extracted; no vcfR involved
if (!is.null(dp_matrix_file)) {
  depth_data <- load_dp_matrix(dp_matrix_file)
  depth_data[depth_data == 0] <- NA
} else {
  stop("No DP matrix provided. This script requires a pre-extracted DP matrix.")
}

sample_ids <- colnames(depth_data)
locus_ids  <- rownames(depth_data)
n_loci     <- length(locus_ids)
n_samples  <- length(sample_ids)
cat("Dimensions:", n_loci, "loci x", n_samples, "samples\n")

# Ploidy map
if (!is.null(ploidy_map_file)) {
  tryCatch({
    ploidy_map <- read.table(ploidy_map_file, header = FALSE,
                             col.names = c("sample", "ploidy"),
                             stringsAsFactors = FALSE, comment.char = "#")
    cat("Loaded ploidy info for", nrow(ploidy_map), "samples\n")
  }, error = function(e) {
    cat("Error reading ploidy map:", e$message, "\n")
    ploidy_map <<- NULL
  })
}

# PCA VCF — small pre-filtered subset only; vcfR only touches this file
if (have_pca_vcf) {
  cat("Loading PCA subset VCF via vcfR:", pca_vcf_file, "\n")
  pca_vcf <- read.vcfR(pca_vcf_file, verbose = FALSE)
  cat("PCA VCF loaded:", nrow(pca_vcf@fix), "sites\n")
} else {
  pca_vcf <- NULL
  cat("No PCA VCF available — PCA will be skipped\n")
}

#### Optional companion files ####
site_qc_file     <- find_companion_file(output_prefix, c(".site_qc.tsv.gz", ".site_qc.tsv"))
sample_qc_file   <- find_companion_file(output_prefix, c(".sample_qc.tsv.gz", ".sample_qc.tsv"))
freq_file        <- find_companion_file(output_prefix, c(".freq.tsv"))
missing_indv_file <- find_companion_file(output_prefix, c(".missing_indv.tsv"))

site_qc_tbl      <- safe_read_tsv(site_qc_file)
sample_qc_tbl    <- safe_read_tsv(sample_qc_file)
freq_tbl         <- safe_read_tsv(freq_file)
missing_indv_tbl <- safe_read_tsv(missing_indv_file)

#### Core QC summaries from DP matrix ####
miss_mat <- is.na(depth_data)

fallback_locus_qc <- tibble(
  locus                     = locus_ids,
  pct_samples_missing       = rowMeans(miss_mat),
  number_samples_with_locus = rowSums(!miss_mat),
  mean_depth                = rowMeans(depth_data, na.rm = TRUE)
)

fallback_sample_qc <- tibble(
  sample            = sample_ids,
  pct_loci_missing  = colMeans(miss_mat),
  number_loci       = colSums(!miss_mat),
  mean_depth_called = colMeans(depth_data, na.rm = TRUE)
)

rm(miss_mat); gc()

#### Build locus QC table ####
if (!is.null(site_qc_tbl)) {
  cat("Using site_qc companion table\n")
  locus_qc <- site_qc_tbl %>%
    mutate(locus = paste(chromo, position, sep = "_"),
           af    = if ("AF" %in% names(.)) extract_first_af(AF) else NA_real_) %>%
    select(any_of(c("locus","chromo","position","REF","ALT","QUAL","FILTER","NS","DP","af"))) %>%
    left_join(fallback_locus_qc, by = "locus")
} else {
  locus_qc <- fallback_locus_qc %>%
    mutate(af = NA_real_)
}

locus_qc <- locus_qc %>%
  mutate(
    pct_samples_missing       = replace_na(pct_samples_missing, 1),
    number_samples_with_locus = replace_na(number_samples_with_locus, 0),
    mean_depth                = replace_na(mean_depth, 0),
    af                        = suppressWarnings(as.numeric(af)),
    maf                       = pmin(af, 1 - af)
  )

#### Build sample QC table ####
if (!is.null(sample_qc_tbl)) {
  cat("Using sample_qc companion table\n")
  sample_qc <- sample_qc_tbl %>%
    rename(sample = any_of(c("sample", "INDV"))) %>%
    mutate(
      pct_loci_missing  = if ("f_missing"     %in% names(.)) f_missing     else pct_loci_missing,
      number_loci       = if ("sites_called"  %in% names(.)) sites_called  else number_loci,
      mean_depth_called = if ("mean_dp_called" %in% names(.)) mean_dp_called else mean_depth_called
    ) %>%
    select(any_of(c("sample","sites_total","sites_called","sites_missing",
                    "pct_loci_missing","number_loci","mean_depth_called",
                    "het","hom_ref","hom_alt")))
} else {
  sample_qc <- fallback_sample_qc
}

if (!is.null(missing_indv_tbl) && all(c("INDV","F_MISS") %in% names(missing_indv_tbl))) {
  f_miss_vals <- suppressWarnings(as.numeric(missing_indv_tbl$F_MISS))
  if (any(!is.na(f_miss_vals))) {
    sample_qc <- sample_qc %>%
      left_join(missing_indv_tbl %>%
                  rename(sample = INDV, f_mv = F_MISS) %>%
                  mutate(f_mv = suppressWarnings(as.numeric(f_mv))),
                by = "sample") %>%
      mutate(pct_loci_missing = if_else(!is.na(f_mv), f_mv, pct_loci_missing)) %>%
      select(-f_mv)
  }
}

sample_qc <- sample_qc %>%
  mutate(
    pct_loci_missing  = replace_na(pct_loci_missing, 1),
    number_loci       = replace_na(number_loci, 0),
    mean_depth_called = replace_na(mean_depth_called, 0)
  ) %>%
  arrange(desc(pct_loci_missing), sample) %>%
  mutate(
    sample_rank       = row_number(),
    flag_high_missing = pct_loci_missing > 0.3,
    flag_low_loci     = number_loci < quantile(number_loci, 0.05, na.rm = TRUE),
    flag_low_depth    = mean_depth_called < quantile(
                          mean_depth_called[mean_depth_called > 0], 0.05, na.rm = TRUE)
  )

#### MAF from freq table ####
if (!is.null(freq_tbl) && all(c("CHROM","POS","{ALLELE:FREQ}") %in% names(freq_tbl))) {
  maf_tbl <- freq_tbl %>%
    mutate(
      locus            = paste(CHROM, POS, sep = "_"),
      allele_freq_list = stringr::str_split(`{ALLELE:FREQ}`, "\t"),
      maf_from_freq    = purrr::map_dbl(allele_freq_list, function(x) {
        freqs <- suppressWarnings(as.numeric(stringr::str_extract(x, "(?<=:)[0-9.]+")))
        freqs <- freqs[!is.na(freqs)]
        if (length(freqs) == 0) NA_real_ else min(freqs)
      })
    ) %>%
    select(locus, maf_from_freq)
  locus_qc <- locus_qc %>%
    left_join(maf_tbl, by = "locus") %>%
    mutate(maf = coalesce(maf_from_freq, maf)) %>%
    select(-maf_from_freq)
}

#### Save derived tables ####
readr::write_tsv(sample_qc, paste0(output_prefix, "_sample_qc_derived.tsv"))
readr::write_tsv(locus_qc,  paste0(output_prefix, "_locus_qc_derived.tsv"))

#### Thresholds ####
missing_thresholds     <- c(0.1, 0.2, 0.3)
shared_threshold       <- median(locus_qc$number_samples_with_locus, na.rm = TRUE)
sample_depth_threshold <- quantile(
  sample_qc$mean_depth_called[sample_qc$mean_depth_called > 0], 0.1, na.rm = TRUE)
locus_depth_low  <- quantile(locus_qc$mean_depth[locus_qc$mean_depth > 0], 0.10, na.rm = TRUE)
locus_depth_high <- quantile(locus_qc$mean_depth[locus_qc$mean_depth > 0], 0.95, na.rm = TRUE)

#### QC Plots ####

locus_missingness_plot <- locus_qc %>%
  ggplot(aes(x = pct_samples_missing)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = missing_thresholds, linetype = "dashed") +
  scale_x_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(x = "% of samples missing locus", y = "Number of Loci",
       subtitle = paste("Median:", percent(median(locus_qc$pct_samples_missing, na.rm=TRUE), accuracy=0.1),
                        "| >30%:", comma(sum(locus_qc$pct_samples_missing > 0.3, na.rm=TRUE))))

sample_missingness_plot <- sample_qc %>%
  ggplot(aes(x = pct_loci_missing)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = missing_thresholds, linetype = "dashed") +
  scale_x_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(x = "% of loci missing in sample", y = "Number of Samples",
       subtitle = paste("Median:", percent(median(sample_qc$pct_loci_missing, na.rm=TRUE), accuracy=0.1),
                        "| >30%:", comma(sum(sample_qc$pct_loci_missing > 0.3, na.rm=TRUE))))

# Depth histogram: sample at most 5M values to stay well under 2^31
depth_plot <- {
  dv <- as.vector(depth_data)
  dv <- dv[!is.na(dv) & dv > 0]
  sampled <- length(dv) > 5e6
  if (sampled) dv <- sample(dv, 5e6)
  data.frame(depth = dv) %>%
    ggplot(aes(x = depth)) +
    geom_histogram(bins = 60) +
    scale_x_continuous(transform = scales::log10_trans(),
                       labels = scales::comma_format(),
                       guide  = guide_axis_logticks()) +
    scale_y_continuous(transform = scales::log10_trans(),
                       labels = scales::comma_format(),
                       guide  = guide_axis_logticks()) +
    labs(x = "Sequencing Depth", y = "Number of Genotypes",
         subtitle = if (sampled) "Called genotypes (5M sample)" else "Called genotypes only")
}

meanDepth_plot <- locus_qc %>%
  filter(mean_depth > 0) %>%
  ggplot(aes(x = mean_depth)) +
  geom_histogram(bins = 60) +
  geom_vline(xintercept = c(locus_depth_low, locus_depth_high), linetype = "dashed") +
  scale_x_continuous(transform = scales::log10_trans(),
                     labels = scales::comma_format(),
                     guide  = guide_axis_logticks()) +
  scale_y_continuous(transform = scales::log10_trans(),
                     labels = scales::comma_format(),
                     guide  = guide_axis_logticks()) +
  labs(x = "Mean Sequencing Depth", y = "Number of Loci",
       subtitle = paste("Low / high reference:", round(locus_depth_low,1), "/", round(locus_depth_high,1)))

distribution_shared_loci <- locus_qc %>%
  ggplot(aes(x = number_samples_with_locus)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = shared_threshold, linetype = "dashed") +
  scale_x_continuous(labels = scales::comma_format()) +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(x = "Number of Samples with Locus", y = "Number of Loci",
       subtitle = paste("Median shared samples:", comma(round(shared_threshold))))

ranked_sample_qc_plot <- sample_qc %>%
  arrange(number_loci, sample) %>%
  mutate(rank = row_number()) %>%
  ggplot(aes(y = rank, x = number_loci)) +
  geom_segment(aes(yend = rank, x = 0, xend = number_loci), linewidth = 0.3) +
  geom_point(aes(color = flag_high_missing), size = 1.5, show.legend = FALSE) +
  scale_color_manual(values = c(`TRUE` = "firebrick", `FALSE` = "black")) +
  scale_y_continuous(labels = scales::comma_format()) +
  scale_x_continuous(labels = scales::comma_format()) +
  labs(y = "Samples ranked by called loci", x = "Number of Loci",
       subtitle = "Worst samples flagged by >30% missing loci")

snp_summary_stats_plot <-
  locus_missingness_plot + sample_missingness_plot +
  depth_plot + meanDepth_plot +
  distribution_shared_loci + ranked_sample_qc_plot +
  plot_layout(ncol = 2) + plot_annotation(title = "Raw SNP Summary Plots") &
  theme_classic(base_size = 16) + theme(panel.background = element_rect(colour = "black"))

ggsave(paste0(output_prefix, "_summary_plots.png"), snp_summary_stats_plot, height=15, width=10)

#### Additional QC plots ####

sample_mean_depth_plot <- sample_qc %>%
  arrange(mean_depth_called) %>%
  mutate(rank = row_number()) %>%
  ggplot(aes(x = rank, y = mean_depth_called)) +
  geom_segment(aes(xend = rank, y = 0, yend = mean_depth_called), linewidth = 0.3) +
  geom_point(aes(color = flag_low_depth), size = 1.5, show.legend = FALSE) +
  scale_color_manual(values = c(`TRUE` = "firebrick", `FALSE` = "black")) +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(x = "Samples ranked by mean called depth", y = "Mean Called Depth",
       subtitle = paste("10th percentile reference:", round(sample_depth_threshold, 2)))

maf_plot <- locus_qc %>%
  filter(!is.na(maf), maf >= 0, maf <= 0.5) %>%
  ggplot(aes(x = maf)) +
  geom_histogram(bins = 60) +
  geom_vline(xintercept = c(0.01, 0.05, 0.1), linetype = "dashed") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.5)) +
  scale_y_continuous(labels = scales::comma_format()) +
  labs(x = "Minor Allele Frequency", y = "Number of Loci",
       subtitle = "From site AF when available")

# Hex bin heatmap — no sampling needed, handles 16M+ loci without overplotting
locus_depth_missingness_plot <- locus_qc %>%
  filter(mean_depth > 0) %>%
  ggplot(aes(x = pct_samples_missing, y = mean_depth)) +
  geom_hex(bins = 80) +
  scale_fill_viridis_c(
    name   = "Loci\n(log10)",
    trans  = "log10",
    labels = scales::comma_format(),
    option = "magma"
  ) +
  geom_vline(xintercept = missing_thresholds, linetype = "dashed", colour = "grey40") +
  scale_x_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  scale_y_continuous(transform = scales::log10_trans(),
                     labels = scales::comma_format(),
                     guide  = guide_axis_logticks()) +
  labs(x = "% of samples missing locus", y = "Mean Sequencing Depth",
       subtitle = "Hex-bin density (log10 count fill)")

advanced_summary_plot <-
  sample_mean_depth_plot + maf_plot + locus_depth_missingness_plot + plot_spacer() +
  plot_layout(ncol = 2) + plot_annotation(title = "Additional Raw SNP QC Plots") &
  theme_classic(base_size = 16) + theme(panel.background = element_rect(colour = "black"))

ggsave(paste0(output_prefix, "_summary_plots_extra.png"), advanced_summary_plot, height=10, width=10)

#### PCA — only runs if PCA VCF was provided ####
if (!is.null(pca_vcf)) {

  if (!is.null(ploidy_map)) {
    raw_genlight <- ploidy_map %>%
      arrange(match(sample, colnames(pca_vcf@gt)[-1])) %>%
      pull(ploidy) %>%
      vcfR2genind(pca_vcf, ploidy = .) %>%
      genind2genlight()
  } else {
    raw_genlight <- vcfR2genlight(pca_vcf)
    ploidy(raw_genlight) <- max(ploidy(raw_genlight), na.rm = TRUE)
  }

  to_remove <- is.na(glMean(raw_genlight, alleleAsUnit = FALSE))
  if (any(to_remove)) {
    raw_genlight <- raw_genlight[, !to_remove]
    cat("Removed", sum(to_remove), "all-NA loci before PCA\n")
  }

  pca_success <- TRUE
  pca_out <- tryCatch(
    glPca(raw_genlight, nf = 2, loadings = FALSE),
    error = function(e) { cat("PCA failed:", e$message, "\n"); pca_success <<- FALSE; NULL }
  )

  if (pca_success && !is.null(pca_out)) {
    pct_var <- scales::percent(pca_out$eig / sum(pca_out$eig), accuracy = 0.1)

    raw_pca_scores <- pca_out$scores %>%
      as_tibble(rownames = "sample") %>%
      flag_outliers(pca_out) %>%
      left_join(sample_qc, by = "sample")

    readr::write_tsv(raw_pca_scores, paste0(output_prefix, "_pca_scores.tsv"))

    use_text_labels <- !is.null(ploidy_map) && is.data.frame(ploidy_map) &&
                       nrow(ploidy_map) > 0 && nrow(ploidy_map) < 10

    raw_snp_pca <- raw_pca_scores %>%
      ggplot(aes(x = PC1, y = PC2, color = pct_loci_missing, size = mean_depth_called)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dashed")

    if (use_text_labels) {
      raw_snp_pca <- raw_snp_pca + geom_text(aes(label = sample))
    } else {
      raw_snp_pca <- raw_snp_pca + geom_point(alpha = 0.85)
      label_data  <- raw_pca_scores %>%
        filter(outlier | robust_outlier | flag_high_missing | flag_low_depth)
      if (nrow(label_data) > 0)
        raw_snp_pca <- raw_snp_pca +
          geom_text_repel(data = label_data, aes(label = sample),
                          hjust = "inward", vjust = "inward", max.overlaps = 50)
    }

    raw_snp_pca <- raw_snp_pca +
      scale_color_viridis_c(labels = scales::percent_format(), option = "D", na.value = "grey50") +
      scale_size_continuous(range = c(1.5, 5)) +
      labs(x = str_c("PC1 (", pct_var[1], ")"), y = str_c("PC2 (", pct_var[2], ")"),
           color = "% missing", size = "Mean\ncalled depth", title = "Raw SNP PCA",
           subtitle = "Colored by sample missingness; labels flag outliers / poor-QC samples") +
      theme_classic(base_size = 16) +
      theme(panel.background = element_rect(colour = "black"))

    ggsave(paste0(output_prefix, "_pca.png"), raw_snp_pca, height = 5, width = 6)

  } else {
    ggplot() +
      annotate("text", x=0.5, y=0.5, label="PCA could not be computed", size=8) +
      theme_void() + theme(panel.border = element_rect(fill = NA)) %>%
      ggsave(paste0(output_prefix, "_pca.png"), ., height=5, width=6)
  }

} else {
  ggplot() +
    annotate("text", x=0.5, y=0.5, label="PCA not available\n(no sites passed MAF filter)", size=8) +
    theme_void() + theme(panel.border = element_rect(fill = NA)) %>%
    ggsave(paste0(output_prefix, "_pca.png"), ., height=5, width=6)
}

#### Flagged summary outputs ####
readr::write_tsv(
  sample_qc %>% arrange(desc(pct_loci_missing), mean_depth_called) %>% slice_head(n = 20),
  paste0(output_prefix, "_worst_samples.tsv"))
readr::write_tsv(
  locus_qc  %>% arrange(desc(pct_samples_missing), desc(mean_depth)) %>% slice_head(n = 20),
  paste0(output_prefix, "_worst_loci.tsv"))

cat("QC plots complete!\n")