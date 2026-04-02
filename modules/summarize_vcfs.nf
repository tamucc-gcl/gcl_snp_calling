process SUMMARIZE_VCFS {
    tag "${vcf.baseName}"
    publishDir "${params.output_dir}/snp_qc", mode: 'copy'

    input:
    path vcf
    path ploidy_map
    val caller

    output:
    tuple val(vcf.simpleName),
          path("${vcf.simpleName}.stats.txt"),
          path("${vcf.simpleName}.site_qc.tsv.gz"),
          path("${vcf.simpleName}.sample_qc.tsv.gz"),
          path("${vcf.simpleName}.missing_site.tsv"),
          path("${vcf.simpleName}.missing_indv.tsv"),
          path("${vcf.simpleName}.freq.tsv"),
          path("${vcf.simpleName}.standardized_summary.txt"),
          path("${vcf.simpleName}_summary_plots.png"),
          path("${vcf.simpleName}_summary_plots_extra.png"),
          path("${vcf.simpleName}_pca.png"),
          path("${vcf.simpleName}_sample_qc_derived.tsv"),
          path("${vcf.simpleName}_locus_qc_derived.tsv"),
          path("${vcf.simpleName}_worst_samples.tsv"),
          path("${vcf.simpleName}_worst_loci.tsv"),
          emit: report_inputs

    script:
    def has_ploidy_map = ploidy_map.name != 'NO_FILE'
    def ploidy_arg     = has_ploidy_map ? "${ploidy_map}" : "NO_PLOIDY_MAP"
    """
    set -euo pipefail

    echo "=== SUMMARIZE_VCFS ===" > ${vcf.simpleName}.standardized_summary.txt
    echo "VCF: ${vcf}"           >> ${vcf.simpleName}.standardized_summary.txt
    echo "Caller: ${caller}"     >> ${vcf.simpleName}.standardized_summary.txt
    echo                         >> ${vcf.simpleName}.standardized_summary.txt

    # -----------------------------------------------------------------------
    # 1. Build the analysis-ready VCF (fill tags; strip/recompute for ANGSD)
    # -----------------------------------------------------------------------
    if [ "${caller}" = "angsd" ]; then
        echo "Preparing ANGSD-aware summary VCF..." | tee -a ${vcf.simpleName}.standardized_summary.txt
        bcftools annotate \
            -x INFO/AF,INFO/AC,INFO/AN,INFO/NS,INFO/MAF,INFO/F_MISSING \
            -Oz -o ${vcf.simpleName}.summary_input.vcf.gz "${vcf}"
        bcftools index ${vcf.simpleName}.summary_input.vcf.gz
        bcftools +fill-tags ${vcf.simpleName}.summary_input.vcf.gz \
            -Oz -o ${vcf.simpleName}.summary_ready.vcf.gz \
            -- -t AC,AN,AF,NS,MAF,TYPE,F_MISSING
        rm -f ${vcf.simpleName}.summary_input.vcf.gz \
              ${vcf.simpleName}.summary_input.vcf.gz.csi
    else
        echo "Preparing FreeBayes summary VCF..." | tee -a ${vcf.simpleName}.standardized_summary.txt
        bcftools +fill-tags "${vcf}" \
            -Oz -o ${vcf.simpleName}.summary_ready.vcf.gz \
            -- -t AC,AN,AF,NS,MAF,TYPE,F_MISSING
    fi
    bcftools index ${vcf.simpleName}.summary_ready.vcf.gz
    SUMMARY_VCF=${vcf.simpleName}.summary_ready.vcf.gz

    # -----------------------------------------------------------------------
    # 2. Lightweight per-site QC table (bcftools query — no R needed)
    # -----------------------------------------------------------------------
    {
        printf "chromo\tposition\tREF\tALT\tQUAL\tFILTER\tNS\tDP\tAF\n"
        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/NS\t%INFO/DP\t%INFO/AF\n' \
            \${SUMMARY_VCF}
    } | bgzip -c > ${vcf.simpleName}.site_qc.tsv.gz

    # -----------------------------------------------------------------------
    # 3. Per-sample QC table (Python streaming parser)
    # -----------------------------------------------------------------------
    python3 <<'PY'
import gzip, re

vcf     = "${vcf.simpleName}.summary_ready.vcf.gz"
samples = []

with gzip.open(vcf, "rt") as fh:
    for line in fh:
        if line.startswith("#CHROM"):
            samples = line.rstrip().split("\t")[9:]
            break

stats = {s: {"sites_total": 0, "sites_called": 0, "sites_missing": 0,
             "dp_sum": 0, "dp_n": 0, "het": 0, "hom_ref": 0, "hom_alt": 0}
         for s in samples}

with gzip.open(vcf, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        fields        = line.rstrip().split("\t")
        fmt           = fields[8].split(":")
        sample_fields = fields[9:]
        gt_i = fmt.index("GT") if "GT" in fmt else None
        dp_i = fmt.index("DP") if "DP" in fmt else None
        for s, val in zip(samples, sample_fields):
            st    = stats[s]
            parts = val.split(":")
            gt    = parts[gt_i] if gt_i is not None and gt_i < len(parts) else "./."
            dp    = parts[dp_i] if dp_i is not None and dp_i < len(parts) else "."
            st["sites_total"] += 1
            alleles = re.split(r'[/|]', gt)
            if all(a == '.' for a in alleles):
                st["sites_missing"] += 1
            else:
                st["sites_called"] += 1
                unique_a = set(alleles) - {'.'}
                if len(unique_a) > 1:
                    st["het"] += 1
                elif '0' in unique_a:
                    st["hom_ref"] += 1
                else:
                    st["hom_alt"] += 1
            if dp not in (".", ""):
                try:
                    st["dp_sum"] += int(dp); st["dp_n"] += 1
                except ValueError:
                    pass

with open("${vcf.simpleName}.sample_qc.tsv", "w") as out:
    out.write("sample\tsites_total\tsites_called\tsites_missing\t"
              "f_missing\tmean_dp_called\thet\thom_ref\thom_alt\n")
    for s in samples:
        st      = stats[s]
        f_miss  = st["sites_missing"] / st["sites_total"] if st["sites_total"] else 0
        mean_dp = st["dp_sum"]        / st["dp_n"]        if st["dp_n"]        else 0
        out.write(f"{s}\t{st['sites_total']}\t{st['sites_called']}\t"
                  f"{st['sites_missing']}\t{f_miss:.6f}\t{mean_dp:.4f}\t"
                  f"{st['het']}\t{st['hom_ref']}\t{st['hom_alt']}\n")
PY
    bgzip -f ${vcf.simpleName}.sample_qc.tsv

    # -----------------------------------------------------------------------
    # 4. bcftools stats, vcftools missingness + freq
    # -----------------------------------------------------------------------
    bcftools stats \${SUMMARY_VCF} > ${vcf.simpleName}.stats.txt

    vcftools --gzvcf \${SUMMARY_VCF} --missing-site --out ${vcf.simpleName} \
        >/dev/null 2>&1 || echo -e "CHR\tPOS\tN_DATA\tN_GENOTYPE_FILTERED\tN_MISS\tF_MISS" \
        > ${vcf.simpleName}.lmiss
    vcftools --gzvcf \${SUMMARY_VCF} --missing-indv --out ${vcf.simpleName} \
        >/dev/null 2>&1 || echo -e "INDV\tN_DATA\tN_GENOTYPES_FILTERED\tN_MISS\tF_MISS" \
        > ${vcf.simpleName}.imiss
    mv ${vcf.simpleName}.lmiss ${vcf.simpleName}.missing_site.tsv
    mv ${vcf.simpleName}.imiss ${vcf.simpleName}.missing_indv.tsv

    vcftools --gzvcf \${SUMMARY_VCF} --freq --out ${vcf.simpleName} \
        >/dev/null 2>&1 || echo -e "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}" \
        > ${vcf.simpleName}.frq
    mv ${vcf.simpleName}.frq ${vcf.simpleName}.freq.tsv

    # -----------------------------------------------------------------------
    # 5. Pre-extract DP and GT matrices for R plot script
    # -----------------------------------------------------------------------
    SAMPLES=\$(bcftools query -l \${SUMMARY_VCF} | tr '\n' '\t' | sed 's/\t\$//')

    echo -e "CHROM\tPOS\t\${SAMPLES}" > ${vcf.simpleName}.dp_matrix.tsv
    bcftools query -f '%CHROM\t%POS[\t%DP]\n' \${SUMMARY_VCF} \
        >> ${vcf.simpleName}.dp_matrix.tsv
    bgzip -f ${vcf.simpleName}.dp_matrix.tsv

    # -----------------------------------------------------------------------
    # 6. Pre-select PCA sites in the shell — sample up to 10k sites with
    #    MAF > 0.05 from the site_qc table, then extract a tiny VCF for R.
    #    R never touches the full VCF.
    # -----------------------------------------------------------------------
    python3 <<'PY'
import gzip, random, sys

site_qc = "${vcf.simpleName}.site_qc.tsv.gz"
max_sites = 10000
maf_min   = 0.05
random.seed(42)

candidates = []
with gzip.open(site_qc, "rt") as fh:
    header = fh.readline()   # skip header
    for line in fh:
        parts = line.rstrip().split("\t")
        if len(parts) < 9:
            continue
        chrom, pos, af_str = parts[0], parts[1], parts[8]
        try:
            af  = float(af_str.split(",")[0])
            maf = min(af, 1.0 - af)
        except (ValueError, IndexError):
            continue
        if maf >= maf_min:
            candidates.append((chrom, pos))

if len(candidates) > max_sites:
    candidates = random.sample(candidates, max_sites)

candidates.sort()   # sort for bcftools regions (required)
with open("pca_regions.txt", "w") as out:
    for chrom, pos in candidates:
        out.write(f"{chrom}\t{pos}\n")

print(f"Selected {len(candidates)} sites for PCA subset", file=sys.stderr)
PY

    # Extract the small PCA VCF — only GT field needed, drop INFO to shrink size
    if [ -s pca_regions.txt ]; then
        bcftools view \${SUMMARY_VCF} \
            -R pca_regions.txt \
            -Oz -o ${vcf.simpleName}.pca_subset.vcf.gz
        bcftools index ${vcf.simpleName}.pca_subset.vcf.gz
        PCA_VCF=${vcf.simpleName}.pca_subset.vcf.gz
    else
        # Fallback: no sites passed MAF filter — pass empty sentinel
        PCA_VCF="NO_PCA_VCF"
    fi

    # -----------------------------------------------------------------------
    # 7. Standardized summary text
    # -----------------------------------------------------------------------
    {
        echo "=== STANDARDIZED RAW SNP SUMMARY ==="
        echo "Caller: ${caller}"
        echo "VCF summarized: \${SUMMARY_VCF}"
        echo
        echo "[Variant counts]"
        grep '^SN' ${vcf.simpleName}.stats.txt || true
        echo
        echo "[Ts/Tv]"
        grep '^TSTV' ${vcf.simpleName}.stats.txt || true
        echo
        echo "[Top of site QC table]"
        zcat ${vcf.simpleName}.site_qc.tsv.gz | head -5 || true
        echo
        echo "[Top of sample QC table]"
        zcat ${vcf.simpleName}.sample_qc.tsv.gz | head -5 || true
        echo
        echo "[Worst 10 samples by missingness]"
        tail -n +2 ${vcf.simpleName}.missing_indv.tsv | sort -k5,5gr | head -10 || true
        echo
        echo "[Top of site missingness table]"
        head -10 ${vcf.simpleName}.missing_site.tsv || true
    } >> ${vcf.simpleName}.standardized_summary.txt

    VARIANT_COUNT=\$(bcftools view -H \${SUMMARY_VCF} | wc -l)
    echo "Total variants: \${VARIANT_COUNT}" | tee -a ${vcf.simpleName}.standardized_summary.txt

    # -----------------------------------------------------------------------
    # 8. R QC plots — receives only pre-extracted flat files + small PCA VCF
    # -----------------------------------------------------------------------
    if [ \${VARIANT_COUNT} -gt 10 ]; then
        echo "Running QC plot generation..." | tee -a ${vcf.simpleName}.standardized_summary.txt
        Rscript ${projectDir}/r_scripts/vcf_qc_plots.R \
            \${PCA_VCF} \
            ${vcf.simpleName} \
            ${ploidy_arg} \
            ${vcf.simpleName}.dp_matrix.tsv.gz
    else
        echo "Too few variants (\${VARIANT_COUNT}) for QC plots." \
            | tee -a ${vcf.simpleName}.standardized_summary.txt
        Rscript -e "
        library(ggplot2)
        p <- ggplot() +
          annotate('text', x=0.5, y=0.5,
                   label='Too few variants (n=\${VARIANT_COUNT})\\nfor QC analysis', size=8) +
          theme_void() + theme(panel.border=element_rect(fill=NA))
        for (f in c('${vcf.simpleName}_summary_plots.png',
                    '${vcf.simpleName}_summary_plots_extra.png',
                    '${vcf.simpleName}_pca.png'))
          ggsave(f, p, width=10, height=10)
        empty <- data.frame(message='Too few variants')
        for (f in c('${vcf.simpleName}_sample_qc_derived.tsv',
                    '${vcf.simpleName}_locus_qc_derived.tsv',
                    '${vcf.simpleName}_worst_samples.tsv',
                    '${vcf.simpleName}_worst_loci.tsv'))
          write.table(empty, f, sep='\t', quote=FALSE, row.names=FALSE)
        "
    fi

    # Clean up intermediates not in publishDir outputs
    rm -f ${vcf.simpleName}.summary_ready.vcf.gz \
          ${vcf.simpleName}.summary_ready.vcf.gz.csi \
          ${vcf.simpleName}.dp_matrix.tsv.gz \
          ${vcf.simpleName}.pca_subset.vcf.gz \
          ${vcf.simpleName}.pca_subset.vcf.gz.csi \
          pca_regions.txt

    echo "Summarization complete for ${vcf}" | tee -a ${vcf.simpleName}.standardized_summary.txt
    """
}