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
          path("${vcf.simpleName}_pca_site_selection.txt"),
          emit: report_inputs

    script:
    def has_ploidy_map = ploidy_map.name != 'NO_FILE'
    def ploidy_arg     = has_ploidy_map ? "${ploidy_map}" : "NO_PLOIDY_MAP"
    def pca_min_spacing = params.pca_min_spacing != null ? params.pca_min_spacing :
                        (params.data_type == 'wgs' ? 5000 : 500)
    def pca_ns_frac     = params.pca_ns_frac     != null ? params.pca_ns_frac     : 0.50

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
    cols = ["sample", "sites_total", "sites_called", "sites_missing",
            "f_missing", "mean_dp_called", "het", "hom_ref", "hom_alt"]
    print("\t".join(cols), file=out)
    for s in samples:
        st      = stats[s]
        f_miss  = st["sites_missing"] / st["sites_total"] if st["sites_total"] else 0
        mean_dp = st["dp_sum"]        / st["dp_n"]        if st["dp_n"]        else 0
        row = [
            s,
            str(st["sites_total"]),
            str(st["sites_called"]),
            str(st["sites_missing"]),
            f"{f_miss:.6f}",
            f"{mean_dp:.4f}",
            str(st["het"]),
            str(st["hom_ref"]),
            str(st["hom_alt"]),
        ]
        print("\t".join(row), file=out)
PY
    bgzip -f ${vcf.simpleName}.sample_qc.tsv

    # -----------------------------------------------------------------------
    # 4. bcftools stats, vcftools missingness + freq
    # -----------------------------------------------------------------------
    bcftools stats \${SUMMARY_VCF} > ${vcf.simpleName}.stats.txt

    if ! vcftools --gzvcf \${SUMMARY_VCF} --missing-site --out ${vcf.simpleName} >/dev/null 2>&1; then
        printf "CHR\tPOS\tN_DATA\tN_GENOTYPE_FILTERED\tN_MISS\tF_MISS\n" > ${vcf.simpleName}.lmiss
    fi
    mv ${vcf.simpleName}.lmiss ${vcf.simpleName}.missing_site.tsv

    if ! vcftools --gzvcf \${SUMMARY_VCF} --missing-indv --out ${vcf.simpleName} >/dev/null 2>&1; then
        printf "INDV\tN_DATA\tN_GENOTYPES_FILTERED\tN_MISS\tF_MISS\n" > ${vcf.simpleName}.imiss
    fi
    mv ${vcf.simpleName}.imiss ${vcf.simpleName}.missing_indv.tsv

    if ! vcftools --gzvcf \${SUMMARY_VCF} --freq --out ${vcf.simpleName} >/dev/null 2>&1; then
        printf "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}\n" > ${vcf.simpleName}.frq
    fi
    mv ${vcf.simpleName}.frq ${vcf.simpleName}.freq.tsv

    # -----------------------------------------------------------------------
    # 5. Pre-extract DP matrix for R plot script
    # -----------------------------------------------------------------------
    bcftools query -l \${SUMMARY_VCF} > vcf_samples.txt
    SAMPLES=\$(bcftools query -l \${SUMMARY_VCF} | tr '\n' '\t' | sed 's/\t\$//')

    echo -e "CHROM\tPOS\t\${SAMPLES}" > ${vcf.simpleName}.dp_matrix.tsv
    bcftools query -f '%CHROM\t%POS[\t%DP]\n' \${SUMMARY_VCF} \
        >> ${vcf.simpleName}.dp_matrix.tsv
    bgzip -f ${vcf.simpleName}.dp_matrix.tsv

    # -----------------------------------------------------------------------
    # 6. Pre-select PCA sites
    #
    # Requires BOTH a call-rate floor (INFO/NS) and a MAF floor. MAF alone is
    # not a meaningful threshold because AN depends on how many samples were
    # called at the site.
    # -----------------------------------------------------------------------
    python3 <<'PY'
import gzip, sys

site_qc     = "${vcf.simpleName}.site_qc.tsv.gz"
max_sites   = 10000
maf_min     = 0.05
min_spacing = ${pca_min_spacing}
ns_frac     = ${pca_ns_frac}

# Relax the call-rate floor rather than emit a PCA built on nothing.
ns_frac_ladder = [f for f in (ns_frac, 0.25, 0.10, 0.0) if f <= ns_frac] or [0.0]
min_sites_ok   = 1000

# site_qc.tsv.gz columns:
#   0 chromo  1 position  2 REF  3 ALT  4 QUAL  5 FILTER  6 NS  7 DP  8 AF
COL_CHROM, COL_POS, COL_NS, COL_AF = 0, 1, 6, 8

with open("vcf_samples.txt") as fh:
    n_samples = sum(1 for line in fh if line.strip())

if n_samples == 0:
    sys.exit("PCA site selection: no samples found in vcf_samples.txt")

log = []
log.append(f"Samples in callset: {n_samples}")
log.append(f"MAF floor: {maf_min}")
log.append(f"Minimum site spacing: {min_spacing} bp")


def parse(parts):
    """Return (chrom, pos, ns, maf) or None if the row is unusable."""
    if len(parts) <= COL_AF:
        return None
    try:
        ns = int(float(parts[COL_NS]))
    except (ValueError, IndexError):
        return None
    try:
        af = float(parts[COL_AF].split(",")[0])
    except (ValueError, IndexError):
        return None
    return parts[COL_CHROM], parts[COL_POS], ns, min(af, 1.0 - af)


# --- pass 1: how many sites survive at each rung of the ladder? ------------
counts = {f: 0 for f in ns_frac_ladder}
n_rows = 0
n_maf_pass = 0

with gzip.open(site_qc, "rt") as fh:
    fh.readline()
    for line in fh:
        n_rows += 1
        rec = parse(line.rstrip("\n").split("\t"))
        if rec is None:
            continue
        _, _, ns, maf = rec
        if maf < maf_min:
            continue
        n_maf_pass += 1
        for f in ns_frac_ladder:
            if ns >= f * n_samples:
                counts[f] += 1

log.append(f"Sites read: {n_rows}")
log.append(f"Sites passing MAF alone: {n_maf_pass}")
for f in ns_frac_ladder:
    log.append(f"  + NS >= {f:.2f} x {n_samples} = {int(f * n_samples)}: {counts[f]}")

chosen = next((f for f in ns_frac_ladder if counts[f] >= min_sites_ok),
              ns_frac_ladder[-1])
ns_min = int(chosen * n_samples)

if chosen != ns_frac_ladder[0]:
    log.append(f"WARNING: call-rate floor relaxed from {ns_frac_ladder[0]:.2f} "
               f"to {chosen:.2f} to reach {min_sites_ok} sites. "
               f"Interpret the PCA with care.")
if counts[chosen] < min_sites_ok:
    log.append(f"WARNING: only {counts[chosen]} sites available even with no "
               f"call-rate floor. The PCA is unlikely to be informative.")

log.append(f"Selected NS floor: {ns_min} ({chosen:.2f} of {n_samples})")


# --- pass 2: collect, enforcing spacing on the fly ------------------------
# site_qc.tsv.gz is written in VCF order, so it is already coordinate-sorted
# and a single greedy sweep is sufficient.
candidates = []
last_chrom = None
last_pos = None
n_dropped_spacing = 0

with gzip.open(site_qc, "rt") as fh:
    fh.readline()
    for line in fh:
        rec = parse(line.rstrip("\n").split("\t"))
        if rec is None:
            continue
        chrom, pos, ns, maf = rec
        if maf < maf_min or ns < ns_min:
            continue

        if min_spacing > 0:
            try:
                ipos = int(pos)
            except ValueError:
                continue
            if chrom == last_chrom and last_pos is not None \
                    and (ipos - last_pos) < min_spacing:
                n_dropped_spacing += 1
                continue
            last_chrom, last_pos = chrom, ipos

        candidates.append((chrom, pos))

log.append(f"Sites dropped by spacing: {n_dropped_spacing}")
log.append(f"Sites after spacing: {len(candidates)}")


# --- thin to max_sites, preserving even spread across the genome ----------
# Systematic every-nth thinning rather than random sampling, so the subset
# stays spread across contigs instead of clumping by chance.
if len(candidates) > max_sites:
    step = len(candidates) / float(max_sites)
    candidates = [candidates[int(i * step)] for i in range(max_sites)]
    log.append(f"Thinned to {len(candidates)} sites (every ~{step:.1f}th)")

with open("pca_regions.txt", "w") as out:
    for chrom, pos in candidates:
        print(chrom, pos, sep="\t", file=out)

with open("pca_site_selection.txt", "w") as out:
    print("=== PCA SITE SELECTION ===", file=out)
    for line in log:
        print(line, file=out)

for line in log:
    print(line, file=sys.stderr)
print(f"Selected {len(candidates)} sites for PCA subset", file=sys.stderr)
PY

    cat pca_site_selection.txt >> ${vcf.simpleName}.standardized_summary.txt
    echo >> ${vcf.simpleName}.standardized_summary.txt


    # Extract the small PCA VCF
    if [ -s pca_regions.txt ]; then
        bcftools view \${SUMMARY_VCF} \
            -R pca_regions.txt \
            -Oz -o ${vcf.simpleName}.pca_subset.vcf.gz
        bcftools index ${vcf.simpleName}.pca_subset.vcf.gz
        PCA_VCF=${vcf.simpleName}.pca_subset.vcf.gz
    else
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
    # 8. R QC plots
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

    # Clean up intermediates
    rm -f ${vcf.simpleName}.summary_ready.vcf.gz \
          ${vcf.simpleName}.summary_ready.vcf.gz.csi \
          ${vcf.simpleName}.dp_matrix.tsv.gz \
          ${vcf.simpleName}.pca_subset.vcf.gz \
          ${vcf.simpleName}.pca_subset.vcf.gz.csi

    mv pca_site_selection.txt ${vcf.simpleName}_pca_site_selection.txt
    echo "Summarization complete for ${vcf}" | tee -a ${vcf.simpleName}.standardized_summary.txt
    """
}