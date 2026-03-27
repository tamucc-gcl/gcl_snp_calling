process SUMMARIZE_VCFS {
    tag "${vcf.baseName}"
    publishDir params.output_dir, mode: 'copy'

    input:
    path vcf
    path ploidy_map
    val caller

    output:
    path("${vcf.simpleName}.stats.txt")
    path("${vcf.simpleName}.site_qc.tsv.gz")
    path("${vcf.simpleName}.sample_qc.tsv.gz")
    path("${vcf.simpleName}.missing_site.tsv")
    path("${vcf.simpleName}.missing_indv.tsv")
    path("${vcf.simpleName}.freq.tsv")
    path("${vcf.simpleName}.standardized_summary.txt")
    path("${vcf.simpleName}_summary_plots.png")
    path("${vcf.simpleName}_summary_plots_extra.png")
    path("${vcf.simpleName}_pca.png")
    path("${vcf.simpleName}_sample_qc_derived.tsv")
    path("${vcf.simpleName}_locus_qc_derived.tsv")
    path("${vcf.simpleName}_worst_samples.tsv")
    path("${vcf.simpleName}_worst_loci.tsv")

    script:
    def has_ploidy_map = ploidy_map.name != 'NO_FILE'
    def ploidy_arg = has_ploidy_map ? "${ploidy_map}" : "NO_PLOIDY_MAP"
    """
    set -euo pipefail

    CALLER="${caller}"
    VCF="${vcf}"
    PREFIX="${vcf.simpleName}"

    echo "=== SUMMARIZE_VCFS ===" > \${PREFIX}.standardized_summary.txt
    echo "VCF: ${vcf}" >> \${PREFIX}.standardized_summary.txt
    echo "Caller: ${caller}" >> \${PREFIX}.standardized_summary.txt
    echo >> \${PREFIX}.standardized_summary.txt

    # ------------------------------------------------------------------
    # 1) Header-aware normalization for summaries
    # ------------------------------------------------------------------
    # FreeBayes VCFs are usually already bcftools-friendly.
    # ANGSD VCFs often need INFO tags recalculated from GT and their original
    # INFO/AF should be preserved separately because it is GL-based.
    # We therefore create a lightweight summary VCF with standardized tags.

    if [ "${caller}" = "angsd" ]; then
        echo "Preparing ANGSD-aware summary VCF..." | tee -a \${PREFIX}.standardized_summary.txt

        bcftools annotate \
            -x INFO/AF,INFO/AC,INFO/AN,INFO/NS,INFO/MAF,INFO/F_MISSING \
            -Oz -o \${PREFIX}.summary_input.vcf.gz \
            "${vcf}"
        tabix -p vcf \${PREFIX}.summary_input.vcf.gz

        bcftools +fill-tags \${PREFIX}.summary_input.vcf.gz \
            -Oz -o \${PREFIX}.summary_ready.vcf.gz \
            -- -t AC,AN,AF,NS,MAF,TYPE,F_MISSING
        tabix -p vcf \${PREFIX}.summary_ready.vcf.gz

        SUMMARY_VCF=\${PREFIX}.summary_ready.vcf.gz

        # Preserve ANGSD site-level INFO fields directly from the raw VCF.
        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/NS\t%INFO/DP\t%INFO/AF\n' \
            "${vcf}" | bgzip -c > \${PREFIX}.site_qc.tsv.gz

        # Build sample table from GT/DP missingness in the summary-ready VCF.
        python3 << PY
import gzip

vcf = "${PREFIX}.summary_ready.vcf.gz"
samples = []
with gzip.open(vcf, 'rt') as fh:
    for line in fh:
        if line.startswith('#CHROM'):
            samples = line.rstrip().split('\t')[9:]
            break

stats = {s: {'sites_total': 0, 'sites_called': 0, 'sites_missing': 0,
             'dp_sum': 0, 'dp_n': 0, 'het': 0, 'hom_ref': 0, 'hom_alt': 0}
         for s in samples}

with gzip.open(vcf, 'rt') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        fields = line.rstrip().split('\t')
        fmt = fields[8].split(':')
        sample_fields = fields[9:]
        try:
            gt_i = fmt.index('GT')
        except ValueError:
            gt_i = None
        try:
            dp_i = fmt.index('DP')
        except ValueError:
            dp_i = None

        for s, val in zip(samples, sample_fields):
            stats[s]['sites_total'] += 1
            parts = val.split(':')
            gt = parts[gt_i] if gt_i is not None and gt_i < len(parts) else './.'
            dp = parts[dp_i] if dp_i is not None and dp_i < len(parts) else '.'

            if gt in ('./.', '.|.'):
                stats[s]['sites_missing'] += 1
            else:
                stats[s]['sites_called'] += 1
                if gt in ('0/1', '1/0', '0|1', '1|0'):
                    stats[s]['het'] += 1
                elif gt in ('0/0', '0|0'):
                    stats[s]['hom_ref'] += 1
                elif gt in ('1/1', '1|1'):
                    stats[s]['hom_alt'] += 1

            if dp not in ('.', ''):
                try:
                    dpi = int(dp)
                    stats[s]['dp_sum'] += dpi
                    stats[s]['dp_n'] += 1
                except ValueError:
                    pass

with open("${PREFIX}.sample_qc.tsv", 'w') as out:
    out.write('sample\tsites_total\tsites_called\tsites_missing\tf_missing\tmean_dp_called\thet\thom_ref\thom_alt\n')
    for s in samples:
        st = stats[s]
        f_missing = st['sites_missing'] / st['sites_total'] if st['sites_total'] else 0
        mean_dp = st['dp_sum'] / st['dp_n'] if st['dp_n'] else 0
        out.write(f"{s}\t{st['sites_total']}\t{st['sites_called']}\t{st['sites_missing']}\t{f_missing:.6f}\t{mean_dp:.4f}\t{st['het']}\t{st['hom_ref']}\t{st['hom_alt']}\n")
PY
        bgzip -f \${PREFIX}.sample_qc.tsv

    else
        echo "Preparing FreeBayes summary VCF..." | tee -a \${PREFIX}.standardized_summary.txt

        bcftools +fill-tags "${vcf}" \
            -Oz -o \${PREFIX}.summary_ready.vcf.gz \
            -- -t AC,AN,AF,NS,MAF,TYPE,F_MISSING
        tabix -p vcf \${PREFIX}.summary_ready.vcf.gz

        SUMMARY_VCF=\${PREFIX}.summary_ready.vcf.gz

        bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/NS\t%INFO/DP\t%INFO/AF\n' \
            \${PREFIX}.summary_ready.vcf.gz | bgzip -c > \${PREFIX}.site_qc.tsv.gz

        python3 << PY
import gzip

vcf = "${PREFIX}.summary_ready.vcf.gz"
samples = []
with gzip.open(vcf, 'rt') as fh:
    for line in fh:
        if line.startswith('#CHROM'):
            samples = line.rstrip().split('\t')[9:]
            break

stats = {s: {'sites_total': 0, 'sites_called': 0, 'sites_missing': 0,
             'dp_sum': 0, 'dp_n': 0, 'het': 0, 'hom_ref': 0, 'hom_alt': 0}
         for s in samples}

with gzip.open(vcf, 'rt') as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        fields = line.rstrip().split('\t')
        fmt = fields[8].split(':')
        sample_fields = fields[9:]
        try:
            gt_i = fmt.index('GT')
        except ValueError:
            gt_i = None
        try:
            dp_i = fmt.index('DP')
        except ValueError:
            dp_i = None

        for s, val in zip(samples, sample_fields):
            stats[s]['sites_total'] += 1
            parts = val.split(':')
            gt = parts[gt_i] if gt_i is not None and gt_i < len(parts) else './.'
            dp = parts[dp_i] if dp_i is not None and dp_i < len(parts) else '.'

            if gt in ('./.', '.|.'):
                stats[s]['sites_missing'] += 1
            else:
                stats[s]['sites_called'] += 1
                if gt in ('0/1', '1/0', '0|1', '1|0'):
                    stats[s]['het'] += 1
                elif gt in ('0/0', '0|0'):
                    stats[s]['hom_ref'] += 1
                elif gt in ('1/1', '1|1'):
                    stats[s]['hom_alt'] += 1

            if dp not in ('.', ''):
                try:
                    dpi = int(dp)
                    stats[s]['dp_sum'] += dpi
                    stats[s]['dp_n'] += 1
                except ValueError:
                    pass

with open("${PREFIX}.sample_qc.tsv", 'w') as out:
    out.write('sample\tsites_total\tsites_called\tsites_missing\tf_missing\tmean_dp_called\thet\thom_ref\thom_alt\n')
    for s in samples:
        st = stats[s]
        f_missing = st['sites_missing'] / st['sites_total'] if st['sites_total'] else 0
        mean_dp = st['dp_sum'] / st['dp_n'] if st['dp_n'] else 0
        out.write(f"{s}\t{st['sites_total']}\t{st['sites_called']}\t{st['sites_missing']}\t{f_missing:.6f}\t{mean_dp:.4f}\t{st['het']}\t{st['hom_ref']}\t{st['hom_alt']}\n")
PY
        bgzip -f \${PREFIX}.sample_qc.tsv
    fi

    # ------------------------------------------------------------------
    # 2) Standardized summaries shared by both callers
    # ------------------------------------------------------------------
    bcftools stats \${SUMMARY_VCF} > \${PREFIX}.stats.txt

    vcftools --gzvcf \${SUMMARY_VCF} --missing-site --out \${PREFIX} >/dev/null 2>&1
    vcftools --gzvcf \${SUMMARY_VCF} --missing-indv --out \${PREFIX} >/dev/null 2>&1
    mv \${PREFIX}.lmiss \${PREFIX}.missing_site.tsv
    mv \${PREFIX}.imiss \${PREFIX}.missing_indv.tsv

    vcftools --gzvcf \${SUMMARY_VCF} --freq --out \${PREFIX} >/dev/null 2>&1
    mv \${PREFIX}.frq \${PREFIX}.freq.tsv

    {
        echo "=== STANDARDIZED RAW SNP SUMMARY ==="
        echo "Caller: ${caller}"
        echo "VCF summarized: \${SUMMARY_VCF}"
        echo
        echo "[Variant counts]"
        grep '^SN' \${PREFIX}.stats.txt || true
        echo
        echo "[Ts/Tv]"
        grep '^TSTV' \${PREFIX}.stats.txt || true
        echo
        echo "[Top of site QC table]"
        zcat \${PREFIX}.site_qc.tsv.gz | head -5 || true
        echo
        echo "[Top of sample QC table]"
        zcat \${PREFIX}.sample_qc.tsv.gz | head -5 || true
        echo
        echo "[Worst 10 samples by missingness]"
        tail -n +2 \${PREFIX}.missing_indv.tsv | sort -k5,5gr | head -10 || true
        echo
        echo "[Top of site missingness table]"
        head -10 \${PREFIX}.missing_site.tsv || true
    } >> \${PREFIX}.standardized_summary.txt

    # ------------------------------------------------------------------
    # 3) Plot workflow using standardized companion files
    # ------------------------------------------------------------------
    VARIANT_COUNT=\$(bcftools view -H \${SUMMARY_VCF} | wc -l)
    echo "Total variants in summary VCF: \${VARIANT_COUNT}" | tee -a \${PREFIX}.standardized_summary.txt
    echo "Has ploidy map: ${has_ploidy_map}" | tee -a \${PREFIX}.standardized_summary.txt
    echo "Ploidy argument to R: ${ploidy_arg}" | tee -a \${PREFIX}.standardized_summary.txt

    if [ \${VARIANT_COUNT} -gt 10 ]; then
        echo "Running QC plot generation..." | tee -a \${PREFIX}.standardized_summary.txt
        Rscript ${projectDir}/r_scripts/vcf_qc_plots.R \${SUMMARY_VCF} \${PREFIX} ${ploidy_arg}
    else
        echo "Too few variants (\${VARIANT_COUNT}) for meaningful QC plots. Creating placeholder images..." | tee -a \${PREFIX}.standardized_summary.txt
        Rscript -e "
        library(ggplot2)
        placeholder <- ggplot() +
            annotate('text', x = 0.5, y = 0.5,
                    label = 'Too few variants (n=\${VARIANT_COUNT})\\nfor QC analysis',
                    size = 8) +
            theme_void() +
            theme(panel.border = element_rect(fill=NA))
        ggsave('\${PREFIX}_summary_plots.png', placeholder, width=10, height=10)
        ggsave('\${PREFIX}_summary_plots_extra.png', placeholder, width=10, height=10)
        ggsave('\${PREFIX}_pca.png', placeholder, width=5, height=5)
        write.table(data.frame(message='Too few variants for derived sample QC'),
                    file='\${PREFIX}_sample_qc_derived.tsv', sep='\\t', quote=FALSE, row.names=FALSE)
        write.table(data.frame(message='Too few variants for derived locus QC'),
                    file='\${PREFIX}_locus_qc_derived.tsv', sep='\\t', quote=FALSE, row.names=FALSE)
        write.table(data.frame(message='Too few variants for worst sample ranking'),
                    file='\${PREFIX}_worst_samples.tsv', sep='\\t', quote=FALSE, row.names=FALSE)
        write.table(data.frame(message='Too few variants for worst locus ranking'),
                    file='\${PREFIX}_worst_loci.tsv', sep='\\t', quote=FALSE, row.names=FALSE)
        "
    fi

    echo "Summarization complete for ${vcf}" | tee -a \${PREFIX}.standardized_summary.txt
    """
}
