process REPORT_SNP_CALLING_SUMMARY {
    tag "${prefix}"
    publishDir params.output_dir, mode: 'copy'

    input:
    val prefix
    val caller
    path raw_vcf
    path stats_txt
    path standardized_summary
    path summary_plot
    path extra_plot
    path pca_plot
    path worst_samples
    path worst_loci
    path sample_qc_derived
    path locus_qc_derived
    path missing_indv
    path missing_site
    path freq_tsv
    path site_qc_tsv
    path sample_qc_tsv

    output:
    path("${prefix}_snp_calling_report.md")

    script:
    """
    set -euo pipefail

    python3 - \
        "${prefix}" \
        "${caller}" \
        "${raw_vcf}" \
        "${stats_txt}" \
        "${standardized_summary}" \
        "${summary_plot}" \
        "${extra_plot}" \
        "${pca_plot}" \
        "${worst_samples}" \
        "${worst_loci}" \
        "${sample_qc_derived}" \
        "${locus_qc_derived}" \
        "${missing_indv}" \
        "${missing_site}" \
        "${freq_tsv}" \
        "${site_qc_tsv}" \
        "${sample_qc_tsv}" <<'PY'
import csv
import gzip
import math
import sys
from pathlib import Path

(
    prefix,
    caller,
    raw_vcf,
    stats_txt,
    standardized_summary,
    summary_plot,
    extra_plot,
    pca_plot,
    worst_samples,
    worst_loci,
    sample_qc_derived,
    locus_qc_derived,
    missing_indv,
    missing_site,
    freq_tsv,
    site_qc_tsv,
    sample_qc_tsv,
) = sys.argv[1:18]

raw_vcf = Path(raw_vcf)
stats_txt = Path(stats_txt)
standardized_summary = Path(standardized_summary)
summary_plot = Path(summary_plot)
extra_plot = Path(extra_plot)
pca_plot = Path(pca_plot)
worst_samples = Path(worst_samples)
worst_loci = Path(worst_loci)
sample_qc_derived = Path(sample_qc_derived)
locus_qc_derived = Path(locus_qc_derived)
missing_indv = Path(missing_indv)
missing_site = Path(missing_site)
freq_tsv = Path(freq_tsv)
site_qc_tsv = Path(site_qc_tsv)
sample_qc_tsv = Path(sample_qc_tsv)

out_md = Path(f"{prefix}_snp_calling_report.md")

# summarize_vcfs outputs are now published under params.output_dir/snp_qc
QC_REL_DIR = Path("snp_qc")


def rel_link_qc(path: Path) -> str:
    return str(QC_REL_DIR / path.name)


def rel_link_root(path: Path) -> str:
    return path.name


def pct(x, digits=1):
    if x is None:
        return "NA"
    try:
        return f"{100 * float(x):,.{digits}f}%"
    except Exception:
        return "NA"


def num(x, digits=2):
    if x is None:
        return "NA"
    try:
        xf = float(x)
        if not math.isfinite(xf):
            return "NA"
        if digits == 0:
            return f"{int(round(xf)):,}"
        return f"{xf:,.{digits}f}"
    except Exception:
        return str(x)


def open_text_auto(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "r")


def read_tsv(path: Path):
    with open_text_auto(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader)


def safe_read_tsv(path: Path):
    try:
        return read_tsv(path)
    except Exception:
        return []


def md_table(rows, columns, max_rows=10):
    rows = rows[:max_rows]
    if not rows or not columns:
        return "_No rows available._\n"
    header = "| " + " | ".join(columns) + " |\n"
    sep = "| " + " | ".join(["---"] * len(columns)) + " |\n"
    body = ""
    for r in rows:
        body += "| " + " | ".join(str(r.get(c, "")) for c in columns) + " |\n"
    return header + sep + body


# ------------------------------------------------------------
# Parse bcftools stats summary
# ------------------------------------------------------------
summary_numbers = {}
ts_tv = None

with open(stats_txt, "r") as fh:
    for line in fh:
        if line.startswith("SN\t"):
            parts = line.rstrip().split("\t")
            if len(parts) >= 4:
                summary_numbers[parts[2].rstrip(":" )] = parts[3]
        elif line.startswith("TSTV\t"):
            parts = line.rstrip().split("\t")
            if len(parts) >= 5:
                ts_tv = {
                    "ts": parts[2],
                    "tv": parts[3],
                    "ratio": parts[4],
                }

# ------------------------------------------------------------
# Load derived tables
# ------------------------------------------------------------
sample_rows = safe_read_tsv(sample_qc_derived)
locus_rows = safe_read_tsv(locus_qc_derived)
worst_sample_rows = safe_read_tsv(worst_samples)
worst_locus_rows = safe_read_tsv(worst_loci)

# ------------------------------------------------------------
# Compute report metrics
# ------------------------------------------------------------
metrics = {}

if sample_rows and "sample" in sample_rows[0]:
    miss = [float(r["pct_loci_missing"]) for r in sample_rows if r.get("pct_loci_missing") not in (None, "", "NA")]
    mean_dp = [float(r["mean_depth_called"]) for r in sample_rows if r.get("mean_depth_called") not in (None, "", "NA")]
    n_loci = [float(r["number_loci"]) for r in sample_rows if r.get("number_loci") not in (None, "", "NA")]

    metrics["samples_total"] = len(sample_rows)
    metrics["sample_missing_median"] = sorted(miss)[len(miss) // 2] if miss else None
    metrics["sample_missing_gt10"] = sum(x > 0.10 for x in miss) if miss else None
    metrics["sample_missing_gt20"] = sum(x > 0.20 for x in miss) if miss else None
    metrics["sample_missing_gt30"] = sum(x > 0.30 for x in miss) if miss else None
    metrics["sample_mean_dp_median"] = sorted(mean_dp)[len(mean_dp) // 2] if mean_dp else None
    metrics["sample_called_loci_median"] = sorted(n_loci)[len(n_loci) // 2] if n_loci else None

if locus_rows and "locus" in locus_rows[0]:
    locus_missing = [float(r["pct_samples_missing"]) for r in locus_rows if r.get("pct_samples_missing") not in (None, "", "NA")]
    shared = [float(r["number_samples_with_locus"]) for r in locus_rows if r.get("number_samples_with_locus") not in (None, "", "NA")]
    locus_depth = [float(r["mean_depth"]) for r in locus_rows if r.get("mean_depth") not in (None, "", "NA")]
    maf = [float(r["maf"]) for r in locus_rows if r.get("maf") not in (None, "", "NA")]

    metrics["loci_total"] = len(locus_rows)
    metrics["locus_missing_median"] = sorted(locus_missing)[len(locus_missing) // 2] if locus_missing else None
    metrics["locus_missing_gt10"] = sum(x > 0.10 for x in locus_missing) if locus_missing else None
    metrics["locus_missing_gt20"] = sum(x > 0.20 for x in locus_missing) if locus_missing else None
    metrics["locus_missing_gt30"] = sum(x > 0.30 for x in locus_missing) if locus_missing else None
    metrics["samples_per_locus_median"] = sorted(shared)[len(shared) // 2] if shared else None
    metrics["locus_depth_median"] = sorted(locus_depth)[len(locus_depth) // 2] if locus_depth else None
    metrics["maf_median"] = sorted(maf)[len(maf) // 2] if maf else None
    metrics["maf_lt01"] = sum(x < 0.01 for x in maf) if maf else None
    metrics["maf_lt05"] = sum(x < 0.05 for x in maf) if maf else None
    metrics["maf_lt10"] = sum(x < 0.10 for x in maf) if maf else None

# ------------------------------------------------------------
# Interpretation text
# ------------------------------------------------------------
interpretation = []

n_records = summary_numbers.get("number of records")
if n_records is not None:
    interpretation.append(
        f"The raw callset contains **{num(n_records, 0)} variant records** across **{num(summary_numbers.get('number of samples', metrics.get('samples_total', 'NA')), 0)} samples**."
    )

if ts_tv and ts_tv.get("ratio") not in (None, "", "NA"):
    interpretation.append(
        f"The **transition/transversion ratio is {num(ts_tv['ratio'], 2)}**, which is a useful high-level check on callset plausibility."
    )

if metrics.get("sample_missing_gt30") is not None:
    interpretation.append(
        f"At the sample level, **{num(metrics['sample_missing_gt30'], 0)} samples exceed 30% missing data**, with median sample missingness of **{pct(metrics.get('sample_missing_median'))}**."
    )

if metrics.get("locus_missing_gt30") is not None:
    interpretation.append(
        f"At the locus level, **{num(metrics['locus_missing_gt30'], 0)} loci exceed 30% missing data**, with median locus missingness of **{pct(metrics.get('locus_missing_median'))}**."
    )

if metrics.get("samples_per_locus_median") is not None:
    interpretation.append(
        f"The median locus is present in **{num(metrics['samples_per_locus_median'], 0)} samples**, which helps indicate whether the calling settings are retaining broadly shared loci or sparse sites."
    )

if metrics.get("maf_median") is not None:
    interpretation.append(
        f"The median minor allele frequency is **{num(metrics['maf_median'], 3)}**; loci below 0.05 and 0.10 are often worth tracking when tuning raw-call parameters for downstream GWAS filtering."
    )

caller_specific = []
if caller.lower() == "angsd":
    caller_specific.append(
        "Because this callset originates from **ANGSD**, the likelihood-based site summaries and standardized companion tables are more trustworthy for tuning than raw hard-call summaries alone."
    )
else:
    caller_specific.append(
        "Because this callset originates from **FreeBayes**, the genotype-based summaries are usually more directly interpretable as conventional VCF metrics."
    )

# ------------------------------------------------------------
# Build markdown
# ------------------------------------------------------------
with open(out_md, "w") as out:
    out.write(f"# SNP calling summary report: `{prefix}`\n\n")
    out.write(f"**Caller:** `{caller}`  \n")
    out.write(f"**Raw VCF:** [`{raw_vcf.name}`]({rel_link_root(raw_vcf)})  \n")
    out.write(f"**Primary VCF summary file:** [`{stats_txt.name}`]({rel_link_qc(stats_txt)})  \n")
    out.write(f"**Standardized summary text:** [`{standardized_summary.name}`]({rel_link_qc(standardized_summary)})\n\n")

    out.write("## Main takeaways\n\n")
    for line in interpretation + caller_specific:
        out.write(f"- {line}\n")
    out.write("\n")

    out.write("## Core summary metrics\n\n")
    out.write("| Metric | Value |\n")
    out.write("| --- | --- |\n")
    out.write(f"| Samples | {num(summary_numbers.get('number of samples', metrics.get('samples_total', 'NA')), 0)} |\n")
    out.write(f"| Variant records | {num(summary_numbers.get('number of records', metrics.get('loci_total', 'NA')), 0)} |\n")
    out.write(f"| SNPs | {num(summary_numbers.get('number of SNPs', 'NA'), 0)} |\n")
    out.write(f"| Multiallelic sites | {num(summary_numbers.get('number of multiallelic sites', 'NA'), 0)} |\n")
    out.write(f"| Ts/Tv | {num(ts_tv['ratio'], 2) if ts_tv else 'NA'} |\n")
    out.write(f"| Median sample missingness | {pct(metrics.get('sample_missing_median'))} |\n")
    out.write(f"| Samples >30% missing | {num(metrics.get('sample_missing_gt30'), 0)} |\n")
    out.write(f"| Median locus missingness | {pct(metrics.get('locus_missing_median'))} |\n")
    out.write(f"| Loci >30% missing | {num(metrics.get('locus_missing_gt30'), 0)} |\n")
    out.write(f"| Median samples per locus | {num(metrics.get('samples_per_locus_median'), 0)} |\n")
    out.write(f"| Median sample mean depth | {num(metrics.get('sample_mean_dp_median'), 2)} |\n")
    out.write(f"| Median locus mean depth | {num(metrics.get('locus_depth_median'), 2)} |\n")
    out.write(f"| Median MAF | {num(metrics.get('maf_median'), 3)} |\n")
    out.write(f"| Loci with MAF <0.05 | {num(metrics.get('maf_lt05'), 0)} |\n")
    out.write(f"| Loci with MAF <0.10 | {num(metrics.get('maf_lt10'), 0)} |\n\n")

    out.write("## QC plots\n\n")
    out.write("### Summary panel\n\n")
    out.write(f"![Summary plots]({rel_link_qc(summary_plot)})\n\n")
    out.write("### Additional diagnostics\n\n")
    out.write(f"![Additional summary plots]({rel_link_qc(extra_plot)})\n\n")
    out.write("### PCA\n\n")
    out.write(f"![PCA plot]({rel_link_qc(pca_plot)})\n\n")

    out.write("## Candidate samples for review or removal\n\n")
    out.write("These are the worst samples by derived sample QC, intended to guide inspection of missingness and depth before downstream filtering.\n\n")
    sample_cols = [c for c in ["sample", "pct_loci_missing", "number_loci", "mean_depth_called", "flag_high_missing", "flag_low_depth"] if worst_sample_rows and c in worst_sample_rows[0]]
    out.write(md_table(worst_sample_rows, sample_cols, max_rows=15))
    out.write("\n")

    out.write("## Candidate loci / locus classes for review\n\n")
    out.write("These are the worst loci by derived locus QC, useful for understanding whether raw calls are failing because of sparse representation, unusual depth, or both.\n\n")
    locus_cols = [c for c in ["locus", "pct_samples_missing", "number_samples_with_locus", "mean_depth", "maf"] if worst_locus_rows and c in worst_locus_rows[0]]
    out.write(md_table(worst_locus_rows, locus_cols, max_rows=15))
    out.write("\n")

    out.write("## Files to inspect next\n\n")
    out.write("| File | Why it matters |\n")
    out.write("| --- | --- |\n")
    out.write(f"| [`{sample_qc_derived.name}`]({rel_link_qc(sample_qc_derived)}) | Derived per-sample QC summary used for ranking poor samples. |\n")
    out.write(f"| [`{locus_qc_derived.name}`]({rel_link_qc(locus_qc_derived)}) | Derived per-locus QC summary used for missingness/depth/MAF interpretation. |\n")
    out.write(f"| [`{missing_indv.name}`]({rel_link_qc(missing_indv)}) | Raw per-sample missingness from vcftools. |\n")
    out.write(f"| [`{missing_site.name}`]({rel_link_qc(missing_site)}) | Raw per-site missingness from vcftools. |\n")
    out.write(f"| [`{freq_tsv.name}`]({rel_link_qc(freq_tsv)}) | Frequency table for checking the MAF distribution of raw SNP calls. |\n")
    out.write(f"| [`{site_qc_tsv.name}`]({rel_link_qc(site_qc_tsv)}) | Standardized site-level QC table; especially important for ANGSD outputs. |\n")
    out.write(f"| [`{sample_qc_tsv.name}`]({rel_link_qc(sample_qc_tsv)}) | Standardized sample-level QC table. |\n")
    out.write(f"| [`{worst_samples.name}`]({rel_link_qc(worst_samples)}) | Quick list of the worst samples to investigate first. |\n")
    out.write(f"| [`{worst_loci.name}`]({rel_link_qc(worst_loci)}) | Quick list of the worst loci to investigate first. |\n\n")

    out.write("## How to use this report for SNP-calling decisions\n\n")
    out.write("1. Check whether the **sample missingness tail** is driven by a few poor samples or is broadly spread across the dataset.  \n")
    out.write("2. Check whether **locus missingness** is acceptable for the intended downstream analysis.  \n")
    out.write("3. Compare **mean depth distributions** to spot low-depth noise or high-depth repeat-like loci.  \n")
    out.write("4. Review the **MAF distribution** to decide whether raw settings are retaining too many marginal alleles.  \n")
    out.write("5. Use the **PCA** to see whether raw SNP structure is plausible or dominated by poor-quality samples.\n\n")

    out.write("## Reproducibility notes\n\n")
    out.write(f"- Caller: `{caller}`\n")
    out.write(f"- Raw VCF: `{raw_vcf.name}`\n")
    out.write(f"- Input stats summary: `{stats_txt.name}`\n")
    out.write(f"- Generated report: `{out_md.name}`\n")
PY
    """
}
