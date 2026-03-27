#!/usr/bin/env python3
import csv
import gzip
import math
import statistics
import sys
from pathlib import Path
from typing import Iterable, Optional

EXPECTED_ARG_COUNT = 17
QC_REL_DIR = Path("snp_qc")


def fail(msg: str) -> None:
    raise SystemExit(msg)


def rel_link_qc(path: Path) -> str:
    return str(QC_REL_DIR / path.name)


def rel_link_root(path: Path) -> str:
    return path.name


def pct(x: Optional[float], digits: int = 1) -> str:
    if x is None:
        return "NA"
    try:
        return f"{100 * float(x):,.{digits}f}%"
    except Exception:
        return "NA"


def num(x, digits: int = 2) -> str:
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
    return open(path, "r", encoding="utf-8")


def read_tsv(path: Path):
    with open_text_auto(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader)


def safe_read_tsv(path: Path):
    try:
        return read_tsv(path)
    except Exception:
        return []


def to_float(value) -> Optional[float]:
    if value in (None, "", "NA", "nan", "NaN"):
        return None
    try:
        x = float(value)
    except Exception:
        return None
    return x if math.isfinite(x) else None


def numeric_values(rows, column: str) -> list[float]:
    values: list[float] = []
    for row in rows:
        x = to_float(row.get(column))
        if x is not None:
            values.append(x)
    return values


def median_or_none(values: Iterable[float]) -> Optional[float]:
    values = list(values)
    if not values:
        return None
    return statistics.median(values)


def count_gt(values: Iterable[float], threshold: float) -> int:
    return sum(x > threshold for x in values)


def count_lt(values: Iterable[float], threshold: float) -> int:
    return sum(x < threshold for x in values)


def md_table(rows, columns, max_rows: int = 10) -> str:
    rows = rows[:max_rows]
    if not rows or not columns:
        return "_No rows available._\n"
    header = "| " + " | ".join(columns) + " |\n"
    sep = "| " + " | ".join(["---"] * len(columns)) + " |\n"
    body = ""
    for row in rows:
        body += "| " + " | ".join(str(row.get(c, "")) for c in columns) + " |\n"
    return header + sep + body


def parse_bcftools_stats(stats_txt: Path):
    summary_numbers = {}
    ts_tv = None
    with open(stats_txt, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("SN\t"):
                parts = line.rstrip().split("\t")
                if len(parts) >= 4:
                    summary_numbers[parts[2].rstrip(":")] = parts[3]
            elif line.startswith("TSTV\t"):
                parts = line.rstrip().split("\t")
                if len(parts) >= 5:
                    ts_tv = {"ts": parts[2], "tv": parts[3], "ratio": parts[4]}
    return summary_numbers, ts_tv


def build_report(
    prefix: str,
    caller: str,
    raw_vcf_link: Path,
    stats_txt: Path,
    standardized_summary: Path,
    summary_plot: Path,
    extra_plot: Path,
    pca_plot: Path,
    worst_samples: Path,
    worst_loci: Path,
    sample_qc_derived: Path,
    locus_qc_derived: Path,
    missing_indv: Path,
    missing_site: Path,
    freq_tsv: Path,
    site_qc_tsv: Path,
    sample_qc_tsv: Path,
) -> Path:
    out_md = Path(f"{prefix}_snp_calling_report.md")

    summary_numbers, ts_tv = parse_bcftools_stats(stats_txt)

    sample_rows = safe_read_tsv(sample_qc_derived)
    locus_rows = safe_read_tsv(locus_qc_derived)
    worst_sample_rows = safe_read_tsv(worst_samples)
    worst_locus_rows = safe_read_tsv(worst_loci)

    metrics = {}

    if sample_rows and "sample" in sample_rows[0]:
        miss = numeric_values(sample_rows, "pct_loci_missing")
        mean_dp = numeric_values(sample_rows, "mean_depth_called")
        n_loci = numeric_values(sample_rows, "number_loci")

        metrics["samples_total"] = len(sample_rows)
        metrics["sample_missing_median"] = median_or_none(miss)
        metrics["sample_missing_gt10"] = count_gt(miss, 0.10) if miss else None
        metrics["sample_missing_gt20"] = count_gt(miss, 0.20) if miss else None
        metrics["sample_missing_gt30"] = count_gt(miss, 0.30) if miss else None
        metrics["sample_mean_dp_median"] = median_or_none(mean_dp)
        metrics["sample_called_loci_median"] = median_or_none(n_loci)

    if locus_rows and "locus" in locus_rows[0]:
        locus_missing = numeric_values(locus_rows, "pct_samples_missing")
        shared = numeric_values(locus_rows, "number_samples_with_locus")
        locus_depth = numeric_values(locus_rows, "mean_depth")
        maf = numeric_values(locus_rows, "maf")

        metrics["loci_total"] = len(locus_rows)
        metrics["locus_missing_median"] = median_or_none(locus_missing)
        metrics["locus_missing_gt10"] = count_gt(locus_missing, 0.10) if locus_missing else None
        metrics["locus_missing_gt20"] = count_gt(locus_missing, 0.20) if locus_missing else None
        metrics["locus_missing_gt30"] = count_gt(locus_missing, 0.30) if locus_missing else None
        metrics["samples_per_locus_median"] = median_or_none(shared)
        metrics["locus_depth_median"] = median_or_none(locus_depth)
        metrics["maf_median"] = median_or_none(maf)
        metrics["maf_lt01"] = count_lt(maf, 0.01) if maf else None
        metrics["maf_lt05"] = count_lt(maf, 0.05) if maf else None
        metrics["maf_lt10"] = count_lt(maf, 0.10) if maf else None

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

    sample_cols = [
        c
        for c in [
            "sample",
            "pct_loci_missing",
            "number_loci",
            "mean_depth_called",
            "flag_high_missing",
            "flag_low_depth",
        ]
        if worst_sample_rows and c in worst_sample_rows[0]
    ]
    locus_cols = [
        c
        for c in [
            "locus",
            "pct_samples_missing",
            "number_samples_with_locus",
            "mean_depth",
            "maf",
        ]
        if worst_locus_rows and c in worst_locus_rows[0]
    ]

    report = f"""# SNP calling summary report: `{prefix}`

**Caller:** `{caller}`  
**Raw VCF:** [`{raw_vcf_link.name}`]({rel_link_root(raw_vcf_link)})  
**Primary VCF summary file:** [`{stats_txt.name}`]({rel_link_qc(stats_txt)})  
**Standardized summary text:** [`{standardized_summary.name}`]({rel_link_qc(standardized_summary)})

## Main takeaways

"""

    for line in interpretation + caller_specific:
        report += f"- {line}\n"

    report += f"""

## Core summary metrics

| Metric | Value |
| --- | --- |
| Samples | {num(summary_numbers.get('number of samples', metrics.get('samples_total', 'NA')), 0)} |
| Variant records | {num(summary_numbers.get('number of records', metrics.get('loci_total', 'NA')), 0)} |
| SNPs | {num(summary_numbers.get('number of SNPs', 'NA'), 0)} |
| Multiallelic sites | {num(summary_numbers.get('number of multiallelic sites', 'NA'), 0)} |
| Ts/Tv | {num(ts_tv['ratio'], 2) if ts_tv else 'NA'} |
| Median sample missingness | {pct(metrics.get('sample_missing_median'))} |
| Samples >30% missing | {num(metrics.get('sample_missing_gt30'), 0)} |
| Median locus missingness | {pct(metrics.get('locus_missing_median'))} |
| Loci >30% missing | {num(metrics.get('locus_missing_gt30'), 0)} |
| Median samples per locus | {num(metrics.get('samples_per_locus_median'), 0)} |
| Median sample mean depth | {num(metrics.get('sample_mean_dp_median'), 2)} |
| Median locus mean depth | {num(metrics.get('locus_depth_median'), 2)} |
| Median MAF | {num(metrics.get('maf_median'), 3)} |
| Loci with MAF <0.05 | {num(metrics.get('maf_lt05'), 0)} |
| Loci with MAF <0.10 | {num(metrics.get('maf_lt10'), 0)} |

## QC plots

### Summary panel

![Summary plots]({rel_link_qc(summary_plot)})

### Additional diagnostics

![Additional summary plots]({rel_link_qc(extra_plot)})

### PCA

![PCA plot]({rel_link_qc(pca_plot)})

## Candidate samples for review or removal

These are the worst samples by derived sample QC, intended to guide inspection of missingness and depth before downstream filtering.

{md_table(worst_sample_rows, sample_cols, max_rows=15)}
## Candidate loci / locus classes for review

These are the worst loci by derived locus QC, useful for understanding whether raw calls are failing because of sparse representation, unusual depth, or both.

{md_table(worst_locus_rows, locus_cols, max_rows=15)}
## Files to inspect next

| File | Why it matters |
| --- | --- |
| [`{sample_qc_derived.name}`]({rel_link_qc(sample_qc_derived)}) | Derived per-sample QC summary used for ranking poor samples. |
| [`{locus_qc_derived.name}`]({rel_link_qc(locus_qc_derived)}) | Derived per-locus QC summary used for missingness/depth/MAF interpretation. |
| [`{missing_indv.name}`]({rel_link_qc(missing_indv)}) | Raw per-sample missingness from vcftools. |
| [`{missing_site.name}`]({rel_link_qc(missing_site)}) | Raw per-site missingness from vcftools. |
| [`{freq_tsv.name}`]({rel_link_qc(freq_tsv)}) | Frequency table for checking the MAF distribution of raw SNP calls. |
| [`{site_qc_tsv.name}`]({rel_link_qc(site_qc_tsv)}) | Standardized site-level QC table; especially important for ANGSD outputs. |
| [`{sample_qc_tsv.name}`]({rel_link_qc(sample_qc_tsv)}) | Standardized sample-level QC table. |
| [`{worst_samples.name}`]({rel_link_qc(worst_samples)}) | Quick list of the worst samples to investigate first. |
| [`{worst_loci.name}`]({rel_link_qc(worst_loci)}) | Quick list of the worst loci to investigate first. |

## How to use this report for SNP-calling decisions

1. Check whether the **sample missingness tail** is driven by a few poor samples or is broadly spread across the dataset.  
2. Check whether **locus missingness** is acceptable for the intended downstream analysis.  
3. Compare **mean depth distributions** to spot low-depth noise or high-depth repeat-like loci.  
4. Review the **MAF distribution** to decide whether raw settings are retaining too many marginal alleles.  
5. Use the **PCA** to see whether raw SNP structure is plausible or dominated by poor-quality samples.

## Reproducibility notes

- Caller: `{caller}`
- Raw VCF: `{raw_vcf_link.name}`
- Input stats summary: `{stats_txt.name}`
- Generated report: `{out_md.name}`
"""

    with open(out_md, "w", encoding="utf-8") as out:
        out.write(report)

    return out_md


def main(argv: list[str]) -> int:
    if len(argv) != EXPECTED_ARG_COUNT:
        fail(
            f"Expected {EXPECTED_ARG_COUNT} arguments, received {len(argv)}: {argv}"
        )

    (
        prefix,
        caller,
        raw_vcf_link,
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
    ) = argv

    build_report(
        prefix=prefix,
        caller=caller,
        raw_vcf_link=Path(raw_vcf_link),
        stats_txt=Path(stats_txt),
        standardized_summary=Path(standardized_summary),
        summary_plot=Path(summary_plot),
        extra_plot=Path(extra_plot),
        pca_plot=Path(pca_plot),
        worst_samples=Path(worst_samples),
        worst_loci=Path(worst_loci),
        sample_qc_derived=Path(sample_qc_derived),
        locus_qc_derived=Path(locus_qc_derived),
        missing_indv=Path(missing_indv),
        missing_site=Path(missing_site),
        freq_tsv=Path(freq_tsv),
        site_qc_tsv=Path(site_qc_tsv),
        sample_qc_tsv=Path(sample_qc_tsv),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
