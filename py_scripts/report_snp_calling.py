#!/usr/bin/env python3
import csv
import gzip
import math
import statistics
import sys
from pathlib import Path
from typing import Iterable, Optional

EXPECTED_ARG_COUNT = 21
QC_REL_DIR = Path("snp_qc")

NO_FILE_SENTINELS = {"NO_FILE", "NO_DECOMPOSED_VCF", "", "null", "None"}

PERCENT_COLUMNS = {
    "pct_loci_missing",
    "pct_samples_missing",
    "pct_subset_missing",
}

INTEGER_COLUMNS = {
    "number_loci",
    "number_samples_with_locus",
}

THREE_DECIMAL_COLUMNS = {
    "maf",
}

TWO_DECIMAL_COLUMNS = {
    "mean_depth_called",
    "mean_depth",
}


def fail(msg: str) -> None:
    raise SystemExit(msg)


def is_missing(path: Path) -> bool:
    """True when a path is a placeholder rather than a real file."""
    return str(path.name) in NO_FILE_SENTINELS or not path.exists()

def is_sentinel(path) -> bool:
    """For link strings passed as val, which are never staged on disk."""
    return str(Path(path).name) in NO_FILE_SENTINELS

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
        return gzip.open(path, "rt", encoding="utf-8")
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


def to_int(value) -> Optional[int]:
    x = to_float(value)
    return int(x) if x is not None else None


def numeric_values(rows, column: str) -> list:
    values = []
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


def format_table_value(column: str, value) -> str:
    if value in (None, "", "NA"):
        return "NA"

    if isinstance(value, str) and value.upper() in {"TRUE", "FALSE"}:
        return value

    x = to_float(value)
    if x is None:
        return str(value)

    if column in PERCENT_COLUMNS:
        return pct(x, digits=1)
    if column in INTEGER_COLUMNS:
        return num(x, digits=0)
    if column in THREE_DECIMAL_COLUMNS:
        return num(x, digits=3)
    if column in TWO_DECIMAL_COLUMNS:
        return num(x, digits=2)

    if math.isclose(x, round(x), abs_tol=1e-12):
        return num(x, digits=0)
    return num(x, digits=2)


def md_table(rows, columns, max_rows: int = 10) -> str:
    rows = rows[:max_rows]
    if not rows or not columns:
        return "_No rows available._\n"
    header = "| " + " | ".join(columns) + " |\n"
    sep = "| " + " | ".join(["---"] * len(columns)) + " |\n"
    body = ""
    for row in rows:
        formatted = [format_table_value(c, row.get(c, "")) for c in columns]
        body += "| " + " | ".join(formatted) + " |\n"
    return header + sep + body


def derive_angsd_companion_paths(raw_vcf_link: Path):
    path_str = str(raw_vcf_link)
    if path_str.endswith('.vcf.gz'):
        base = path_str[:-7]
    elif path_str.endswith('.vcf'):
        base = path_str[:-4]
    else:
        base = str(raw_vcf_link.with_suffix(''))
    return Path(f"{base}.beagle.gz"), Path(f"{base}.mafs.gz")


def parse_bcftools_stats(stats_txt: Path):
    """Return (SN summary dict, TSTV dict). Empty on any failure."""
    summary_numbers = {}
    ts_tv = None
    if is_missing(stats_txt):
        return summary_numbers, ts_tv
    try:
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
    except Exception:
        pass
    return summary_numbers, ts_tv


def parse_key_value_tsv(path: Path) -> dict:
    """Parse the two-column decompose_summary.tsv."""
    out = {}
    if is_missing(path):
        return out
    try:
        with open(path, "r", encoding="utf-8") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    out[parts[0]] = parts[1]
    except Exception:
        pass
    return out


# ---------------------------------------------------------------------------
# Callset files section
# ---------------------------------------------------------------------------

def build_callset_files_section(caller: str,
                                raw_vcf_link: Path,
                                decomposed_vcf_link: Path,
                                raw_summary: dict,
                                decomp_counts: dict) -> str:
    """Describe both published callsets and which one to use for what."""
    if is_sentinel(decomposed_vcf_link):
        # ANGSD path, or FreeBayes run without decomposition.
        return ""

    raw_records = raw_summary.get("number of records")
    decomp_records = decomp_counts.get("decomposed_records")

    lines = ["## Callset files\n"]
    lines.append(
        "Two callsets are published from this run. They contain the same "
        "underlying evidence in different representations, and they are not "
        "interchangeable.\n"
    )
    lines.append("| File | Records | Representation | Use it for |")
    lines.append("| --- | --- | --- | --- |")
    lines.append(
        f"| [`{raw_vcf_link.name}`]({rel_link_root(raw_vcf_link)}) "
        f"| {num(raw_records, 0)} "
        f"| FreeBayes-native. Multiallelic, MNP and complex records intact. "
        f"| Inspecting the collapsed-paralog signal; re-decomposing under "
        f"different settings. |"
    )
    lines.append(
        f"| [`{decomposed_vcf_link.name}`]({rel_link_root(decomposed_vcf_link)}) "
        f"| {num(decomp_records, 0)} "
        f"| Multiallelics split, MNPs and complex alleles atomized to single "
        f"SNPs, duplicates removed, INFO tags recomputed. "
        f"| All downstream filtering and analysis. Every QC panel in this "
        f"report is computed from this file. |"
    )
    lines.append("")
    lines.append(
        "Record counts differ in both directions. Splitting a multiallelic "
        "site and atomizing an MNP both increase the count, while exact "
        "deduplication decreases it, so the totals are not a simple ratio.\n"
    )
    lines.append(
        "Note that the missingness figures in this report are not comparable "
        "to a run that summarized the raw callset, because the denominator is "
        "now the decomposed record count.\n"
    )
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Collapsed-paralog signal section
# ---------------------------------------------------------------------------

def build_paralog_section(raw_summary: dict, decomp_counts: dict) -> str:
    """Surface the multiallelic / MNP / complex composition of the raw callset.

    In a post-whole-genome-duplication genome assembled from short reads,
    paralogous copies collapse into a single reference locus. Reads from both
    copies then pile onto it, and fixed differences between the copies present
    as multiallelic sites that are heterozygous in nearly every sample. That
    signal only exists per-site in the raw callset; decomposition flattens it.
    """
    if not raw_summary and not decomp_counts:
        return ""

    n_records = to_int(raw_summary.get("number of records")) \
        or to_int(decomp_counts.get("raw_records"))
    n_multi = to_int(raw_summary.get("number of multiallelic sites")) \
        or to_int(decomp_counts.get("raw_multiallelic"))
    n_multi_snp = to_int(raw_summary.get("number of multiallelic SNP sites"))
    n_mnp = to_int(raw_summary.get("number of MNPs")) \
        or to_int(decomp_counts.get("raw_mnp"))
    n_indel = to_int(raw_summary.get("number of indels")) \
        or to_int(decomp_counts.get("raw_indel"))
    n_other = to_int(raw_summary.get("number of others")) \
        or to_int(decomp_counts.get("raw_complex"))

    if n_records is None:
        return ""

    def frac(x):
        if x is None or not n_records:
            return "NA"
        return pct(x / n_records)

    lines = ["## Collapsed-paralog signal (raw callset)\n"]
    lines.append(
        "In a post-duplication genome assembled from short reads, paralogous "
        "copies collapse into one reference locus. Reads from both copies pile "
        "onto it, and fixed differences between the copies appear as "
        "multiallelic sites that are heterozygous in nearly every sample. "
        "These counts come from the raw callset, because decomposition "
        "flattens the per-site evidence.\n"
    )
    lines.append("| Record class | Count | % of raw records |")
    lines.append("| --- | --- | --- |")
    lines.append(f"| Total records | {num(n_records, 0)} | 100.0% |")
    lines.append(f"| Multiallelic sites | {num(n_multi, 0)} | {frac(n_multi)} |")
    if n_multi_snp is not None:
        lines.append(f"| Multiallelic SNP sites | {num(n_multi_snp, 0)} | {frac(n_multi_snp)} |")
    lines.append(f"| MNPs | {num(n_mnp, 0)} | {frac(n_mnp)} |")
    lines.append(f"| Indels | {num(n_indel, 0)} | {frac(n_indel)} |")
    lines.append(f"| Complex / other | {num(n_other, 0)} | {frac(n_other)} |")
    lines.append("")

    notes = []
    if n_multi is not None and n_records:
        multi_frac = n_multi / n_records
        if multi_frac > 0.10:
            notes.append(
                f"**{pct(multi_frac)} of raw records are multiallelic.** That is "
                f"high enough to warrant treating collapsed paralogy as a "
                f"first-order concern rather than a tail effect. Sites that are "
                f"multiallelic, high-depth, and heterozygous in nearly every "
                f"sample are the diagnostic combination."
            )
        elif multi_frac > 0.03:
            notes.append(
                f"{pct(multi_frac)} of raw records are multiallelic, which is a "
                f"moderate level. Worth checking whether they cluster on "
                f"particular contigs before deciding how to treat them."
            )
        else:
            notes.append(
                f"{pct(multi_frac)} of raw records are multiallelic, which is "
                f"low and does not by itself suggest widespread collapsed "
                f"paralogy."
            )

    notes.append(
        "Multiallelic records were split rather than discarded, so they are "
        "still present in the decomposed callset as separate biallelic rows. "
        "Filtering them out requires going back to the raw callset to identify "
        "which positions they came from."
    )
    notes.append(
        "Allele balance at high-depth heterozygous sites distinguishes "
        "collapsed paralogs from true polyploidy: paralogy is locus-specific, "
        "high-depth, and heterozygous across nearly all samples, whereas a "
        "ploidy difference is genome-wide and consistent within a sample."
    )

    for n in notes:
        lines.append(f"- {n}")
    lines.append("")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Ts/Tv sweep section
# ---------------------------------------------------------------------------

def build_sweep_section(tstv_sweep: Path) -> str:
    """Render the Ts/Tv vs stringency grid and interpret its trajectory."""
    if is_missing(tstv_sweep):
        return ""

    rows = safe_read_tsv(tstv_sweep)
    if not rows:
        return ""

    lines = ["## Ts/Tv versus filter stringency\n"]
    lines.append(
        "Each base has one transition partner and two transversion partners, "
        "so a callset dominated by sequencing error converges on Ts/Tv = 0.5 "
        "regardless of organism. Real polymorphism sits well above it. This "
        "floor is universal; there is no universal ceiling, so no target value "
        "is assumed here. Human whole-genome data runs near 2.0 and exome data "
        "near 3.0 to 3.3, the latter because coding sequence is CpG-rich and "
        "5-methylcytosine deamination drives transitions. Reduced-"
        "representation data anchored on a CpG-containing cut site shows the "
        "same elevation.\n"
    )
    lines.append(
        "**Read the trajectory, not any single cell.** Ts/Tv rising with "
        "stringency means the callset is noisy but recoverable by filtering. "
        "Ts/Tv flat near 0.5 even among the most widely shared, common-allele "
        "sites points upstream, to mapping or reference problems that no caller "
        "setting would fix.\n"
    )

    cols = ["min_ns_frac", "min_ns", "min_maf", "n_sites", "n_ts", "n_tv",
            "tstv", "reliable"]
    header_labels = ["Min call rate", "Min NS", "Min MAF", "Sites",
                     "Ts", "Tv", "Ts/Tv", "Reliable"]

    lines.append("| " + " | ".join(header_labels) + " |")
    lines.append("| " + " | ".join(["---"] * len(header_labels)) + " |")
    for row in rows:
        ns_frac = to_float(row.get("min_ns_frac"))
        cells = [
            pct(ns_frac, digits=0) if ns_frac is not None else "NA",
            num(row.get("min_ns"), 0),
            str(row.get("min_maf", "NA")),
            num(row.get("n_sites"), 0),
            num(row.get("n_ts"), 0),
            num(row.get("n_tv"), 0),
            str(row.get("tstv", "NA")),
            str(row.get("reliable", "NA")),
        ]
        lines.append("| " + " | ".join(cells) + " |")
    lines.append("")

    # --- interpretation ----------------------------------------------------
    reliable = [r for r in rows if r.get("reliable") == "yes"]
    unreliable = [r for r in rows if r.get("reliable") == "no"]

    notes = []

    unfiltered = next(
        (r for r in rows
         if to_float(r.get("min_ns_frac")) == 0.0
         and to_float(r.get("min_maf")) == 0.0),
        None,
    )

    if unfiltered:
        base_tstv = to_float(unfiltered.get("tstv"))
        if base_tstv is not None:
            notes.append(
                f"Unfiltered biallelic SNPs: Ts/Tv = **{num(base_tstv, 2)}** on "
                f"{num(unfiltered.get('n_sites'), 0)} sites."
            )
            n_ts = to_float(unfiltered.get("n_ts")) or 0.0
            n_tv = to_float(unfiltered.get("n_tv")) or 0.0
            total = n_ts + n_tv
            if total > 0:
                ts_frac = n_ts / total
                estimates = []
                for assumed in (1.5, 2.0, 2.5, 3.0):
                    real_ts = assumed / (assumed + 1.0)
                    denom = real_ts - (1.0 / 3.0)
                    if denom <= 0:
                        continue
                    f_real = (ts_frac - (1.0 / 3.0)) / denom
                    f_real = max(0.0, min(1.0, f_real))
                    estimates.append((assumed, f_real, f_real * total))
                if estimates:
                    lo = min(e[1] for e in estimates)
                    hi = max(e[1] for e in estimates)
                    notes.append(
                        f"Treating the callset as a mixture of error sites at "
                        f"Ts/Tv 0.5 and real sites at an assumed 1.5 to 3.0, the "
                        f"implied real fraction is **{pct(lo)} to {pct(hi)}**, or "
                        f"roughly {num(min(e[2] for e in estimates), 0)} to "
                        f"{num(max(e[2] for e in estimates), 0)} real SNPs. This "
                        f"is a rough decomposition and the spread across assumed "
                        f"values is the uncertainty, not a confidence interval."
                    )

    if reliable:
        best = max(reliable, key=lambda r: to_float(r.get("tstv")) or -1.0)
        best_tstv = to_float(best.get("tstv"))
        best_ns_frac = to_float(best.get("min_ns_frac"))
        notes.append(
            f"Highest Ts/Tv among sufficiently populated cells: "
            f"**{num(best_tstv, 2)}** at a "
            f"{pct(best_ns_frac, digits=0)} call-rate floor and MAF >= "
            f"{best.get('min_maf')}, on {num(best.get('n_sites'), 0)} sites."
        )
        if best_tstv is not None:
            if best_tstv >= 1.5:
                notes.append(
                    "This is comfortably clear of the 0.5 error floor, so real "
                    "polymorphism is present and the callset is a filtering "
                    "problem rather than a calling problem. Re-running variant "
                    "calling with different parameters is unlikely to help."
                )
            elif best_tstv >= 1.0:
                notes.append(
                    "This is above the error floor but below what a diverse "
                    "outbreeding organism would typically show, so real signal "
                    "is present but diluted. Worth pushing the stringency "
                    "further before concluding anything about the data."
                )
            else:
                notes.append(
                    "**This stays close to the 0.5 error floor even in the most "
                    "stringent populated cell.** That is the signature of a "
                    "problem upstream of variant calling. Investigate mapping, "
                    "reference quality and chimeric reads before adjusting any "
                    "caller parameters, because no parameter set fixes it."
                )

    if unreliable:
        notes.append(
            f"{len(unreliable)} of {len(rows)} grid cells fell below the site-"
            f"count threshold and are marked unreliable. Ts/Tv on a few hundred "
            f"sites carries roughly plus or minus 0.2 to 0.3, and below about "
            f"200 sites it stops being informative."
        )

    notes.append(
        "The MAF axis here is a **diagnostic** axis, not a filtering "
        "recommendation. A MAF floor is appropriate for structure and PCA and "
        "destructive for anything that reads the site frequency spectrum, so "
        "the pipeline deliberately publishes no pre-filtered callset."
    )

    for n in notes:
        lines.append(f"- {n}")
    lines.append("")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# Main report
# ---------------------------------------------------------------------------

def build_report(
    prefix: str,
    caller: str,
    raw_vcf_link: Path,
    decomposed_vcf_link: Path,
    raw_stats: Path,
    tstv_sweep: Path,
    decompose_summary: Path,
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
    raw_summary, raw_ts_tv = parse_bcftools_stats(raw_stats)
    decomp_counts = parse_key_value_tsv(decompose_summary)

    sample_rows = safe_read_tsv(sample_qc_derived)
    locus_rows = safe_read_tsv(locus_qc_derived)
    worst_sample_rows = safe_read_tsv(worst_samples)
    worst_locus_rows = safe_read_tsv(worst_loci)

    has_decomposed = not is_sentinel(decomposed_vcf_link)

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
        scope = "decomposed callset" if has_decomposed else "raw callset"
        interpretation.append(
            f"The {scope} contains **{num(n_records, 0)} variant records** across "
            f"**{num(summary_numbers.get('number of samples', metrics.get('samples_total', 'NA')), 0)} samples**."
        )

    if ts_tv and ts_tv.get("ratio") not in (None, "", "NA"):
        interpretation.append(
            f"The **transition/transversion ratio is {num(ts_tv['ratio'], 2)}**. "
            f"An error-dominated callset converges on 0.5, so this is the "
            f"single most informative summary number here. See the Ts/Tv sweep "
            f"section for how it behaves under filtering."
        )

    if metrics.get("sample_missing_gt30") is not None:
        interpretation.append(
            f"At the sample level, **{num(metrics['sample_missing_gt30'], 0)} samples exceed 30% missing data**, with median sample missingness of **{pct(metrics.get('sample_missing_median'))}**."
        )
        if metrics.get("samples_total") and \
                metrics["sample_missing_gt30"] == metrics["samples_total"]:
            interpretation.append(
                "Every sample exceeds that threshold, which usually means the "
                "denominator rather than the samples is the problem: a callset "
                "still dominated by low-call-rate false positives makes every "
                "sample look poor. Do not cull samples on this number. Judge "
                "them on a core set instead, using the Ts/Tv sweep to pick one."
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

    beagle_path = None
    mafs_path = None
    raw_call_links = [f"[`{raw_vcf_link.name}`]({rel_link_root(raw_vcf_link)})"]
    if has_decomposed:
        raw_call_links.append(
            f"[`{decomposed_vcf_link.name}`]({rel_link_root(decomposed_vcf_link)})"
        )
    if caller.lower() == "angsd":
        beagle_path, mafs_path = derive_angsd_companion_paths(raw_vcf_link)
        raw_call_links.append(f"[`{beagle_path.name}`]({rel_link_root(beagle_path)})")
        raw_call_links.append(f"[`{mafs_path.name}`]({rel_link_root(mafs_path)})")
    raw_call_files_line = " | ".join(raw_call_links)

    reproducibility_notes = [
        f"- Caller: `{caller}`",
        f"- Raw VCF: `{raw_vcf_link.name}`",
    ]
    if has_decomposed:
        reproducibility_notes.append(
            f"- Decomposed VCF: `{decomposed_vcf_link.name}` (summarized in this report)"
        )
    if beagle_path is not None and mafs_path is not None:
        reproducibility_notes.append(f"- Beagle file: `{beagle_path.name}`")
        reproducibility_notes.append(f"- MAF file: `{mafs_path.name}`")
    reproducibility_notes.append(f"- Input stats summary: `{stats_txt.name}`")
    if not is_missing(raw_stats):
        reproducibility_notes.append(f"- Raw-callset stats: `{raw_stats.name}`")
    if not is_missing(tstv_sweep):
        reproducibility_notes.append(f"- Ts/Tv sweep: `{tstv_sweep.name}`")
    reproducibility_notes.append(f"- Generated report: `{out_md.name}`")
    reproducibility_notes_text = "\n".join(reproducibility_notes)

    summarized_note = (
        f"**Summarized file:** [`{decomposed_vcf_link.name}`]"
        f"({rel_link_root(decomposed_vcf_link)}) (decomposed)  \n"
        if has_decomposed else ""
    )

    report = f"""# SNP calling summary report: `{prefix}`

**Caller:** `{caller}`  
**Raw call files:** {raw_call_files_line}  
{summarized_note}**Primary VCF summary file:** [`{stats_txt.name}`]({rel_link_qc(stats_txt)})  
**Standardized summary text:** [`{standardized_summary.name}`]({rel_link_qc(standardized_summary)})

## Main takeaways

"""

    for line in interpretation + caller_specific:
        report += f"- {line}\n"

    report += "\n\n"
    report += build_callset_files_section(
        caller, raw_vcf_link, decomposed_vcf_link, raw_summary, decomp_counts
    )
    report += build_sweep_section(tstv_sweep)
    report += build_paralog_section(raw_summary, decomp_counts)

    report += f"""## Core summary metrics

Computed on the {'decomposed' if has_decomposed else 'raw'} callset.

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

Missingness and depth in this panel are measured on the PCA site subset, not
the full callset, because `glPca` mean-imputes missing genotypes and it is
subset missingness that pulls samples toward the origin. Site selection
requires both a call-rate floor and a MAF floor; the selection log is recorded
in the standardized summary.

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
"""

    if not is_missing(tstv_sweep):
        report += (
            f"| [`{tstv_sweep.name}`]({rel_link_qc(tstv_sweep)}) "
            f"| Ts/Tv across the threshold grid. Start here when choosing "
            f"filtering thresholds. |\n"
        )

    report += f"""| [`{sample_qc_derived.name}`]({rel_link_qc(sample_qc_derived)}) | Derived per-sample QC summary used for ranking poor samples. |
| [`{locus_qc_derived.name}`]({rel_link_qc(locus_qc_derived)}) | Derived per-locus QC summary used for missingness/depth/MAF interpretation. |
| [`{missing_indv.name}`]({rel_link_qc(missing_indv)}) | Raw per-sample missingness from vcftools. |
| [`{missing_site.name}`]({rel_link_qc(missing_site)}) | Raw per-site missingness from vcftools. |
| [`{freq_tsv.name}`]({rel_link_qc(freq_tsv)}) | Frequency table for checking the MAF distribution of raw SNP calls. |
| [`{site_qc_tsv.name}`]({rel_link_qc(site_qc_tsv)}) | Standardized site-level QC table; especially important for ANGSD outputs. |
| [`{sample_qc_tsv.name}`]({rel_link_qc(sample_qc_tsv)}) | Standardized sample-level QC table. |
| [`{worst_samples.name}`]({rel_link_qc(worst_samples)}) | Quick list of the worst samples to investigate first. |
| [`{worst_loci.name}`]({rel_link_qc(worst_loci)}) | Quick list of the worst loci to investigate first. |

## How to use this report for SNP-calling decisions

1. Read the **Ts/Tv sweep** first. It answers whether this is a filtering problem or an upstream problem, and nothing else in the report is interpretable until that is settled.  
2. Check whether the **sample missingness tail** is driven by a few poor samples or is broadly spread. If every sample looks poor, suspect the denominator rather than the samples.  
3. Check whether **locus missingness** is acceptable for the intended downstream analysis, remembering that reduced-representation data has structural dropout that whole-genome data does not.  
4. Compare **mean depth distributions** to spot low-depth noise and high-depth repeat-like loci. A bimodal locus-depth distribution usually separates real loci from noise.  
5. Review the **MAF distribution**, keeping in mind that MAF is uninterpretable without a call-rate floor beneath it, since AN scales with the number of samples called.  
6. Use the **PCA** to see whether structure is plausible or still dominated by data quantity. Axes that track missingness are an artifact, not a result.  
7. Check the **collapsed-paralog signal** before setting depth ceilings, since paralogous pileups and genuinely high-coverage loci require different treatment.

## Reproducibility notes

{reproducibility_notes_text}
"""

    with open(out_md, "w", encoding="utf-8") as out:
        out.write(report)

    return out_md


def main(argv: list) -> int:
    if len(argv) != EXPECTED_ARG_COUNT:
        fail(
            f"Expected {EXPECTED_ARG_COUNT} arguments, received {len(argv)}: {argv}"
        )

    (
        prefix,
        caller,
        raw_vcf_link,
        decomposed_vcf_link,
        raw_stats,
        tstv_sweep,
        decompose_summary,
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
        decomposed_vcf_link=Path(decomposed_vcf_link),
        raw_stats=Path(raw_stats),
        tstv_sweep=Path(tstv_sweep),
        decompose_summary=Path(decompose_summary),
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
