// modules/tstv_sweep.nf - Ts/Tv as a function of filter stringency
//
// Diagnostic, not a filter. This process emits no filtered VCF, on purpose:
// a MAF floor is appropriate for structure and PCA and destructive for
// anything reading the site frequency spectrum, so the thresholds cannot be
// chosen once at the pipeline level. It measures whether real polymorphism
// exists under the false-positive mass and leaves the decision downstream.
//
// Why Ts/Tv works as the readout: each base has one transition partner and
// two transversion partners, so an error-dominated callset converges on
// Ts/Tv = 0.5 regardless of organism. Real polymorphism runs well above it.
// The 0.5 floor is universal; the ceiling is not. Human WGS is ~2.0, human
// exome ~3.0-3.3 because coding sequence is CpG-rich and 5mC deamination
// drives transitions. ddRAD anchored on a CpG-containing cut site shows the
// same elevation. No target value is hardcoded anywhere here.
//
// Single pass over the VCF. The threshold grid is applied in awk against an
// extracted table rather than by re-filtering the VCF once per cell, which
// at WGS scale is the difference between minutes and days.

process TSTV_SWEEP {
    tag "tstv_sweep"
    publishDir "${params.output_dir}/snp_qc", mode: 'copy'

    input:
    path vcf
    path vcf_index
    val data_type

    output:
    path "${vcf.simpleName}.tstv_sweep.tsv",       emit: sweep
    path "${vcf.simpleName}.sweep_sites.tsv.gz",   emit: sites
    path "${vcf.simpleName}.tstv_sweep.txt",       emit: report

    script:
    // Grids differ by data type. ddRAD has structural locus dropout from
    // restriction-site polymorphism and allele dropout, so call rates below
    // 50% are normal and the grid has to reach down there. WGS at high sample
    // counts should not, and its dominant false-positive class is collapsed-
    // repeat pileups rather than dropout, so the grid stays high and the MAF
    // axis extends lower to keep rare variants in view.
    def ns_default = data_type == 'wgs' ?
        '0,0.50,0.80,0.90,0.95' :
        '0,0.25,0.50,0.70,0.80,0.90'
    def maf_default = data_type == 'wgs' ?
        '0,0.001,0.01,0.05' :
        '0,0.01,0.05,0.10'

    def ns_grid  = params.tstv_ns_grid  != null ? params.tstv_ns_grid  : ns_default
    def maf_grid = params.tstv_maf_grid != null ? params.tstv_maf_grid : maf_default

    // Rows built from fewer than this many sites carry too much sampling
    // error to read. Flagged rather than suppressed.
    def min_reliable = params.tstv_min_sites != null ? params.tstv_min_sites : 200

    """
    set -euo pipefail

    N_SAMPLES=\$(bcftools query -l "${vcf}" | wc -l)
    echo "\${N_SAMPLES}" > n_samples.txt
    echo "Samples: \${N_SAMPLES}"
    echo "Data type: ${data_type}"
    echo "NS grid (fractions): ${ns_grid}"
    echo "MAF grid: ${maf_grid}"

    if [ "\${N_SAMPLES}" -eq 0 ]; then
        echo "ERROR: no samples in ${vcf}" >&2
        exit 1
    fi

    # -----------------------------------------------------------------------
    # 1. Single extraction pass
    #
    # Restricted to biallelic SNPs. After decomposition every record is
    # already biallelic, so -m2 -M2 mainly guards against star alleles if
    # --atom-overlaps '*' was used upstream.
    #
    # DP is carried through unused by the sweep so a depth-axis sweep can be
    # done later from sweep_sites.tsv.gz without another pass over the VCF.
    # That is the natural WGS extension: mean depth per called sample
    # (DP/NS) is the collapsed-repeat detector.
    # -----------------------------------------------------------------------
    bcftools view -v snps -m2 -M2 --threads ${task.cpus} -Ou "${vcf}" \
      | bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t%INFO/NS\\t%INFO/MAF\\t%INFO/DP\\n' \
      > sweep_sites.tsv

    N_SITES=\$(wc -l < sweep_sites.tsv)
    echo "Biallelic SNPs extracted: \${N_SITES}"

    # -----------------------------------------------------------------------
    # 2. Bin against the grid
    # -----------------------------------------------------------------------
    awk -F'\\t' \
        -v ns_grid="${ns_grid}" \
        -v maf_grid="${maf_grid}" \
        -v n_samples="\${N_SAMPLES}" \
        -v min_reliable="${min_reliable}" \
        -v out_tsv="${vcf.simpleName}.tstv_sweep.tsv" '
    function is_transition(r, a) {
        if ((r == "A" && a == "G") || (r == "G" && a == "A")) return 1
        if ((r == "C" && a == "T") || (r == "T" && a == "C")) return 1
        return 0
    }
    BEGIN {
        n_ns  = split(ns_grid,  NSF, ",")
        n_maf = split(maf_grid, MAF, ",")
        for (i = 1; i <= n_ns; i++) NSMIN[i] = NSF[i] * n_samples
        skipped = 0
    }
    {
        ref = toupper(\$3); alt = toupper(\$4)
        # Only unambiguous single-base substitutions are classifiable.
        if (length(ref) != 1 || length(alt) != 1) { skipped++; next }
        if (ref !~ /^[ACGT]\$/ || alt !~ /^[ACGT]\$/) { skipped++; next }

        ns = \$5 + 0
        maf = (\$6 == "." || \$6 == "") ? -1 : \$6 + 0
        if (maf < 0) { skipped++; next }

        ts = is_transition(ref, alt)

        for (i = 1; i <= n_ns; i++) {
            if (ns < NSMIN[i]) continue
            for (j = 1; j <= n_maf; j++) {
                if (maf < MAF[j] + 0) continue
                if (ts) TS[i, j]++; else TV[i, j]++
            }
        }
    }
    END {
        printf "min_ns_frac\\tmin_ns\\tmin_maf\\tn_sites\\tn_ts\\tn_tv\\ttstv\\treliable\\n" > out_tsv
        printf "%-12s %-9s %-9s %14s %12s %12s %8s %10s\\n",
               "MIN_NS_FRAC", "MIN_NS", "MIN_MAF", "n_sites", "n_Ts", "n_Tv",
               "Ts/Tv", "reliable"

        for (i = 1; i <= n_ns; i++) {
            for (j = 1; j <= n_maf; j++) {
                nts = TS[i, j] + 0
                ntv = TV[i, j] + 0
                tot = nts + ntv
                if (ntv > 0) { tstv = sprintf("%.3f", nts / ntv) } else { tstv = "NA" }
                rel = (tot >= min_reliable) ? "yes" : "no"

                printf "%.4f\\t%d\\t%s\\t%d\\t%d\\t%d\\t%s\\t%s\\n",
                       NSF[i] + 0, NSMIN[i], MAF[j], tot, nts, ntv, tstv, rel > out_tsv
                printf "%-12.4f %-9d %-9s %14d %12d %12d %8s %10s\\n",
                       NSF[i] + 0, NSMIN[i], MAF[j], tot, nts, ntv, tstv, rel
            }
        }
        printf "\\nUnclassifiable records skipped: %d\\n", skipped
    }' sweep_sites.tsv > sweep_console.txt

    cat sweep_console.txt

    # -----------------------------------------------------------------------
    # 3. Narrative summary
    #
    # The mixture estimate below assumes error sites sit at Ts/Tv = 0.5 and
    # real sites at an assumed value, then solves for the real fraction. It is
    # a rough decomposition, sensitive to the assumed real Ts/Tv, and is
    # reported with that sensitivity rather than as a point estimate.
    # -----------------------------------------------------------------------
    python3 <<'PY' > ${vcf.simpleName}.tstv_sweep.txt
import csv, sys

sweep_path = "${vcf.simpleName}.tstv_sweep.tsv"
n_samples  = int(open("n_samples.txt").read().strip())

with open(sweep_path) as fh:
    rows = list(csv.DictReader(fh, delimiter="\\t"))

print("=== Ts/Tv VS FILTER STRINGENCY ===")
print(f"Data type: ${data_type}")
print(f"Samples: {n_samples}")
print()
print("An error-dominated callset converges on Ts/Tv = 0.5, because each base")
print("has one transition partner and two transversion partners. Real")
print("polymorphism runs well above that. Rising Ts/Tv with stringency means")
print("the callset is noisy but recoverable by filtering. Flat Ts/Tv near 0.5")
print("even in the most-shared, common-allele core points upstream, to")
print("mapping or reference problems that no caller setting would fix.")
print()

def f(row, key, default=0.0):
    try:
        return float(row[key])
    except (ValueError, KeyError, TypeError):
        return default

unfiltered = next((r for r in rows
                   if f(r, "min_ns_frac") == 0.0 and f(r, "min_maf") == 0.0), None)

if unfiltered and f(unfiltered, "n_sites") > 0:
    nts, ntv = f(unfiltered, "n_ts"), f(unfiltered, "n_tv")
    tot = nts + ntv
    ts_frac = nts / tot if tot else 0.0
    print(f"[Unfiltered biallelic SNPs] n={int(tot):,}  "
          f"Ts/Tv={unfiltered['tstv']}  Ts fraction={ts_frac:.4f}")

    # Ts fraction: error = 1/3 (Ts/Tv 0.5); real = r/(r+1) for Ts/Tv = r.
    print()
    print("Rough mixture estimate of the real-variant fraction, assuming")
    print("error sites at Ts/Tv 0.5:")
    for assumed in (1.5, 2.0, 2.5, 3.0):
        real_ts = assumed / (assumed + 1.0)
        denom = real_ts - (1.0 / 3.0)
        if denom <= 0:
            continue
        frac = (ts_frac - (1.0 / 3.0)) / denom
        frac = max(0.0, min(1.0, frac))
        print(f"  if real Ts/Tv = {assumed:.1f}  ->  {frac*100:5.1f}% real  "
              f"(~{int(frac*tot):,} sites)")
    print()
    print("Treat the spread across assumed values as the uncertainty. If the")
    print("estimate is near zero, the mixture model is not the right lens and")
    print("the problem is more likely upstream than statistical.")
else:
    print("No unfiltered row available; skipping mixture estimate.")

print()
print("[Trajectory across the grid]")
reliable = [r for r in rows if r.get("reliable") == "yes"]
if reliable:
    best = max(reliable, key=lambda r: f(r, "tstv", -1))
    print(f"Highest Ts/Tv among rows with >= ${min_reliable} sites: "
          f"{best['tstv']} at NS>={best['min_ns']} "
          f"({float(best['min_ns_frac'])*100:.0f}% call rate), "
          f"MAF>={best['min_maf']}, n={int(f(best,'n_sites')):,}")
else:
    print("No grid cell reached the reliability threshold of ${min_reliable} sites.")
    print("Either the callset is very small or the grid is too stringent for it.")

n_unreliable = sum(1 for r in rows if r.get("reliable") == "no")
if n_unreliable:
    print()
    print(f"{n_unreliable} of {len(rows)} grid cells fell below "
          f"${min_reliable} sites and are marked unreliable. Ts/Tv on a few "
          f"hundred sites carries roughly +/- 0.2-0.3; below ~200 it stops "
          f"being informative.")

print()
print("[Full grid]")
with open(sweep_path) as fh:
    sys.stdout.write(fh.read())
PY

    cat ${vcf.simpleName}.tstv_sweep.txt

    bgzip -f -@ ${task.cpus} sweep_sites.tsv
    mv sweep_sites.tsv.gz ${vcf.simpleName}.sweep_sites.tsv.gz
    """
}
