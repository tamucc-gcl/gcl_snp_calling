// modules/angsd_chunk.nf - ANGSD genotyping for genomic chunks

process ANGSD_CHUNK {
    tag "angsd_${chunk_id}"

    // let Nextflow control threads; use cpus in the process definition
    cpus 8

    input:
    tuple val(chunk_id), val(regions_string)
    path reference
    path reference_fai
    path bams        // BAM files
    path bais        // BAI files
    path config      // ANGSD configuration JSON (may be missing / null)
    path sites_file  // Optional sites file

    output:
    tuple val(chunk_id), path("${chunk_id}.*"), emit: chunk_files

    script:
    """
    echo "Processing chunk ${chunk_id} with ANGSD"
    # echo "Regions: ${regions_string}"
    echo "Reference: ${reference}"
    touch ${reference}.fai

    # --- Basic BAM checks ----------------------------------------------------
    BAM_COUNT=\$(ls -1 *.bam 2>/dev/null | wc -l)
    echo "Total BAM files: \$BAM_COUNT"

    if [ "\$BAM_COUNT" -eq 0 ]; then
        echo "ERROR: No BAM files found!"
        exit 1
    fi

    echo "Verifying BAM indices:"
    for bam in *.bam; do
        if [ -f "\$bam" ]; then
            if [ -f "\${bam}.bai" ]; then
                echo "  \$bam: index present"
            else
                echo "  \$bam: ERROR - index missing!"
                exit 1
            fi
        fi
    done

    # Create BAM list
    ls -1 *.bam > bam.list
    echo "BAM list created with \$BAM_COUNT files"

    # --- Regions file from regions_string -----------------------------------
    echo "${regions_string}" | tr ',' '\\n' > regions.txt

    > angsd_regions.txt
    while read -r region; do
        [ -n "\$region" ] && echo "\$region" >> angsd_regions.txt
    done < regions.txt

    echo "ANGSD regions:"
    head angsd_regions.txt
    echo "..."

    # --- Build base ANGSD command -------------------------------------------
    ANGSD_CMD="angsd -bam bam.list"
    ANGSD_CMD="\$ANGSD_CMD -ref ${reference}"
    ANGSD_CMD="\$ANGSD_CMD -rf angsd_regions.txt"
    ANGSD_CMD="\$ANGSD_CMD -out ${chunk_id}_raw"
    ANGSD_CMD="\$ANGSD_CMD -nThreads ${task.cpus}"

    # Optional sites file
    if [ -s "${sites_file}" ] && [ "${sites_file}" != "null" ]; then
        echo "Using sites file: ${sites_file}"
        ANGSD_CMD="\$ANGSD_CMD -sites ${sites_file}"
    fi

    # Decide whether we have a usable config JSON
    USE_CONFIG=true
    if [ ! -f "${config}" ] || [ "${config}" = "null" ]; then
        USE_CONFIG=false
    fi

    if \$USE_CONFIG; then
        echo "Loading ANGSD parameters from config file: ${config}"

        python3 << 'PYTHON_CONFIG' > angsd_params.sh
import json, os, sys

config_path = r"${config}"

try:
    if (not config_path) or (config_path == "null") or (not os.path.exists(config_path)):
        raise FileNotFoundError("No config file found")

    with open(config_path) as f:
        config = json.load(f)

    basic   = config.get("basic_options", {}) or {}
    filters = config.get("filters", {}) or {}
    glopts  = config.get("genotype_likelihoods", {}) or {}
    snp     = config.get("snp_calling", {}) or {}
    pooled  = config.get("pooled_sequencing", {}) or {}
    outopt  = config.get("output_options", {}) or {}

    params = []

    # Track which key groups have been set by user
    gl_set          = False
    doGlf_set       = False
    doCounts_set    = False
    doGeno_set      = False
    doPost_set      = False
    doMaf_set       = False
    doMajorMinor_set= False
    dobcf_set       = False
    ignore_rg_set   = False

    def has(d, k):
        return (k in d) and (d[k] is not None)

    # ------------- BASIC OPTIONS -------------
    if has(basic, "doCounts"):
        params.append(f"-doCounts {int(basic['doCounts'])}")
        doCounts_set = True
    if has(basic, "remove_bads"):
        params.append(f"-remove_bads {int(basic['remove_bads'])}")
    if has(basic, "unique_only"):
        params.append(f"-uniqueOnly {int(basic['unique_only'])}")
    if has(basic, "only_proper_pairs"):
        params.append(f"-only_proper_pairs {int(basic['only_proper_pairs'])}")
    if has(basic, "trim"):
        params.append(f"-trim {int(basic['trim'])}")
    if has(basic, "C"):
        params.append(f"-C {basic['C']}")
    if has(basic, "baq"):
        params.append(f"-baq {basic['baq']}")
    if has(basic, "checkBamHeaders"):
        params.append(f"-checkBamHeaders {int(basic['checkBamHeaders'])}")

    # ------------- FILTERS -------------------
    if has(filters, "minQ"):
        params.append(f"-minQ {filters['minQ']}")
    if has(filters, "minMapQ"):
        params.append(f"-minMapQ {filters['minMapQ']}")
    if has(filters, "minInd"):
        params.append(f"-minInd {filters['minInd']}")
    if has(filters, "setMinDepth"):
        params.append(f"-setMinDepth {filters['setMinDepth']}")
    if has(filters, "setMaxDepth"):
        params.append(f"-setMaxDepth {filters['setMaxDepth']}")
    if has(filters, "setMinDepthInd"):
        params.append(f"-setMinDepthInd {filters['setMinDepthInd']}")
    if has(filters, "setMaxDepthInd"):
        params.append(f"-setMaxDepthInd {filters['setMaxDepthInd']}")

    # ------------- GENOTYPE LIKELIHOODS ------
    if has(glopts, "GL"):
        gl_model = glopts["GL"]
        params.append(f"-GL {gl_model}")
        gl_set = True
    if has(glopts, "doGlf"):
        params.append(f"-doGlf {glopts['doGlf']}")
        doGlf_set = True
    if has(glopts, "doGeno"):
        params.append(f"-doGeno {glopts['doGeno']}")
        doGeno_set = True
    if has(glopts, "doPost"):
        params.append(f"-doPost {glopts['doPost']}")
        doPost_set = True
    if has(glopts, "geno_minDepth"):
        params.append(f"-geno_minDepth {glopts['geno_minDepth']}")
    if has(glopts, "geno_minMM"):
        params.append(f"-geno_minMM {glopts['geno_minMM']}")

    # ------------- SNP CALLING ---------------
    if has(snp, "doMajorMinor"):
        params.append(f"-doMajorMinor {snp['doMajorMinor']}")
        doMajorMinor_set = True
    if has(snp, "doMaf"):
        params.append(f"-doMaf {snp['doMaf']}")
        doMaf_set = True
    if has(snp, "SNP_pval"):
        params.append(f"-SNP_pval {snp['SNP_pval']}")
    if has(snp, "rmTriallelic"):
        params.append(f"-rmTriallelic {int(snp['rmTriallelic'])}")
    if has(snp, "skipTriallelic"):
        params.append(f"-skipTriallelic {int(snp['skipTriallelic'])}")
    if has(snp, "minMaf"):
        params.append(f"-minMaf {snp['minMaf']}")
    if has(snp, "minLRT"):
        params.append(f"-minLRT {snp['minLRT']}")

    # ------------- POOLED SEQUENCING ---------
    if has(pooled, "doSaf"):
        params.append(f"-doSaf {pooled['doSaf']}")
    if has(pooled, "fold"):
        params.append(f"-fold {int(pooled['fold'])}")
    if has(pooled, "anc"):
        params.append(f"-anc {pooled['anc']}")

    # ------------- OUTPUT OPTIONS ------------
    if has(outopt, "dumpCounts"):
        params.append(f"-dumpCounts {outopt['dumpCounts']}")
    if has(outopt, "doDepth"):
        params.append(f"-doDepth {outopt['doDepth']}")
    if has(outopt, "maxDepth"):
        params.append(f"-maxDepth {outopt['maxDepth']}")

    # ------------- Minimal Beagle + BCF enforcement -------------
    # We ONLY auto-add things required to ensure Beagle and BCF exist.
    if not gl_set:
        params.append("-GL 1")
    if not doGlf_set:
        params.append("-doGlf 2")            # Beagle output
    if not doCounts_set:
        params.append("-doCounts 1")
        doCounts_set = True
    if not doMajorMinor_set:
        params.append("-doMajorMinor 1")
        doMajorMinor_set = True
    if not doMaf_set:
        params.append("-doMaf 1")
        doMaf_set = True
    if not doPost_set:
        params.append("-doPost 1")
        doPost_set = True
    if not doGeno_set:
        params.append("-doGeno 1")
        doGeno_set = True
    if not dobcf_set:
        params.append("-dobcf 1")            # BCF output
        dobcf_set = True
    # ignore-RG recommended by ANGSD for BCF output
    params.append("--ignore-RG 0")

    # Emit shell code to set ANGSD_PARAMS and some logging
    joined = " ".join(params)
    print(f"ANGSD_PARAMS='{joined}'")
    print(f"echo 'Loaded {len(params)} parameters from config (including minimal Beagle/BCF defaults as needed)'")

except Exception as e:
    # On any error, fall back to minimal defaults that guarantee Beagle & BCF
    sys.stderr.write(f"# Error parsing config: {e}\\n")
    print("echo 'WARNING: Failed to parse config file, using minimal defaults for Beagle/BCF'")
    print("ANGSD_PARAMS='-GL 1 -doGlf 2 -doMajorMinor 1 -doMaf 1 -doCounts 1 -doGeno 1 -doPost 1 -dobcf 1 --ignore-RG 0'")
PYTHON_CONFIG

        # Bring ANGSD_PARAMS into this shell
        source angsd_params.sh

        if [ -n "\$ANGSD_PARAMS" ]; then
            ANGSD_CMD="\$ANGSD_CMD \$ANGSD_PARAMS"
        fi
    else
        echo "No ANGSD config file provided - using ANGSD defaults + minimal Beagle/BCF options"
        ANGSD_CMD="\$ANGSD_CMD -GL 1 -doGlf 2 -doMajorMinor 1 -doMaf 1 -doCounts 1 -doGeno 1 -doPost 1 -dobcf 1 --ignore-RG 0"
    fi

    echo "Full ANGSD command:"
    echo "  \$ANGSD_CMD"
    echo "----------------------------------------"

    # Run ANGSD
    \$ANGSD_CMD
    ANGSD_EXIT=\$?

    echo "----------------------------------------"
    echo "ANGSD exit status: \$ANGSD_EXIT"

    # Quick summary of outputs
    if [ -f "${chunk_id}_raw.beagle.gz" ]; then
        SITES=\$(zcat ${chunk_id}_raw.beagle.gz 2>/dev/null | tail -n +2 | wc -l)
        echo "Sites with genotype likelihoods in Beagle: \$SITES"
    else
        echo "No Beagle file generated"
    fi

    if [ -f "${chunk_id}_raw.bcf" ]; then
        VARIANTS=\$(bcftools view -H ${chunk_id}_raw.bcf 2>/dev/null | wc -l || echo "0")
        echo "Variants in BCF: \$VARIANTS"
    else
        echo "No BCF file generated"
    fi

    # Reheader Beagle sample names using bam.list
    if [ -f "${chunk_id}_raw.beagle.gz" ] && [ -f "bam.list" ]; then
        echo "Reheadering Beagle sample names based on bam.list"

        # Build clean sample names:
        : > samples.txt
        while read -r bam; do
            base=\$(basename "\$bam")
            base=\${base%.bam}
            base=\${base%.filtered}
            echo "\$base" >> samples.txt
        done < bam.list

        # Write final Beagle in one shot
        zcat "${chunk_id}_raw.beagle.gz"  | \
        awk -v names="\$(tr '\n' ' ' < samples.txt)" '
            NR==1 {
                n = split(names, a, " ")
                # Keep only marker allele1 allele2, then each sample three times
                out = \$1 FS \$2 FS \$3
                for (i = 1; i <= n; i++) {
                    out = out FS a[i] FS a[i] FS a[i]
                }
                print out
                next
            }
            { print }
        ' | bgzip > "${chunk_id}.beagle.gz"

    else
        echo "Skipping Beagle reheader (missing ${chunk_id}_raw.beagle.gz or bam.list)"
    fi
    
    # Join pos, counts, hwe, snpStat files
    # Reheader counts sample names and join with positions using bam.list
    if [ -f "${chunk_id}_raw.counts.gz" ] && [ -f "bam.list" ]; then
        echo "Reheadering Beagle sample names based on bam.list"

        # Build clean sample names:
        : > samples.txt
        while read -r bam; do
            base=\$(basename "\$bam")
            base=\${base%.bam}
            base=\${base%.filtered}
            echo "\$base" >> samples.txt
        done < bam.list

        if [ -f "${chunk_id}_raw.hwe.gz" ] && [ -f "${chunk_id}_raw.snpStat.gz" ]; then
            paste \
            <(zcat "${chunk_id}_raw.pos.gz") \
            <(
                zcat "${chunk_id}_raw.counts.gz"  | \
                awk -v names="\$(sed 's/\$/_depth/' samples.txt | tr '\n' ' ')" '
                    BEGIN { FS = OFS = "\t" }
                    NR==1 {
                        n = split(names, a, " ")
                        # rename header to be the individual names
                        out = a[1]
                        for (i = 2; i <= n; i++) {
                            out = out OFS a[i]
                        }
                        print out
                        next
                    }
                    { print }
                ' \
                | sed 's/\t*\$//'
            ) \
            <(zcat "${chunk_id}_raw.hwe.gz" | cut -f3-) \
            <(zcat "${chunk_id}_raw.snpStat.gz" | cut -f3- | sed 's/\t*\$//') \
            | bgzip > "${chunk_id}.counts.gz"
        else
           paste \
            <(zcat "${chunk_id}_raw.pos.gz") \
            <(
                zcat "${chunk_id}_raw.counts.gz"  | \
                awk -v names="\$(sed 's/\$/_depth/' samples.txt | tr '\n' ' ')" '
                    BEGIN { FS = OFS = "\t" }
                    NR==1 {
                        n = split(names, a, " ")
                        # rename header to be the individual names
                        out = a[1]
                        for (i = 2; i <= n; i++) {
                            out = out OFS a[i]
                        }
                        print out
                        next
                    }
                    { print }
                ' \
                | sed 's/\t*\$//'
            ) \
            | bgzip > "${chunk_id}.counts.gz"
        fi


    else
        echo "Skipping Count reheader (missing ${chunk_id}_raw.counts.gz or bam.list)"
    fi

    # Beagle is already handled separately (reheadered from RAW_PREFIX.beagle.gz)
    # Now map the rest of the ANGSD outputs to the clean names
    for ext in mafs.gz bcf geno.gz depthSample depthGlobal arg; do
        src="${chunk_id}_raw.\${ext}"
        dest="${chunk_id}.\${ext}"
        if [ -f "\$src" ]; then
            echo "Creating final \$dest from \$src"
            cp "\$src" "\$dest"
        else
            echo "Note: \$src not found, skipping"
        fi
    done

    if [ \$ANGSD_EXIT -ne 0 ]; then
        echo "ERROR: ANGSD failed for chunk ${chunk_id}"
        exit 1
    fi

    echo "Chunk ${chunk_id} processing complete"
    """
}
