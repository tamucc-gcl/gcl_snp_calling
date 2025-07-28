# gcl_snp_calling
![](pipeline_dag.png)

This pipeline takes a directory of mapped and indexed BAM files, filters them, and performs SNP calling using [FreeBayes](https://github.com/freebayes/freebayes). It is written in [Nextflow](https://www.nextflow.io/) DSL2 and is designed to be run either directly or via SLURM using the provided sbatch script.

---

## 🔧 Tools Used

| Tool | Purpose | Link |
|------|---------|------|
| [samtools](http://www.htslib.org/) | BAM file filtering and indexing | https://www.htslib.org/ |
| [freebayes](https://github.com/freebayes/freebayes) | Variant calling | https://github.com/freebayes/freebayes |
| [bcftools](http://www.htslib.org/doc/bcftools.html) | VCF filtering and manipulation | https://samtools.github.io/bcftools/ |

---

## Folder Structure
```
📁 project-root/
├── 📁 gcl_snp_calling/             # ✅ Must be present (cloned repo)
│   ├── 📄 main.nf                  # Entry point for the pipeline
│   ├── 📄 nextflow.config          # Execution and profile config
│   ├── 📄 run_snpCall.sbatch            # SLURM submission script
│   └── 📁 modules/                 # Individual module processes
│       ├── combine_vcfs.nf
│       ├── create_chunks.nf
│       ├── filter_bam.nf
│       ├── freebayes_chunk.nf
│       ├── summarize_vcf.nf
├── 📁 data/                          
│   ├── 📁 bam/                       # ✅ Must be present (unfiltered BAM files - can be output by [QC Pipeline](https://github.com/tamucc-gcl/gcl_illumina_qc))
│   ├── 📁 snps/                      # 🚀 Created by pipeline
├── 📁 genome/                        # ✅ Must be present Downloaded and indexed reference (can be created by [QC Pipeline](https://github.com/tamucc-gcl/gcl_illumina_qc))
├── 📁 logs/                          # SLURM/Nextflow job logs
```

## Parameter Setting Files
- Example: [BAM Filter Settings](config_files/bam_filter_parameters.json)
- Example: [freebayes Settings](config_files/freebayes_parameters.json)

## To Run

1. Clone repo into your directory.
2. Run the pipeline from the project root:

```bash
nextflow run gcl_snp_calling/main.nf \
    -profile slurm \
    -resume \
    --bams "./data/bam/*.bam" \
    --reference "./genome/genome.fa" \
    --freebayes_config {freebayes_parameters} \
    --bam_filter_config {bam_filter_parameters} \
    --num_chunks {nChunk} \
    --output_dir "data/snps" \
    --output_vcf "all_samples.vcf.gz"
```
3. Or Run in SLURM `sbatch run_snpCall.sbatch {freebayes_parameters.json} {bam_filter_parameters} 25`

Edit `gcl_snp_calling/run_snpCall.sbatch` to customize resources used by orchestra conductor job (e.g., CPUs, memory, partition).

Edit `gcl_illumina_qc/nextflow.config` to customize resources used by each QC stage (e.g., CPUs, memory, partition).