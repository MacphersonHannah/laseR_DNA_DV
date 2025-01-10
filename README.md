<h1>Long-Read Structural Variant (SV) Calling Workflow</h1>  

This repository contains a Snakemake-based pipeline for analyzing long-read sequencing data. The workflow supports structural variant (SV) calling, phasing, and haplotype-specific genome reference construction using tools such as Minimap2, DeepVariant, Longshot, Samtools, Whatshap, and others. Some of these tools do the same things through different methods, making some rules redundant. This is purposeful as it allows comparison between different tools.

<h2>Getting Started</h2>
<h3>Prerequisites</h3>

Install Required Tools:

- Miniconda/Conda: To manage environments.
- Snakemake: Workflow manager.
- Singularity: For containerized tools like DeepVariant.

<h3>Clone the Repository:</h3>

```bash
git clone https://github.com/MacphersonHannah/laseR_DNA_DV.git
cd laseR_DNA_DV
```

<h3>Environment Setup:</h3>

```bash
conda env create -f envs/environment.yml
conda activate laseR_DNA_DV
```

<h3>Input Configuration</h3>
Modify the config.yaml file to specify the following parameters:

| Parameter             | Description                            | Example                             |
|-----------------------|----------------------------------------|-------------------------------------|
| `workdir`             | Directory for results                 | `/path/to/output/`                 |
| `sample_name`         | Prefix for output files               | `sample1`                          |
| `genome`              | Reference genome in FASTA format      | `/path/to/genome.fasta`            |
| `fastq`               | Path to input FASTQ file              | `/path/to/reads.fastq`             |
| `threads`             | Number of CPU threads to use          | `50`                               |


Example config.yaml:

```bash
workdir: "/path/to/output"
sample_name: "sample1"
genome: "/path/to/genome.fa"
fastq: "/path/to/reads.fastq" (multiple fastq files can be concatenated into one file prior to pipeline usage)
region: "chr4:10078847-140078843"
threads: 50
deepvariant_bin_version: "1.5.0"
num_shards: 30
```


<h3>Running the Workflow</h3>

Run the pipeline using the following command:

```bash
snakemake --use-conda -j 50 all --verbose 2>&1 | tee snakemake.log
```

Where:

```--use-conda``` Ensures Conda environments are used for dependency management.
</br>```-j 50``` Specifies 50 parallel jobs (adjust based on available resources).
</br>```all``` Targets the final rule to execute the entire workflow. To run just one rule and its prerequisites, change 'all' to any other rule name.
</br>```--verbose``` Enables detailed logging.
</br>```snakemake.log``` Log file will automatically be created in snakemake folder, but this can be changed to be created in any location.

<h3>Workflow Overview</h3>

The pipeline comprises the following steps:

<b>Read Preparation:</b>

Truncates read names for compatibility (truncate_read_names).
Aligns reads to the reference genome using Minimap2 (align).

<b>Variant Calling:</b>

Calls SNPs using Longshot (longshot).
Performs variant calling using DeepVariant (run_deepvariant).

<b>Filtering and Phasing:</b>

Filters BAM and VCF files for the region of interest (filter_bam, filter_vcf_for_region).
Phases variants using Whatshap (phase_vcf).

<b>Haplotype Reference Construction:</b>
Creates haplotype-specific genome references from phased variants (create_haplotyped_personalised_reference).

<b>Future Updates:</b>
Detect structural variants (SVs) using Sniffles or copy number variations (CNVs) using NanoVar (currently commented out for customization).

<b>Output</b>
The pipeline generates the following outputs:

| Output Type           | Description                                | Example Path                                     |
|-----------------------|--------------------------------------------|-------------------------------------------------|
| Aligned BAM           | Reads aligned to the reference genome.     | `mapping/sample1_sorted.bam`                   |
| Longshot VCF          | SNP calls from Longshot.                   | `longshot/sample1.vcf.gz`                      |
| DeepVariant VCF       | Variant calls from DeepVariant.            | `deepvariant/sample1.vcf`                      |
| Phased VCFs           | Phased variants for haplotypes.            | `phased_vcf/sample1_phased_region_longshot.vcf.gz` |
| Haplotype References  | Haplotype-specific genome references.      | `personalised_genomes/sample1_hap1_longshot_sorted.bam` |


<h3>Troubleshooting</h3>

<b>Missing Dependencies:</b>

Ensure all tools are installed in the Conda environment. Recreate the environment if necessary:
```bash
conda env create -f envs/environment.yml
```

<b>Debugging Failed Jobs:</b>

Check the snakemake.log file or rule-specific logs for errors.
Re-run individual rules for debugging:

```bash
snakemake --use-conda -j 1 <rule_name>
```

<b>Singularity Issues:</b>

Ensure Singularity is installed and configured if using containerized tools like DeepVariant.
