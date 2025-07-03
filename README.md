<h1>laseR_DNA_DV: Long-Read Variant Calling Pipeline</h1>

<p><strong>laseR_DNA_DV</strong> is a modular Snakemake pipeline for analyzing long-read DNA sequencing data, with a focus on accurate variant calling using <strong>DeepVariant</strong>. It supports alignment, variant calling via both DeepVariant and Longshot, and optional region-based filtering.</p>

<p>Haplotype reference generation and phasing steps are included in the codebase but currently commented out for simplicity. This makes the current version faster and more lightweight.</p>

<hr>

<h2>🧰 Tools Used</h2>
<ul>
  <li><strong>Minimap2</strong> — long-read aligner</li>
  <li><strong>DeepVariant</strong> — deep learning-based variant caller</li>
  <li><strong>Longshot</strong> — haplotype-aware SNP caller</li>
  <li><strong>Samtools</strong> — BAM manipulation</li>
  <li><strong>Bcftools</strong> — optional VCF filtering</li>
</ul>

<hr>

<h2>🚀 Getting Started</h2>

<h3>1. Prerequisites</h3>
<ul>
  <li>Miniconda or Conda</li>
  <li>Snakemake (tested on ≥7.18)</li>
  <li>Singularity (for DeepVariant)</li>
</ul>

<h3>2. Clone the Repository</h3>
<pre><code class="language-bash">
git clone https://github.com/MacphersonHannah/laseR_DNA_DV.git
cd laseR_DNA_DV
</code></pre>

<h3>3. Set Up the Environment</h3>
<pre><code class="language-bash">
conda env create -f envs/environment.yml
conda activate laseR_DNA_DV
</code></pre>

<hr>

<h2>⚙️ Input Configuration</h2>

<p>Edit the <code>config.yaml</code> file with your sample-specific details:</p>

<table>
  <thead>
    <tr><th>Parameter</th><th>Description</th><th>Example</th></tr>
  </thead>
  <tbody>
    <tr><td><code>workdir</code></td><td>Output directory</td><td><code>/path/to/output/</code></td></tr>
    <tr><td><code>sample_name</code></td><td>Prefix for output files</td><td><code>sample1</code></td></tr>
    <tr><td><code>genome</code></td><td>Reference genome (FASTA)</td><td><code>/path/to/genome.fa</code></td></tr>
    <tr><td><code>fastq</code></td><td>FASTQ file (already demultiplexed)</td><td><code>/path/to/reads.fastq</code></td></tr>
    <tr><td><code>region</code></td><td>Optional region for targeted analysis</td><td><code>chr4:10000000-14000000</code></td></tr>
    <tr><td><code>use_region</code></td><td>Enable region-based filtering</td><td><code>true</code> or <code>false</code></td></tr>
    <tr><td><code>threads</code></td><td>Number of threads</td><td><code>50</code></td></tr>
    <tr><td><code>deepvariant_bin_version</code></td><td>DeepVariant container version</td><td><code>1.5.0</code></td></tr>
    <tr><td><code>num_shards</code></td><td>DeepVariant parallel threads</td><td><code>30</code></td></tr>
    <tr><td><code>root_workdir</code></td><td>Path mount for Singularity</td><td><code>/absolute/path/to/project</code></td></tr>
  </tbody>
</table>

<p><strong>Note:</strong> If using barcoded reads, you must demultiplex them before running this workflow. Only one sample per FASTQ file is expected.</p>

<hr>

<h2>🏃 Running the Pipeline</h2>

<p>Run the full workflow:</p>

<pre><code class="language-bash">
snakemake --use-conda -j 50 all --verbose 2>&1 | tee snakemake.log
</code></pre>

<p>To dry run without executing:</p>

<pre><code class="language-bash">
snakemake --use-conda -j 50 all --dry-run --verbose 2>&1 | tee snakemake_dryrun.log
</code></pre>

<hr>

<h2>📋 Workflow Overview</h2>

<ul>
  <li><strong>truncate_read_names</strong> — Simplifies FASTQ headers</li>
  <li><strong>align</strong> — Aligns reads using Minimap2</li>
  <li><strong>run_deepvariant</strong> — Runs DeepVariant in a Singularity container</li>
  <li><strong>prepare_bam</strong> — Filters BAM file by region (if <code>use_region: true</code>)</li>
  <li><strong>longshot</strong> — Optional SNP calling with Longshot</li>
  <li><strong>filter_vcf_for_region</strong> — Optional VCF trimming by region</li>
</ul>

<p>Phasing and haplotype-reference rules are commented out for simplicity and speed. They can be re-enabled if needed in the Snakefile.</p>

<hr>

<h2>🧾 Output</h2>

<table>
  <thead>
    <tr><th>Output Type</th><th>Description</th><th>Example Path</th></tr>
  </thead>
  <tbody>
    <tr><td>Aligned BAM</td><td>Reads aligned to the reference</td><td><code>mapping/sample1_sorted.bam</code></td></tr>
    <tr><td>DeepVariant VCF</td><td>Variant calls</td><td><code>deepvariant/sample1.vcf</code></td></tr>
    <tr><td>DeepVariant GVCF</td><td>GVCF for genotyping</td><td><code>deepvariant/sample1.gvcf</code></td></tr>
    <tr><td>Longshot VCF</td><td>Optional SNP calls</td><td><code>longshot/sample1.vcf.gz</code></td></tr>
    <tr><td>Filtered VCF</td><td>If region filtering is enabled</td><td><code>longshot/sample1_region.vcf.gz</code></td></tr>
  </tbody>
</table>

<hr>

<h2>🐛 Troubleshooting</h2>

<h3>Missing Dependencies</h3>
<p>Recreate the Conda environment if something is missing:</p>
<pre><code class="language-bash">
conda env create -f envs/environment.yml
</code></pre>

<h3>Debugging a Rule</h3>
<pre><code class="language-bash">
snakemake --use-conda -j 1 <rule_name>
</code></pre>

<h3>DeepVariant & Singularity</h3>
<p>Make sure Singularity is correctly installed and the <code>root_workdir</code> path is bindable.</p>

---
