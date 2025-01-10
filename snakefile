# conda: long-read_SV_calling
# snakemake --use-conda -j 50 all --verbose 2>&1 | tee snakemake.log


from snakemake.utils import min_version, validate
from pathlib import Path
from os import path
import os
import yaml
import shutil
 
min_version("7.18")

configfile: "config.yaml"
shutil.copy2(SNAKEDIR + "/config.yaml", WORKDIR)
workdir: config["workdir"]


WORKDIR = config["workdir"]
SNAKEDIR = path.dirname(workflow.snakefile)

sample = config["sample_name"]
region = config["region"]

target_list = [
    "longshot/" + sample + ".vcf.gz",
    "deepvariant/" + sample + ".vcf",
    "deepvariant/" + sample + ".gvcf",
    "personalised_genomes/" + sample + "_hap1_longshot_sorted.bam",
    "personalised_genomes/" + sample + "_hap2_longshot_sorted.bam",
    "personalised_genomes/" + sample + "_hap1_deepvariant_sorted.bam",
    "personalised_genomes/" + sample + "_hap2_deepvariant_sorted.bam",
    "deepvariant/" + sample + ".vcf"
]



# Align reads -------------------------------------------------------------
rule truncate_read_names:
    input:
        fq_concat = config['fastq']
    output:
        fq_truncated = path.join("processed_reads", f"{sample}_truncated_reads.fastq")
    threads: config["threads"]
    conda:
        "envs/environment.yml"
    shell: """
        awk 'BEGIN {{FS=" "}} /^@/{{print $1; next}}1' {input.fq_concat} > {output.fq_truncated}
    """

rule align:
    input:
        fq = rules.truncate_read_names.output.fq_truncated,
        ref = config["genome"]
    output:
        bam = path.join("mapping", f"{sample}_sorted.bam"),
        bai = path.join("mapping", f"{sample}_sorted.bam.bai")
    threads: config["threads"]
    params:
        minimap2_opts = config["minimap2_opts"],
        sample_name = sample
    conda:
        "envs/environment.yml"
    shell:
        """
        # Align reads
        minimap2 {params.minimap2_opts} -y -x map-ont -t {threads} -a --eqx -k 17 -K 5g {input.ref} {input.fq} > temp_unsorted.sam

        # Add read group information
        samtools addreplacerg -@ {threads} -r 'ID:{params.sample_name} SM:{params.sample_name}' -o temp_unsorted_rg.bam temp_unsorted.sam

        # Sort and index the BAM file
        samtools sort -@ {threads} -O BAM -o {output.bam} temp_unsorted_rg.bam
        samtools index {output.bam}

        # Clean up intermediate files
        rm temp_unsorted.sam temp_unsorted_rg.bam
        """




rule run_deepvariant:
    input:
        bam_file = rules.align.output.bam
    output:
        vcf = path.join("deepvariant", f"{sample}.vcf"),
        gvcf = path.join("deepvariant", f"{sample}.gvcf"),
    params:
        genome = config["genome"],
        region = config["region"],
        deepvariant_bin_version = config["deepvariant_bin_version"],
        num_shards = config["num_shards"],
        root_workdir = config["root_workdir"],
    log:
        "deepvariant/run_deepvariant.log"
    shell:
        """

        singularity run -B /usr/lib/locale/:/usr/lib/locale/,"{params.root_workdir}":"{params.root_workdir}" \
        docker://google/deepvariant:"{params.deepvariant_bin_version}" \
        /opt/deepvariant/bin/run_deepvariant \
        --model_type=ONT_R104 \
        --ref="{params.genome}" \
        --reads="{input.bam_file}" \
        --regions "{params.region}" \
        --output_vcf="{output.vcf}" \
        --output_gvcf="{output.gvcf}" \
        --num_shards="{params.num_shards}" 2>> {log}
        
        """



# Filter BAM for region of interest ---------------------------------------
rule filter_bam:
    input:
        bam = rules.align.output.bam,
        bai = rules.align.output.bai
    output:
        filtered_bam = path.join("mapping", f"{sample}_region.bam")
    threads: config["threads"]
    conda:
        "envs/environment.yml"
    shell:
        """
        samtools view -b {input.bam} {config[region]} > {output.filtered_bam}
        samtools index {output.filtered_bam}
        """


# Sort and index the filtered BAM -----------------------------------------
rule sort_and_index:
    input:
        bam = rules.filter_bam.output.filtered_bam
    
    output:
        bam = path.join("mapping", f"{sample}.bam"),
        bai = path.join("mapping", f"{sample}.bam.bai")
    
    threads: config["threads"]
    conda:
        "envs/environment.yml"
    
    shell:
        """
        samtools sort -@ {threads} -O BAM -o {output.bam} {input.bam}
        samtools index {output.bam}
        """


# Call SNPs ---------------------------------------------------------------
rule longshot:
    input:
        bam = rules.sort_and_index.output.bam
    output:
        vcf = path.join("longshot", f"{sample}.vcf.gz")
    threads: config["threads"]
    conda:
        "envs/environment.yml"
    shell: """
        set -euo pipefail
        longshot --bam {input.bam} --ref {config[genome]} --out {output.vcf}.temp.vcf
        bgzip -c {output.vcf}.temp.vcf > {output.vcf}
        tabix -p vcf {output.vcf}
        rm {output.vcf}.temp.vcf
    """


# Filter VCF for region of interest ---------------------------------------
rule filter_vcf_for_region:
    input:
        vcf = rules.longshot.output.vcf
    output:
        filtered_vcf = path.join("longshot", f"{sample}_region.vcf.gz")
    threads: config["threads"]
    conda:
        "envs/bcftools.yml"
    shell:
        """
        bcftools view --regions {config[region]} {input.vcf} -Oz -o {output.filtered_vcf}
        tabix -p vcf {output.filtered_vcf}
        """


# Phase VCF ---------------------------------------------------------------
rule phase_vcf:
    input:
        bam = rules.sort_and_index.output.bam,
        longshot_vcf = rules.filter_vcf_for_region.output.filtered_vcf,
        deepvariant_vcf = rules.run_deepvariant.output.vcf,
        ref = config["genome"]
    output:
        longshot_phased_vcf = path.join("phased_vcf", f"{sample}_phased_region_longshot.vcf.gz"),
        deepvariant_phased_vcf = path.join("phased_vcf", f"{sample}_phased_region_deepvariant.vcf.gz")
    threads: config["threads"]
    conda:
        "envs/whatshap.yml"
    shell:
        """
        whatshap phase --reference {input.ref} --output {output.longshot_phased_vcf} {input.longshot_vcf} {input.bam}
        tabix -p vcf {output.longshot_phased_vcf}

        whatshap phase --reference {input.ref} --output {output.deepvariant_phased_vcf} {input.deepvariant_vcf} {input.bam}
        tabix -p vcf {output.deepvariant_phased_vcf}
        """


# Create haplotyped personalised reference --------------------------------
rule create_haplotyped_personalised_reference:
    input:
        ref = config["genome"],
        longshot_vcf = rules.phase_vcf.output.longshot_phased_vcf,
        deepvariant_vcf = rules.phase_vcf.output.deepvariant_phased_vcf
    output:
        hap1_ref_longshot = path.join("personalised_genomes", f"{sample}_longshot_hap1_region.fa.gz"),
        hap2_ref_longshot = path.join("personalised_genomes", f"{sample}_longshot_hap2_region.fa.gz"),
        hap1_ref_deepvariant = path.join("personalised_genomes", f"{sample}_deepvariant_hap1_region.fa.gz"),
        hap2_ref_deepvariant = path.join("personalised_genomes", f"{sample}_deepvariant_hap2_region.fa.gz")
    conda:
        "envs/bcftools.yml"
    shell:
        """
        # Extract the region of interest from the reference
        samtools faidx {input.ref} {config[region]} > temp_region.fa

        # Generate consensus sequences using Longshot VCF
        bcftools consensus -H 1 -f temp_region.fa {input.longshot_vcf} | bgzip -c > {output.hap1_ref_longshot}
        bcftools consensus -H 2 -f temp_region.fa {input.longshot_vcf} | bgzip -c > {output.hap2_ref_longshot}

        # Generate consensus sequences using deepvariant VCF
        bcftools consensus -H 1 -f temp_region.fa {input.deepvariant_vcf} | bgzip -c > {output.hap1_ref_deepvariant}
        bcftools consensus -H 2 -f temp_region.fa {input.deepvariant_vcf} | bgzip -c > {output.hap2_ref_deepvariant}

        # Clean up
        rm temp_region.fa
        """




rule align_haplotypes_to_hg38:
    input:
        # From the previous rule
        hap1_ref_longshot = path.join("personalised_genomes", f"{sample}_longshot_hap1_region.fa.gz"),
        hap2_ref_longshot = path.join("personalised_genomes", f"{sample}_longshot_hap2_region.fa.gz"),
        hap1_ref_deepvariant = path.join("personalised_genomes", f"{sample}_deepvariant_hap1_region.fa.gz"),
        hap2_ref_deepvariant = path.join("personalised_genomes", f"{sample}_deepvariant_hap2_region.fa.gz"),
        ref = config["genome"]
    output:
        # Outputs for Longshot
        hap1_bam_longshot = "personalised_genomes/" + sample + "_hap1_longshot_sorted.bam",
        hap1_bai_longshot = "personalised_genomes/" + sample + "_hap1_longshot_sorted.bam.bai",
        hap2_bam_longshot = "personalised_genomes/" + sample + "_hap2_longshot_sorted.bam",
        hap2_bai_longshot = "personalised_genomes/" + sample + "_hap2_longshot_sorted.bam.bai",
        # Outputs for Deepvariant
        hap1_bam_deepvariant = "personalised_genomes/" + sample + "_hap1_deepvariant_sorted.bam",
        hap1_bai_deepvariant = "personalised_genomes/" + sample + "_hap1_deepvariant_sorted.bam.bai",
        hap2_bam_deepvariant = "personalised_genomes/" + sample + "_hap2_deepvariant_sorted.bam",
        hap2_bai_deepvariant = "personalised_genomes/" + sample + "_hap2_deepvariant_sorted.bam.bai"
    threads: config["threads"]
    conda:
        "envs/environment.yml"
    shell:
        """
        # Align haplotypes from Longshot
        minimap2 -a -x asm5 {input.ref} {input.hap1_ref_longshot} | samtools sort -@ {threads} -o {output.hap1_bam_longshot}
        samtools index {output.hap1_bam_longshot}
        
        minimap2 -a -x asm5 {input.ref} {input.hap2_ref_longshot} | samtools sort -@ {threads} -o {output.hap2_bam_longshot}
        samtools index {output.hap2_bam_longshot}

        # Align haplotypes from Deepvariant
        minimap2 -a -x asm5 {input.ref} {input.hap1_ref_deepvariant} | samtools sort -@ {threads} -o {output.hap1_bam_deepvariant}
        samtools index {output.hap1_bam_deepvariant}
        
        minimap2 -a -x asm5 {input.ref} {input.hap2_ref_deepvariant} | samtools sort -@ {threads} -o {output.hap2_bam_deepvariant}
        samtools index {output.hap2_bam_deepvariant}
        """






# rule detect_structural_variations:
#     input:
#         bam = rules.align.output.bam
#     output:
#         vcf = path.join("SV_calls", f"{sample}_SV.vcf")
#     threads: config["threads"]
#     conda:
#         "envs/sv_cnv_detection.yml"
#     shell:
#         """
#         sniffles -i {input.bam} -v {output.vcf} --threads {threads}
#         """

# rule detect_copy_number_variation:
#     input:
#         bam = rules.align.output.bam,
#         ref = config["genome"]
#     output:
#         cnv_vcf = path.join("CNV_calls", f"{sample}_CNV.vcf")
#     threads: config["threads"]
#     conda:
#         "envs/sv_cnv_detection.yml"
#     shell:
#         """
#         nanovar -i {input.bam} -r {input.ref} -o {output.cnv_vcf}
#         """




rule all:
    input:
        target_list



#whatshap phase --reference /home/MRHannahMacpherson/MinaRyten/Emil/references/GRCh38.primary_assembly.genome.fa --output test_phased.vcf.gz /home/MRHannahMacpherson/MinaRyten/hannah/240731_AS_SNCA_MAPT_JE_CELLLINES/A53T_828_21.12_SW_Barcode02/longshot/A53T_828_21.12_SW_Barcode02.vcf.gz /home/MRHannahMacpherson/MinaRyten/hannah/240731_AS_SNCA_MAPT_JE_CELLLINES/A53T_828_21.12_SW_Barcode02/mapping/A53T_828_21.12_SW_Barcode02_region.bam
