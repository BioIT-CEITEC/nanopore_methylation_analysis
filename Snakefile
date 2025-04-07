from pathlib import Path
import pandas as pd

configfile: "config.json"
GLOBAL_REF_PATH = config["globalResources"] 

##### BioRoot utilities - reference #####
module BR:
    snakefile: github("BioIT-CEITEC/bioroots_utilities", path="bioroots_utilities.smk",branch="master")
    config: config

use rule * from BR as other_*
config = BR.load_organism()
sample_tab = BR.load_sample()

print(sample_tab.sample_name)

##### Config processing #####
# Folders
#
#sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")
#reference_path = os.path.join(GLOBAL_REF_PATH,config["organism"], config["reference"], "seq", config["reference"] + ".fa")

METHYLATION = {} # empty dictionary for methylations to work with

if config["6mA_methylation"]:
    METHYLATION["6mA"] = "a"
if config["5mC_methylation"]:
    METHYLATION["5mC"] = "m"
    METHYLATION["5hmC"] = "h"

##### Target rules #####
rule all:
    input:
        expand("methylation/{sample_name}/{sample_name}_filtered_{mod_name}.bed", sample_name = sample_tab.sample_name, mod_name = METHYLATION.keys()),
        expand("methylation/{sample_name}/{sample_name}_6mA_100kb.bed", sample_name = sample_tab.sample_name),
        expand("methylation/{sample_name}/{sample_name}_genes-stats.tsv", sample_name = sample_tab.sample_name)


rule create_100kb_windows:
    input:
        genome = config["organism_fasta"]
    output:
        chr_sizes = "methylation/chrom.sizes",
        chr_windows = "methylation/chr_windows.bed"
    conda: 
        "envs/methylation_change.yaml"
    shell:
        """
        cut -f1,2 {input.genome}.fai > {output.chr_sizes}
        bedtools makewindows -g {output.chr_sizes} -w 100000 > {output.chr_windows}
        """

rule filter_modifications:
    input:      
        bed = "methylation/{sample_name}/{sample_name}_modkit.bed"
    output: 
        filtered = "methylation/{sample_name}/{sample_name}_filtered_{mod_name}.bed"
    params: 
        mod_char = lambda wildcards: METHYLATION[wildcards.mod_name]cd 
    shell:
        """
        awk '$4 == "{params.mod_char}" && $11 > 5 && $12 >= 2' {input.bed} >  {output.filtered}
        """

rule create_bedgraph:
    input:      
        filtered = "methylation/{sample_name}/{sample_name}_filtered_{mod_name}.bed"
    output: 
        bed_graph = "methylation/{sample_name}/{sample_name}_filtered_{mod_name}.bedgraph"
    shell:
        """
        awk 'BEGIN {{OFS="\t"}} {{print $1, $2, $3, $11}}' {input.filtered} >  {output.bed_graph}
        """

rule compute_methylation_in_100kbSwindows:
    input:      
        bed = "methylation/{sample_name}/{sample_name}_filtered_{mod_name}.bed",
        chr_windows = "methylation/chr_windows.bed"
    output:
        bed_100kb = "methylation/{sample_name}/{sample_name}_{mod_name}_100kb.bed"
    conda: 
        "envs/methylation_change.yaml"
    shell:
        """
        bedtools map -a {input.chr_windows} -b  {input.bed} -c 11 -o mean > {output.bed_100kb}
        """

rule create_genes_region_bed:
    input: 
        ref_gtf = config["organism_gtf"]
    output:
       "methylation/genes.bed"
    shell:
        """ 
        awk '$3 == "gene"' {input.ref_gtf} | \
        awk 'BEGIN{{OFS="\t"}} {{split($9,a,";"); print $1, $4-1, $5, a[1], ".", $7}}' > {output}
        """

rule modkit_stats:
    input:
        gene_region = "methylation/genes.bed",
        bed = "methylation/{sample_name}/{sample_name}_modkit.bed"
    output:
        tsv = "methylation/{sample_name}/{sample_name}_genes-stats.tsv"
    conda: 
        "envs/methylation.yaml"
    shell:
        """
        bgzip {input.bed}
        tabix -p bed {input.bed}.gz
        modkit stats {input.bed} --min-coverage 5 --regions {input.gene_region} --out-table {output.tsv}
        """

# ### separate bams for all chromosomes in reference
# rule create_bam_for_chromosome:
#     input: 
#         bam = 'aligned/{sample_name}/{sample_name}_sorted.bam'
#     output:
#         bam = 'aligned/{sample_name}/{sample_name}_sorted.REF_chr1'
#     conda: 
#         "envs/methylation_change.yaml"
#     shell:
#         """
#         bamtools split -in {input.bam} -reference
#         """

# rule find_highly_modified_motifs:
#     input: 
#         bed = "methylation/{sample_name}/{sample_name}_modkit.bed"
#     output: 
#         log = "methylation/{sample_name}/{sample_name}_modkit_find_motifs_log.txt",
#         tsv = "methylation/{sample_name}/{sample_name}_motifs.tsv"
#     params: 
#         genome = config["organism_fasta"],
#     conda: 
#         "envs/methylation.yaml"
#     shell: 
#         """
#         modkit motif search -i {input.bed} -r {params.genome} -o {output.tsv} --threads 32 --log {output.log}
#         """


# ### separate bams for modifications
# rule extract_separate_modifications:
#     input: 
#         bam = 'aligned/{sample_name}/{sample_name}_sorted.REF_{chromosome}.bam'
#     output: 
#         filtered_bam = 'aligned/{sample_name}/{modification}/{sample_name}_sorted.REF_{chromosome}_{modification}.bam'
#     params: 
#         mod5mC_5hmC= lambda wildcards, input: int(config["5mC_methylation"]), # convert to 0/1 for bash
#         mod6mA= lambda wildcards, input: int(config["6mA_methylation"])
#     conda: 
#         "envs/methylation_change.yaml"
#     shell:
#         """
#         for file in *_chr*.bam; do
#             output="5mC/${file%.bam}_5mC.bam"
#             samtools view -h "$file" | grep -E "MM:Z:[^;]*5mC" | samtools view -bS > "$output"
#             echo "Zpracován: $file -> $output"
#         done
#         if [ {params.mod5mC_5hmC} -eq 1 ]; then
#             samtools view -h {input.bam} | grep -E "MM:Z:[^;]*5mC" | samtools view -bS > filtered_6mA.bam
#         """

# #rule QC_after_SV_calling: