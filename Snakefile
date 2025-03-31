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

##### Config processing #####
# Folders
#
#sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")
#reference_path = os.path.join(GLOBAL_REF_PATH,config["organism"], config["reference"], "seq", config["reference"] + ".fa")

METHYLATION = []

if config["6mA_methylation"]:
    METHYLATION.append("6mA")
if config["5mC_methylation"]:
    METHYLATION.extend(["5mC", "5hmC"])

##### Target rules #####
rule all:
    input:
        expand("methylation/{sample_name}/{sample_name}_motifs.tsv", sample_name = sample_tab.sample_name), 
        expand('aligned/{sample_name}/{sample_name}_sorted.REF_chr1', sample_name = sample_tab.sample_name)

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
        bedtools makewindows -g {output} -w 100000 > {output.chr_windows}
        """

rule filter_modifications:
    input:      
        bed = "methylation/{sample_name}/{sample_name}_modkit.bed"
    output: 
        mod6mA = "methylation/{sample_name}/{sample_name}_filtered_6mA.bed", 
        mod5hmC = "methylation/{sample_name}/{sample_name}_filtered_5hmC.bed", 
        mod5mC = "methylation/{sample_name}/{sample_name}_filtered_5mC.bed"
    shell:
        """
        awk '$4 == "a" && $11 > 5 && $12 >= 2' {input.bed} >  {output.mod6mA}
        awk '$4 == "a" && $11 > 5 && $12 >= 2' {input.bed} >  {output.mod5hmC}
        awk '$4 == "a" && $11 > 5 && $12 >= 2' {input.bed} >  {output.mod5mC}
        """

rule filter_modifications:
    input:      
        mod6mA = "methylation/{sample_name}/{sample_name}_filtered_6mA.bed", 
        mod5hmC = "methylation/{sample_name}/{sample_name}_filtered_5hmC.bed", 
        mod5mC = "methylation/{sample_name}/{sample_name}_filtered_5mC.bed"
    output: 
        mod6mA = "methylation/{sample_name}/{sample_name}_filtered_6mA.bedgraph", 
        mod5hmC = "methylation/{sample_name}/{sample_name}_filtered_5hmC.bedgraph", 
        mod5mC = "methylation/{sample_name}/{sample_name}_filtered_5mC.bedgraph"
    shell:
        """
        awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $11}' {input.mod6mA} >  {output.mod6mA}
        awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $11}' {input.mod5hmC} >  {output.mod5hmC}
        awk 'BEGIN {OFS="\t"} {print $1, $2, $3, $11}' {input.mod5mC} >  {output.mod5mC}
        """

rule compute_methylation_in_100kbSwindows:
    input:      
        mod6mA = "methylation/{sample_name}/{sample_name}_filtered_6mA.bed", 
        mod5hmC = "methylation/{sample_name}/{sample_name}_filtered_5hmC.bed", 
        mod5mC = "methylation/{sample_name}/{sample_name}_filtered_5mC.bed",
        chr_windows = "methylation/chr_windows.bed"
    output:
        mod6mA_100kb = "methylation/{sample_name}/{sample_name}_6mA_100kb.bed", 
        mod5hmC_100kb = "methylation/{sample_name}/{sample_name}_5hmC_100kb.bed", 
        mod5mC_100kb = "methylation/{sample_name}/{sample_name}_5mC_100kb.bed",
    conda: 
        "envs/methylation_change.yaml"
    shell:
        """
        bedtools map -a {input.chr_windows} -b  {input.mod6mA} -c 11 -o mean > {output.mod6mA_100kb}
        bedtools map -a {input.chr_windows} -b  {input.mod5hmC} -c 11 -o mean > {output.mod5hmC_100kb}
        bedtools map -a {input.chr_windows} -b  {input.mod5mC} -c 11 -o mean > {output.mod5mC_100kb}
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