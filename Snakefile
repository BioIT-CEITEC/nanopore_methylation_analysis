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

##### Target rules #####
rule all:
    input:
        expand("SV_calling/{sample_name}/variants.vcf", sample_name = sample_tab.sample_name)

### separate bams for all chromosomes in reference
rule create_bam_for_chromosome:
    input: 
        bam = 'aligned/{sample_name}/{sample_name}_sorted.bam'
    output:
        bam = 'aligned/{sample_name}/{sample_name}_sorted.REF_chr1'
    conda: 
        "envs/methylation_change.yaml"
    shell:
        """
        bamtools split -in {input.bam} -reference
        """

rule find_highly_modified_motifs:
    input: 
        bed = "methylation/{sample_name}/{sample_name}_modkit.bed"
    output: 
        log = "methylation/{sample_name}/{sample_name}_modkit_find_motifs_log.txt",
        tsv = "methylation/{sample_name}/{sample_name}_motifs.tsv"
    params: 
        genome = config["organism_fasta"],
    conda: 
        "envs/methylation.yaml"
    shell: 
        """
        modkit motif search -i {input.bed} -r {params.genome} -o {output.tsv} --threads 32 --log {output.log}
        """


### separate bams for modifications
rule extract_separate_modifications:
    input: 
        bam = 'aligned/{sample_name}/{sample_name}_sorted.REF_{chromosome}.bam'
    output: 
        filtered_bam = 'aligned/{sample_name}/{modification}/{sample_name}_sorted.REF_{chromosome}_{modification}.bam'
    params: 
        mod5mC_5hmC= lambda wildcards, input: int(config["5mC_methylation"]), # convert to 0/1 for bash
        mod6mA= lambda wildcards, input: int(config["6mA_methylation"])
    conda: 
        "envs/methylation_change.yaml"
    shell:
        """
        for file in *_chr*.bam; do
            output="5mC/${file%.bam}_5mC.bam"
            samtools view -h "$file" | grep -E "MM:Z:[^;]*5mC" | samtools view -bS > "$output"
            echo "Zpracován: $file -> $output"
        done
        if [ {params.mod5mC_5hmC} -eq 1 ]; then
            samtools view -h {input.bam} | grep -E "MM:Z:[^;]*5mC" | samtools view -bS > filtered_6mA.bam
        """

#rule QC_after_SV_calling:
