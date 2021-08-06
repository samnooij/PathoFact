#Toxin

import glob
import os

# HMM scan
rule run_HMM_tox:
    input:
        hmm=config["pathofact"]["tox_hmm"],
        renamed=os.path.join(DATA_DIR,"{project}/splitted/{sample}/{file_i}.fasta")
    output:
        temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/TOXIN/HMM_toxin/{sample}/{file_i}.hmmscan"))
    log:
        os.path.join(DATA_DIR,"{project}/logs/{sample}/{file_i}.log")
    params:
        runtime=config["pathofact"]["runtime"]["long"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    conda:
        "../../envs/HMMER.yaml"
    threads:
        1
    message: "Executing toxin prediction with {threads} threads on the following sample(s): {wildcards.project} - {wildcards.sample}"
    shell:
        """
        hmmsearch --cpu {threads} --noali --notextw -T {config[pathofact][tox_threshold]} --tblout {output} {input.hmm} {input.renamed} &> {log}
        """

# Adjust HMM results to correct format
rule HMM_correct_format:
    input:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/TOXIN/HMM_toxin/{sample}/{file_i}.hmmscan")
    output:
        temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/TOXIN/HMM_toxin/{sample}/{file_i}.hmm.csv"))
    params:
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    message: "Adjust output format of toxin predictions: {wildcards.project} - {wildcards.sample}"
    shell:
        """
        sed '/^#/ d' {input} | sed 's/ \+/\\t/g' > {output}
        """

def aggregate_hmm(wildcards):
    checkpoint_output = checkpoints.splitting.get(**wildcards).output.splits
    return expand(
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/TOXIN/HMM_toxin/{sample}/{file_i}.hmm.csv"),
        project=wildcards.project,
        sample=wildcards.sample,
        file_i=glob_wildcards(os.path.join(checkpoint_output, "{i}.fasta")).i
    )

rule HMM_correct_format_2:
    input:
        aggregate_hmm
    output:
        temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/TOXIN/HMM_toxin/{sample}.Input_HMM_R_temp.csv"))
    params:
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    message: "Aggregate toxin prediction of the following sample(s): {wildcards.project} - {wildcards.sample}"
    shell:
        """
        cut -f 1,3,5,6 {input} | uniq > {output}
        """

rule HMM_correct_format_3:
    input:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/TOXIN/HMM_toxin/{sample}.Input_HMM_R_temp.csv")
    output:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/TOXIN/HMM_toxin/{sample}.Input_HMM_R.csv")
    params:
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    shell:
        """
        echo "#Toxin" > {config[pathofact][datadir]}/{wildcards.project}/PathoFact_intermediate/TOXIN/HMM_toxin/{wildcards.sample}_header
        cat {config[pathofact][datadir]}/{wildcards.project}/PathoFact_intermediate/TOXIN/HMM_toxin/{wildcards.sample}_header {input} > {output}
        rm -rf {config[pathofact][datadir]}/{wildcards.project}/PathoFact_intermediate/TOXIN/HMM_toxin/{wildcards.sample}_header
        sed -i $'1 i\\\ Query_sequence\\tHMM_Name\\tSignificance_Evalue\\tScore' {output}    
        """
