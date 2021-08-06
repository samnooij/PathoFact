#AMR

import glob
import os

##########################
#     AMR Prediction     #
##########################

## deepARG:
rule run_deepARG:
    input:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/splitted/{sample}/{file_i}.fasta")
    output:
        temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/AMR/deepARG_results/{sample}/{file_i}.out.mapping.ARG"))
    log:
        os.path.join(DATA_DIR,"{project}/logs/{sample}/{file_i}.out.mapping.ARG.log")
    threads: 24
    params:
        runtime=config["pathofact"]["runtime"]["medium"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    conda:
        "../../envs/DeepARG.yaml"
    message: "executing deep-arg on the following sample(s): {wildcards.project} - {wildcards.sample}"
    shell:
        """
        deeparg predict --model LS --model-version v2 --type prot -d {config[pathofact][scripts]}/deeparg_data/deepARG --input {input} --out {config[pathofact][datadir]}/{wildcards.project}/PathoFact_intermediate/AMR/deepARG_results/{wildcards.sample}/{wildcards.file_i}.out &> {log}
        """

def aggregate_AMR(wildcards):
    checkpoint_output = checkpoints.splitting.get(**wildcards).output.splits
    return expand(
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/AMR/deepARG_results/{sample}/{file_i}.out.mapping.ARG"),
        project=wildcards.project,
        sample=wildcards.sample,
        file_i=glob_wildcards(os.path.join(checkpoint_output, "{i}.fasta")).i
    )

rule aggregate_deepARG:
    input:
        aggregate_AMR
    output:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/AMR/deepARG_results/{sample}.out.mapping.ARG")
    params:
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    message: "aggregating deep-arg results on the following sample(s): {wildcards.project} - {wildcards.sample}"
    shell:
        """
        cat {input} > {output}
        rm -rf {config[pathofact][datadir]}/{wildcards.project}/PathoFact_intermediate/AMR/deepARG_results/{wildcards.sample}
        """

# RGI
rule run_RGI:
    input:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/splitted/{sample}/{file_i}.fasta")
    output:
        temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/AMR/RGI_results/{sample}/{file_i}.RGI.txt"))
    log:
        os.path.join(DATA_DIR,"{project}/logs/{sample}/{file_i}.RGI.log")
    params:
        runtime=config["pathofact"]["runtime"]["medium"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    conda:
        "../../envs/rgi.yaml"
    message: "executing RGI on the following sample(s): {wildcards.project} - {wildcards.sample}"
    shell:
        "rgi main --input_sequence {input} --output_file {config[pathofact][datadir]}/{wildcards.project}/PathoFact_intermediate/AMR/RGI_results/{wildcards.sample}/{wildcards.file_i}.RGI --input_type protein --local --clean  &> {log}"

def aggregate_RGI(wildcards):
    checkpoint_output = checkpoints.splitting.get(**wildcards).output.splits
    return expand(
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/AMR/RGI_results/{sample}/{file_i}.RGI.txt"),
        project=wildcards.project,
        sample=wildcards.sample,
        file_i=glob_wildcards(os.path.join(checkpoint_output, "{i}.fasta")).i
    )

rule aggregate_RGI:
    input:
        aggregate_RGI
    output:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/AMR/RGI_results/{sample}.RGI.txt")
    params:
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    message: "Aggregate RGI results on the following sample(s): {wildcards.project} - {wildcards.sample}"
    shell:
        """
        cat {input} > {output}
        rm -rf {config[pathofact][datadir]}/{wildcards.project}/PathoFact_intermediate/AMR/RGI_results/{wildcards.sample}
        """

# Combine DeepARG and RGI results

rule combine_AMR:
    input:
        RGI=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/AMR/RGI_results/{sample}.RGI.txt"),
        DeepARG=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/AMR/deepARG_results/{sample}.out.mapping.ARG"),
        AMR_index= "scripts/PathoFact_AMR_index.tsv"
    output:
        AMR_combined=os.path.join(DATA_DIR,"{project}/AMR/{sample}_AMR_prediction.tsv")
    log:
        os.path.join(DATA_DIR,"{project}/logs/{sample}/combine_AMR_temp.log")
    params:
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    conda: "../../envs/R.yaml"
    message: "Combine AMR prediction of RGI and DeepArg for the following sample: {wildcards.project} - {wildcards.sample}"
    script:
        "../../scripts/AMR.R"
