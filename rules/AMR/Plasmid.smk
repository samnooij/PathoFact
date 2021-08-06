#Plasmid

import glob
import os

##########################
#   Plasmid Prediction   #
##########################

# PlasFlow Preprocessing
rule filter_seq:
    input:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/renamed/{sample}_Contig_ID.fna")
    output:
        temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/plasmid/{sample}_filtered.fna"))
    log:
        os.path.join(DATA_DIR,"{project}/logs/{sample}/plasmid_filtered.log")
    conda:
        "../../envs/Biopython.yaml"
    params:
        runtime=config["pathofact"]["runtime"]["medium"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"],
        minlen=config["pathofact"]["plasflow_minlen"]
    message: "Filter samples on length for PlasFlow predictions: {wildcards.project} - {wildcards.sample}"
    shell:
        "{config[pathofact][scripts]}/filter.pl {params.minlen} {input} > {output} 2> {log}"

checkpoint splitplasmid:
    input:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/plasmid/{sample}_filtered.fna")
    output:
        split=temp(directory(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/MGE/plasmid_splitted/{sample}/")))
    params:
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"],
        split=config["pathofact"]["size_fasta"]
    conda:
        "../../envs/Biopython.yaml"
    shell:
        """
         python {config[pathofact][scripts]}/split.py {input} {params.split} {output.split}
        """

# PlasFlow Plasmid prediction
rule run_PLASMID:
    input:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/MGE/plasmid_splitted/{sample}/{file_i}.fasta")
    output:
        temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/MGE/plasmid/PlasFlow/{sample}/{file_i}_plasflow_prediction.tsv"))
    log:
        os.path.join(DATA_DIR,"{project}/logs/{sample}/{file_i}_plasflow_prediction.log")
    conda:
        "../../envs/PlasFlow.yaml"
    params:
        runtime=config["pathofact"]["runtime"]["long"],
        mem=config["pathofact"]["mem"]["big_mem_per_core_gb"],
        threshold=config["pathofact"]["plasflow_threshold"]
    message: "Executing PlasFLow on the following sample(s): {wildcards.project} - {wildcards.sample}"
    shell:
        """
        PlasFlow.py --input {input} --output {output} --threshold {params.threshold} &> {log}
        """

def aggregate_plasmid_input(wildcards):
    checkpoint_output= checkpoints.splitplasmid.get(**wildcards).output.split
    return expand(
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/MGE/plasmid/PlasFlow/{sample}/{file_i}_plasflow_prediction.tsv"),
        project=wildcards.project,
        sample=wildcards.sample,
        file_i=glob_wildcards(os.path.join(checkpoint_output, "{i}.fasta")).i
    )

rule Plasmid_aggregate:
    input:
        aggregate_plasmid_input
    output:
        temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/MGE/plasmid/PlasFlow/{sample}_plasflow_aggregated.tsv"))
    params:
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    message: "Aggregate PlasFlow results on the following sample(s): {wildcards.project} - {wildcards.sample}"
    shell:
        """
        cat {input} > {output}
        rm -rf {config[pathofact][datadir]}/{wildcards.project}/PathoFact_intermediate/MGE/plasmid/PlasFlow/{wildcards.sample}
        """

rule select:
    input:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/MGE/plasmid/PlasFlow/{sample}_plasflow_aggregated.tsv")
    output:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/MGE/plasmid/PlasFlow/{sample}_plasflow_prediction_final.tsv")
    params:
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    shell:
        """
        cut -f 3,6 {input} > {output}
        """

#rule run_MOBsuite:
#    input: os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/contig_splitted/{sample}/{file_i}.fasta")
#    output:         
#        temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/MGE/plasmid/MOB_suite/{sample}/{file_i}_MOB_suite_prediction.txt"))
#    log:
#        os.path.join(DATA_DIR,"{project}/logs/{sample}/{file_i}_MOB_suite_prediction.log")
#    threads:
#        config["pathofact"]["mem"]["big_mem_cores"]
#    conda:
#       "../../envs/MOB_suite.yaml"
#    message: "Executing MOB_suite with {threads} threads on the following sample(s): {wildcards.project} - {wildcards.sample}"
#    shell: "mob_typer --multi --infile {input} --out_file {output} -n {threads} &> {log}"

#def aggregate_MOBsuite(wildcards):
#    checkpoint_output= checkpoints.splitcontig.get(**wildcards).output.split
#    return expand(
#        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/MGE/plasmid/MOB_suite/{sample}/{file_i}_MOB_suite_prediction.txt"),
#        project=wildcards.project,
#        sample=wildcards.sample,
#        file_i=glob_wildcards(os.path.join(checkpoint_output, "{i}.fasta")).i
#    )

#rule aggregate_MOBsuite:
#    input: aggregate_MOBsuite
#    output:
#        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/MGE/plasmid/MOB_suite/{sample}_MOB_suite_aggregated.tsv")
#    params:
#        runtime=config["pathofact"]["runtime"]["short"],
#        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
#    message: "Aggregate MOB_suite results on the following sample(s): {wildcards.project} - {wildcards.sample}"
#    shell:
#        "cat {input} > {output}"

