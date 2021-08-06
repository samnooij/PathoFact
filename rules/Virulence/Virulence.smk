#Virulence

import glob
import os

# HMM scan
rule run_HMM_vir:
    input:
        hmm=config["pathofact"]["vir_hmm"],
        renamed=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/splitted/{sample}/{file_i}.fasta")
    output:
        temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/HMM_virulence/{sample}/{file_i}.hmmscan"))
    log:
        os.path.join(DATA_DIR,"{project}/logs/{sample}/{file_i}.log")
    message:
        "Run HMM scan on {input.renamed} to generate {output}"
    params:
        runtime=config["pathofact"]["runtime"]["long"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    conda:
        "../../envs/HMMER.yaml"
    threads:
        12
    message: "Executing prediction of virulence factors (HMM) with {threads} threads on the following sample(s): {wildcards.project} - {wildcards.sample}"
    shell:
        """
        hmmsearch --cpu {threads} --noali --notextw --tblout {output} {input.hmm} {input.renamed} &> {log}
        """

# Adjust HMM results to correct format
rule HMM_correct_format_vir:
    input:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/HMM_virulence/{sample}/{file_i}.hmmscan")
    output:
        temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/HMM_virulence/{sample}/{file_i}.hmm.csv"))
    message:
        "Adjust {input} to correct format: {output}"
    params:
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    message: "Format virulence factor predictions (HMM): {wildcards.project} - {wildcards.sample}"
    shell:
        """
        sed '/^#/ d' {input} | sed 's/ \+/\\t/g' > {output}
        """

def aggregate_hmm(wildcards):
    checkpoint_output = checkpoints.splitting.get(**wildcards).output.splits
    return expand(
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/HMM_virulence/{sample}/{file_i}.hmm.csv"),
        project=wildcards.project,
        sample=wildcards.sample,
        file_i=glob_wildcards(os.path.join(checkpoint_output, "{i}.fasta")).i
    )

rule HMM_correct_format_2_vir:
    input:
        aggregate_hmm
    output:
        temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/HMM_virulence/{sample}.Input_HMM_R.csv"))
    message: "Aggregate virulence factor prediction (HMM) on the following sample(s): {wildcards.project} - {wildcards.sample}"
    params:
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    shell:
        """
        cut -f 1,3,5,6 {input} |uniq > {output}
        """

rule HMM_R_VIR:
    input:
        hmmresults=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/HMM_virulence/{sample}.Input_HMM_R.csv"),
        positive=config["pathofact"]["vir_domains"] + "/positive_domains.tsv",
        negative=config["pathofact"]["vir_domains"] + "/negative_domains.tsv",
        shared=config["pathofact"]["vir_domains"] + "/shared_domains.tsv",
        ID=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/renamed/{sample}_translation.tsv")
    output:
        temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/HMM_virulence/{sample}.hmm_results.csv"))
    log:
        os.path.join(DATA_DIR,"{project}/logs/{sample}/hmm_results.log")
    params:
        runtime=config["pathofact"]["runtime"]["medium"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    conda:
        "../../envs/R.yaml"
    script:
        "../../scripts/hmm.R"

# Give pathogenicity prediction
rule HMM_VIR_classification:
    input:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/HMM_virulence/{sample}.hmm_results.csv")
    output:
        non_path=temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/HMM_virulence/{sample}_hmm_non_path.txt")),
        pathogenic=temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/HMM_virulence/{sample}_hmm_pathogenic.txt")),
        unclassified=temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/HMM_virulence/{sample}_hmm_unclassified.txt"))
    params:
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    message: "Classify HMM virulence predictions on the following sample(s): {wildcards.project} - {wildcards.sample}"
    shell:
        """
        awk '$2 == "True" && $3 == "False"  && $4 == "False"' {input} | awk '$5 = "negative"' | sed 's/ /\\t/g' > {output.non_path}
        awk '$3 == "True"' {input} | awk '$5 = "pathogenic"' | sed 's/ /\\t/g' > {output.pathogenic}
        awk '$3 == "False"  && $4 == "True"' {input} | awk '$5 = "unclassified"' | sed 's/ /\\t/g' > {output.unclassified}
        """


rule HMM_VIR_report:
    input:
        non_pathogenic= os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/HMM_virulence/{sample}_hmm_non_path.txt"),
        pathogenic= os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/HMM_virulence/{sample}_hmm_pathogenic.txt"),
        unclassified= os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/HMM_virulence/{sample}_hmm_unclassified.txt")
    output:
        temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/HMM_virulence/{sample}_hmm_prediction_intermediate.tsv"))
    params:
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    shell:
        """
        cat {input} | sort > {output}
        """

rule HMM_VIR_finalformat:
    input:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/HMM_virulence/{sample}_hmm_prediction_intermediate.tsv")
    output:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/HMM_virulence/{sample}_hmm_prediction.tsv")
    params:
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    shell:
        "cut -f 1,5 {input}  > {output}"

##########################
#  Virulence classifier  #
##########################
rule AAC:
    input:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/splitted/{sample}/{file_i}.fasta")
    output:
        temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/classifier_virulence/{sample}/{file_i}_AAC.txt"))
    log:
        os.path.join(DATA_DIR,"{project}/logs/{sample}/{file_i}_AAC.log")
    params:
        runtime=config["pathofact"]["runtime"]["medium"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    conda:
        "../../envs/Biopython.yaml"
    message: "Identify features (AAC) on the following sample(s): {wildcards.project} - {wildcards.sample}"
    shell:
        "python {config[pathofact][scripts]}/AAC.py --file {input} --out {output} &> {log}"

rule DPC:
    input:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/splitted/{sample}/{file_i}.fasta")
    output:
        temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/classifier_virulence/{sample}/{file_i}_DPC.txt"))
    log:
        os.path.join(DATA_DIR,"{project}/logs/{sample}/{file_i}_DPC.log")
    params:
        runtime=config["pathofact"]["runtime"]["medium"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    conda:
        "../../envs/Biopython.yaml"
    message: "Identify features (DPC) on the following sample(s): {wildcards.project} - {wildcards.sample}"
    shell:
        "python {config[pathofact][scripts]}/DPC.py --file {input} --out {output} &> {log}"

rule CTDC:
    input:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/splitted/{sample}/{file_i}.fasta")
    output:
        temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/classifier_virulence/{sample}/{file_i}_CTDC.txt"))
    log:
        os.path.join(DATA_DIR,"{project}/logs/{sample}/{file_i}_CTDC.log")
    params:
        runtime=config["pathofact"]["runtime"]["medium"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    conda:
        "../../envs/Biopython.yaml"
    message: "Identify features (CTDC) on the following sample(s): {wildcards.project} - {wildcards.sample}"
    shell:
        "python {config[pathofact][scripts]}/CTDC.py --file {input} --out {output} &> {log}"

rule CTDT:
    input:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/splitted/{sample}/{file_i}.fasta")
    output:
        temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/classifier_virulence/{sample}/{file_i}_CTDT.txt"))
    log:
        os.path.join(DATA_DIR,"{project}/logs/{sample}/{file_i}_CTDT.log")
    params:
        runtime=config["pathofact"]["runtime"]["medium"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    conda:
        "../../envs/Biopython.yaml"
    message: "Identify features (CTDT) on the following sample(s): {wildcards.project} - {wildcards.sample}"
    shell:
        "python {config[pathofact][scripts]}/CTDT.py --file {input} --out {output} &> {log}"

rule CTDD:
    input:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/splitted/{sample}/{file_i}.fasta")
    output:
        temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/classifier_virulence/{sample}/{file_i}_CTDD.txt"))
    log:
        os.path.join(DATA_DIR,"{project}/logs/{sample}/{file_i}_CTDD.log")
    params:
        runtime=config["pathofact"]["runtime"]["medium"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    conda:
        "../../envs/Biopython.yaml"
    message: "Identify features (CTDD) on the following sample(s): {wildcards.project} - {wildcards.sample}"
    shell:
        "python {config[pathofact][scripts]}/CTDD.py --file {input} --out {output} &> {log}"

rule join_matrix:
    input:
        AAC=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/classifier_virulence/{sample}/{file_i}_AAC.txt"),
        DPC=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/classifier_virulence/{sample}/{file_i}_DPC.txt"),
        CTDC=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/classifier_virulence/{sample}/{file_i}_CTDC.txt"),
        CTDT=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/classifier_virulence/{sample}/{file_i}_CTDT.txt"),
        CTDD=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/classifier_virulence/{sample}/{file_i}_CTDD.txt")
    output:
        temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/classifier_virulence/{sample}/{file_i}_matrix.tsv"))
    params:
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    message: "Join feature matrix on the following sample(s): {wildcards.project} - {wildcards.sample}"
    shell:
        """
        xjoin() {{
            local f
                local srt="sort -k 1b,1"

                if [ "$#" -lt 2 ]; then
                        echo "xjoin: need at least 2 files" >&2
                        return 1
                elif [ "$#" -lt 3 ]; then
                        join -t $'\\t' <($srt "$1") <($srt "$2")
                else
                        f=$1
                        shift
                        join -t $'\\t' <($srt "$f") <(xjoin "$@")
                fi
        }}

        xjoin {input} | sort -k 1n,1 > {output}
        """

rule classifier:
    input:
        input=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/classifier_virulence/{sample}/{file_i}_matrix.tsv"),
        model=config["pathofact"]["scripts"] + "/Virulence_factor_model.sav"
    output:
        temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/classifier_virulence/{sample}/{file_i}_classifier_prediction.tsv"))
    log:
        os.path.join(DATA_DIR,"{project}/logs/{sample}/{file_i}_classifier_prediction.log")
    params:
        runtime=config["pathofact"]["runtime"]["medium"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    conda:
        "../../envs/Biopython.yaml"
    message: "Execute classifier for virulence factor prediction on the following sample(s): {wildcards.project} - {wildcards.sample}"
    shell:
        "python {config[pathofact][scripts]}/virulence_prediction.py {input.input} {output} {input.model} &> {log}"

def aggregate_classifier(wildcards):
    checkpoint_output = checkpoints.splitting.get(**wildcards).output.splits
    return expand(
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/classifier_virulence/{sample}/{file_i}_classifier_prediction.tsv"),
        project=wildcards.project,
        sample=wildcards.sample,
        file_i=glob_wildcards(os.path.join(checkpoint_output, "{i}.fasta")).i
    )

rule format_classifier:
    input:
        aggregate_classifier
    output:
        temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/classifier_virulence/{sample}_classifier_results_format.tsv"))
    params:
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    message: "Aggregate classifier predictions for virulence factors on the following sample(s): {wildcards.project} - {wildcards.sample}"
    shell:
        """
        sed 's/"//g' {input} | cut -f2,3 >{output}
        """

rule format_classifier_2:
    input:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/classifier_virulence/{sample}_classifier_results_format.tsv")
    output:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/classifier_virulence/{sample}_classifier_results_formatted.tsv")
    params:
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    shell:
        """
        awk '{{$1=sprintf("%010d", $1)}}1' {input} | sed 's/ /\\t/g' > {output}
        """
