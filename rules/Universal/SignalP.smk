#SignalP

import glob
import os

# Split sequences for signalp v5.0 (max. 5000 seq)
checkpoint splittingsignalP:
    input:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/renamed/{sample}_ID.faa")
    output:
        splits=temp(directory(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/SignalP/splitted/{sample}_dir/")))
    log:
        os.path.join(DATA_DIR,"{project}/logs/{sample}/SignalP_split.log")
    params:
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"],
        split=config["pathofact"]["size_fasta"]
    conda:                                          
        "../../envs/Biopython.yaml"                                                
    shell:                                                      
        """                                                                 
        python {config[pathofact][scripts]}/split.py {input} 2000 {output} &> {log}
        """


#Run SignalP on split sequence files
rule signalp_gramp:
    input:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/SignalP/splitted/{sample}_dir/{file_i}.fasta")        
    output:
        SignalP_gramP=temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/SignalP/Gram+/{sample}/{file_i}_summary.signalp5"))
    log:
        os.path.join(DATA_DIR,"{project}/logs/{sample}/SignalP_p_{file_i}.log")
    message:
        "Execute signalP on the following sample(s): {wildcards.project} - {wildcards.sample}"
    params:
        runtime=config["pathofact"]["runtime"]["long"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    shell:
        """
	export PATH={config[pathofact][signalp]}:$PATH
        signalp -fasta {input} -org gram+ -prefix $(realpath {output} | sed 's/_summary.signalp5//g') &> {log}
        """

rule signalp_gramn:
    input:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/SignalP/splitted/{sample}_dir/{file_i}.fasta")
    output:
        SignalP_gramN=temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/SignalP/Gram-/{sample}/{file_i}_summary.signalp5"))
    log:
        os.path.join(DATA_DIR,"{project}/logs/{sample}/SignalP_n_{file_i}.log")
    message:
        "Execute signalP on the following sample(s): {wildcards.project} - {wildcards.sample}"
    params:
        runtime=config["pathofact"]["runtime"]["long"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    shell:
        """
        export PATH={config[pathofact][signalp]}:$PATH
        signalp -fasta {input} -org gram- -prefix $(realpath {output} | sed 's/_summary.signalp5//g') &> {log}
        """

def aggregate_signalpP_input(wildcards):
    checkpoint_output= checkpoints.splittingsignalP.get(**wildcards).output.splits
    return expand(
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/SignalP/Gram+/{sample}/{file_i}_summary.signalp5"),
        project=wildcards.project,
        sample=wildcards.sample,
        file_i=glob_wildcards(os.path.join(checkpoint_output, "{i}.fasta")).i
    )

def aggregate_signalpN_input(wildcards):
    checkpoint_output= checkpoints.splittingsignalP.get(**wildcards).output.splits
    return expand(
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/SignalP/Gram-/{sample}/{file_i}_summary.signalp5"),
        project=wildcards.project,
        sample=wildcards.sample,
        file_i=glob_wildcards(os.path.join(checkpoint_output, "{i}.fasta")).i
    )

rule SignalPP_aggregate:
    input:
        aggregate_signalpP_input
    output:
        temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/SignalP/{sample}/gramp_summary.signalp5"))
    params:
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    shell:
        "cat {input} > {output}"

rule SignalPN_aggregate:
    input:
        aggregate_signalpN_input
    output:
        temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/SignalP/{sample}/gramn_summary.signalp5"))
    params:
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    shell:
        "cat {input} > {output}"

rule aggregate_signalP:
    input:
        SignalP_gramP=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/SignalP/{sample}/gramp_summary.signalp5"),
        SignalP_gramN=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/SignalP/{sample}/gramn_summary.signalp5")
    output:
        SignalP_report=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/SignalP/aggregated/{sample}_SignalP_results.tsv")
    message:
        "concatenate multiple split signalP files in a single joined file: {wildcards.project} - {wildcards.sample}"
    log:
        os.path.join(DATA_DIR,"{project}/logs/{sample}/SignalP_temp.log")
    params:
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    conda: "../../envs/R.yaml"    
    script: "../../scripts/SignalP.R"
