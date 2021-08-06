#Plasmid

import glob
import os

##########################
#     Phage Prediction   #
##########################

# VIRSORTER 

rule run_VirSorter:
    input:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/renamed/{sample}_Contig_ID.fna")
    output:
        temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/MGE/phage/{sample}/virsorter/VIRSorter_global-phage-signal.csv"))
    log:
        os.path.join(DATA_DIR,"{project}/logs/{sample}/VIRSorter_global-phage-signal.log")
    params:
        runtime=config["pathofact"]["runtime"]["long"],
        mem=config["pathofact"]["mem"]["big_mem_per_core_gb"]
    conda:
        "../../envs/VirSorter.yaml"
    threads:
        config["pathofact"]["mem"]["big_mem_cores"]
    message: "Executing VirSorter with {threads} threads on the following sample(s): {wildcards.project} - {wildcards.sample}"
    shell:
        """
        wrapper_phage_contigs_sorter_iPlant.pl -f {input} --ncpu {threads} --wdir {config[pathofact][datadir]}/{wildcards.project}/PathoFact_intermediate/MGE/phage/{wildcards.sample}/virsorter --data-dir {config[pathofact][scripts]}/virsorter-data &> {log}
        """

localrules: aggregate_VirSorter
rule aggregate_VirSorter:
    input:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/MGE/phage/{sample}/virsorter/VIRSorter_global-phage-signal.csv")
    output:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/MGE/phage/{sample}_VIRSorter_aggregated.csv")
    message: "VirSorter failsave for empty files: {wildcards.project} - {wildcards.sample}"
    shell:
        """
        if [ -s {output} ]
        then
           cp scripts/virsorter_headers.csv {output}
        else
            mv {input} {output}
        fi
        rm -rf {config[pathofact][datadir]}/{wildcards.project}/PathoFact_intermediate/MGE/phage/{wildcards.sample}/virsorter
        """

# VIRFINDER Prediction
rule run_VirFinder:
    input:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/contig_splitted/{sample}/{file_i}.fasta")
    output:
        temp(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/MGE/phage/{sample}/virfinder/{file_i}.fasta_gt1bp_dvfpred.txt"))
    log:
        os.path.join(DATA_DIR,"{project}/logs/{sample}/{file_i}.fasta_gt1bp_dvfpred.log")
    params:
        runtime=config["pathofact"]["runtime"]["long"],
        mem=config["pathofact"]["mem"]["big_mem_per_core_gb"]
    conda:
        "../../envs/DeepVirFinder.yaml"
    threads: config["pathofact"]["mem"]["big_mem_cores"]
    message: "Executing Deep-VirFinder with {threads} threads on the following sample(s): {wildcards.project} - {wildcards.sample}"
    shell:
        "python {config[pathofact][deepvirfinder]} -i {input} -o {config[pathofact][datadir]}/{wildcards.project}/PathoFact_intermediate/MGE/phage/{wildcards.sample}/virfinder -c {threads} &> {log}"

def aggregate_VirFinder(wildcards):
    checkpoint_output= checkpoints.splitcontig.get(**wildcards).output.split
    return expand(
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/MGE/phage/{sample}/virfinder/{file_i}.fasta_gt1bp_dvfpred.txt"),
        project=wildcards.project,
        sample=wildcards.sample,
        file_i=glob_wildcards(os.path.join(checkpoint_output, "{i}.fasta")).i
    )

rule aggregate_VirFinder:
    input:
        aggregate_VirFinder
    output:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/MGE/phage/{sample}_VirFinder_aggregated.csv")
    params:
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    message: "Aggregate VirFinder predictions on the following sample(s): {wildcards.project} - {wildcards.sample}"
    shell:
        "cat {input} >{output}"
