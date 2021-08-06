#Prepare fasta

import glob
import os

##########################################
#     Generate ORFs and mapping file     #
##########################################

rule Prodigal:
    input:
        os.path.join(DATA_DIR,"{sample}.fna")
    output:
        ORF=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/Prodigal/{sample}.faa"),
        GFF=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/Prodigal/{sample}.gff")
    message:
        "Generates ORFs and gff: {wildcards.project} - {wildcards.sample}"
    conda:
        "../../envs/Prodigal.yaml"
    log: 
        os.path.join(DATA_DIR,"{project}/logs/{sample}/Prodigal.log")
    params:
        runtime=config["pathofact"]["runtime"]["medium"],
        mem=config["pathofact"]["mem"]["big_mem_per_core_gb"]
    shell:
        """
        prodigal -i {input} -o {output.GFF} -a {output.ORF} -f gff -p meta &> {log}
        sed -i 's/[^>]*ID=//;s/;.*//;s/*//' {output.ORF}
        """

rule mapping_file:
    input:
        ORF=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/Prodigal/{sample}.faa"),
        GFF=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/Prodigal/{sample}.gff")
    output:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/Prodigal/{sample}.contig")
    message:
        "Generate mapping file: {wildcards.project} - {wildcards.sample}"
    params:
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    shell:
        """
        sed -i '/^#/d' {input.GFF}
        cut -f 1,9 {input.GFF} |cut -d';' -f1| sed 's/ID=//' > {output}
        """

##############################
#     Modify fasta input     #
##############################

# Generate unique ID number for each sequence
rule generate_ID:
    input:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/Prodigal/{sample}.faa")
    output:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/renamed/{sample}_ID.faa")
    message:
        "Replace fasta headers with unique ID number: {wildcards.project} - {wildcards.sample}"
    params:
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    shell:
        """
        awk 'BEGIN{{zeros="0000000000"}}{{if(substr($1,1,1)==">"){{i+=1;print">"substr(zeros,1,10-length(i))""i}}else{{print$0}}}}' {input}  > {output}
        """

# Generate translation file combining original header with unique ID
rule generate_translation:
    input:
        renamed=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/renamed/{sample}_ID.faa"),
        original=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/Prodigal/{sample}.faa")
    output:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/renamed/{sample}_translation.tsv")
    message:
        "Generate {output} containing original fasta header with corresponding ID number: {wildcards.project} - {wildcards.sample}"
    params:
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    shell:
        """
        paste {input.renamed} {input.original} | awk 'sub(/^>/,"")' OFS='\\t' > {output}
        """

###############################
#   Checkpoint split fasta    #
###############################

# Split fasta file
#checkpoint splitting:
checkpoint splitting:
    input:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/renamed/{sample}_ID.faa")
    output:
        splits=temp(directory(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/splitted/{sample}/")))
    log:
        os.path.join(DATA_DIR,"{project}/logs/{sample}/split_ORF.log")
    params:
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"],
        split=config["pathofact"]["size_fasta"]
    conda:
        "../../envs/Biopython.yaml"
    shell:
        """
        python {config[pathofact][scripts]}/split.py {input} {params.split} {output.splits} &> {log}
        """
