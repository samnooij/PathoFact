#Prepare fasta

import glob
import os

##############################
#     Modify fasta input     #
##############################

# Generate unique ID number for each sequence
rule generate_contigID:
    input:
        os.path.join(DATA_DIR,"{sample}.fna")
    output:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/renamed/{sample}_Contig_ID.fna")
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
rule generate_ContigTranslation:
    input:
        renamed=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/renamed/{sample}_Contig_ID.fna"),
        original=os.path.join(DATA_DIR,"{sample}.fna")
    output:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/renamed/{sample}_Contig_translation.tsv")
    message:
        "Generate {output} containing original fasta header with corresponding ID number: {wildcards.project} - {wildcards.sample}"
    params:
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    shell:
        """
        paste {input.renamed} {input.original} | awk 'sub(/^>/,"")' OFS='\\t' > {output}
        """

checkpoint splitcontig:
    input:
        os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/renamed/{sample}_Contig_ID.fna")
    output:
        split=temp(directory(os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/contig_splitted/{sample}/")))
    params:
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"],
        split=config["pathofact"]["size_fasta"]
    conda:
        "../../envs/Biopython.yaml"
    log:
        os.path.join(DATA_DIR,"{project}/logs/{sample}/split_contig.log")
    shell:
        """
         python {config[pathofact][scripts]}/split.py {input} {params.split} {output.split} &> {log}
        """


