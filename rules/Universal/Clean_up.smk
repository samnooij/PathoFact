#Clean-up
import glob
import os

rule clean_all:
    input: os.path.join(DATA_DIR,"{project}/PathoFact_report/PathoFact_{sample}_predictions.tsv")
    output: os.path.join(DATA_DIR,"{project}/logs/{sample}_compressed.zip")
    shell: """
        zip -rm {output} {config[pathofact][datadir]}/{wildcards.project}/logs/{wildcards.sample}
        rm -rf {config[pathofact][datadir]}/{wildcards.project}/splitted/{wildcards.sample}
        rm -rf {config[pathofact][datadir]}/{wildcards.project}/contig_splitted/{wildcards.sample}_dir
        rm -rf {config[pathofact][datadir]}/{wildcards.project}/PathoFact_intermediate/SignalP/splitted/{wildcards.sample}
        rm -rf {config[pathofact][datadir]}/{wildcards.project}/PathoFact_intermediate/MGE/plasmid_splitted/{wildcards.sample}
        find {config[pathofact][datadir]}/{wildcards.project} -type d -empty -delete
        """

