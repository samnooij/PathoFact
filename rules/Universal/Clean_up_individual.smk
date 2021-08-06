#Clean-up
import glob
import os

rule clean_Toxin_workflow:
    input: os.path.join(DATA_DIR,"{project}/PathoFact_report/Toxin_prediction_{sample}_report.tsv")
    output: os.path.join(DATA_DIR,"{project}/logs/Tox_{sample}_compressed.zip")
    shell: """
        zip -rm {output} {config[pathofact][datadir]}/{wildcards.project}/logs/{wildcards.sample}
        rm -rf {config[pathofact][datadir]}/{wildcards.project}/splitted/{wildcards.sample}
        rm -rf {config[pathofact][datadir]}/{wildcards.project}/PathoFact_intermediate/SignalP/splitted/{wildcards.sample}
        find {config[pathofact][datadir]}/{wildcards.project} -type d -empty -delete
        """

rule clean_VF_workflow:
    input: os.path.join(DATA_DIR,"{project}/PathoFact_report/Virulence_prediction_{sample}_report.tsv")
    output: os.path.join(DATA_DIR,"{project}/logs/VF_{sample}_compressed.zip")
    shell: """
        zip -rm {output} {config[pathofact][datadir]}/{wildcards.project}/logs/{wildcards.sample}
        rm -rf {config[pathofact][datadir]}/{wildcards.project}/splitted/{wildcards.sample}
        rm -rf {config[pathofact][datadir]}/{wildcards.project}/PathoFact_intermediate/SignalP/splitted/{wildcards.sample}
        find {config[pathofact][datadir]}/{wildcards.project} -type d -empty -delete
        """

rule clean_AMR_workflow:
    input: os.path.join(DATA_DIR,"{project}/PathoFact_report/AMR_MGE_prediction_{sample}_report.tsv")
    output: os.path.join(DATA_DIR,"{project}/logs/AMR_{sample}_compressed.zip")
    shell: """
        zip -rm {output} {config[pathofact][datadir]}/{wildcards.project}/logs/{wildcards.sample}
        rm -rf {config[pathofact][datadir]}/{wildcards.project}/splitted/{wildcards.sample}
        rm -rf {config[pathofact][datadir]}/{wildcards.project}/contig_splitted/{wildcards.sample}_dir
        rm -rf {config[pathofact][datadir]}/{wildcards.project}/PathoFact_intermediate/MGE/plasmid_splitted/{wildcards.sample}
        find {config[pathofact][datadir]}/{wildcards.project} -type d -empty -delete
        """
