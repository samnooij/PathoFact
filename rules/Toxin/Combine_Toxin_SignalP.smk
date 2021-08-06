#Toxin

import glob
import os

#################################
# Combine Toxin HMM and SignalP #
#################################

# Put Toxin HMM results in the correct format & join SignalP and Toxin HMM files
rule R_script:
    input:
        input_HMM=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/TOXIN/HMM_toxin/{sample}.Input_HMM_R.csv"),
        translation=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/renamed/{sample}_translation.tsv"),
        signalP=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/SignalP/aggregated/{sample}_SignalP_results.tsv"),
        library=config["pathofact"]["tox_lib"]
    output:
        gene_library=os.path.join(DATA_DIR,"{project}/PathoFact_report/Toxin_gene_library_{sample}_report.tsv"),
        gene_toxic=os.path.join(DATA_DIR,"{project}/PathoFact_report/Toxin_prediction_{sample}_report.tsv")
    log:
        os.path.join(DATA_DIR,"{project}/logs/{sample}/gene_table_library.log")
    params:
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    message: "Generate Toxin report on the following sample(s): {wildcards.project} - {wildcards.sample}"
    conda:
        "../../envs/R.yaml"
    script:
        "../../scripts/ownHMM_library.R"

