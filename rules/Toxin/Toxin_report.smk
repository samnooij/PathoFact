#Toxin

import glob
import os

#################################
# Combine Toxin HMM and SignalP #
#################################

# Put Toxin HMM results in the correct format & join SignalP and Toxin HMM files
rule R_script:
    input:
        input_HMM="{datadir}/{project}/TOXIN/HMM_toxin/{sample}.Input_HMM_R.csv",
        translation="{datadir}/{project}/renamed/{sample}_translation.tsv",
        signalP="{datadir}/{project}/SignalP/aggregated/{sample}_SignalP_results.tsv",
        library=config["pathofact"]["tox_lib"]
    output:
        gene_library="{datadir}/{project}/Toxin_gene_library_{sample}_report.tsv",
        gene_toxic="{datadir}/{project}/Toxin_prediction_{sample}_report.tsv"
    log:
        "{datadir}/{project}/logs/{sample}_gene_table_library.log"
    message:
        "Run external R script to join SignalP and ToxinHMM and create Toxin report (incl. confidence levels)"
    params:
        outdir="{datadir}",
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    conda:
        "../../envs/R.yaml"
    script:
        "../../scripts/ownHMM_library.R"

