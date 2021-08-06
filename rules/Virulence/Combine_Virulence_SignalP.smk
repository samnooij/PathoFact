# Combine Virulence with SignalP

import glob
import os

#################################
# Combine Virulence and SignalP #
#################################

rule merge_SignalPVir:
    input:
        hmm=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/HMM_virulence/{sample}_hmm_prediction.tsv"),
        classifier=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/VIRULENCE/classifier_virulence/{sample}_classifier_results_formatted.tsv"),
        SignalP=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/SignalP/aggregated/{sample}_SignalP_results.tsv"),
        translation=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/renamed/{sample}_translation.tsv")
    output:
        Virulence_report=os.path.join(DATA_DIR,"{project}/PathoFact_report/Virulence_prediction_{sample}_report.tsv")
    log:
        os.path.join(DATA_DIR,"{project}/logs/{sample}/combine_virulence_results.log")
    params:
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    conda: "../../envs/R.yaml"
    message: "Generate report on predicted virulence factors on the following sample(s): {wildcards.project} - {wildcards.sample}"
    script: "../../scripts/Virulence.R"

