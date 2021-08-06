#Prepare fasta

import glob
import os

rule combine_PathoFact:
    input: 
        Virulence_factor= os.path.join(DATA_DIR,"{project}/PathoFact_report/Virulence_prediction_{sample}_report.tsv"),
        Toxins=os.path.join(DATA_DIR,"{project}/PathoFact_report/Toxin_prediction_{sample}_report.tsv"),
        AMR_MGE=os.path.join(DATA_DIR,"{project}/PathoFact_report/AMR_MGE_prediction_{sample}_report.tsv")
    output:
        PathoFact_report= os.path.join(DATA_DIR,"{project}/PathoFact_report/PathoFact_{sample}_predictions.tsv")
    log:
        os.path.join(DATA_DIR,"{project}/logs/{sample}/PathoFact_predictions.log")
    params:
        runtime=config["pathofact"]["runtime"]["medium"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    conda:
        "../../envs/R.yaml"
    message: "Generate PathoFact report on the following sample(s): {wildcards.project} - {wildcards.sample}"
    script:
        "../../scripts/PathoFact.R"
