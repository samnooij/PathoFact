#AMR

import glob
import os

##########################
#     AMR Prediction     #
##########################

rule combine_AMR_plasmid:
    input:
        ORF_translation=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/renamed/{sample}_translation.tsv"),
        Contig_ORFs=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/Prodigal/{sample}.contig"),
        Contig_translation=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/renamed/{sample}_Contig_translation.tsv"),
        AMR=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/AMR/{sample}_AMR_prediction.tsv"),
        PlasFlow=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/MGE/plasmid/PlasFlow/{sample}_plasflow_prediction_final.tsv"),
        DeepVirFinder=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/MGE/phage/{sample}_VirFinder_aggregated.csv"),
        VirSorter=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/MGE/phage/{sample}_VIRSorter_aggregated.csv")
    output:    
        Report_1=os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/AMR/{sample}_AMR_MGE_prediction_detailed.tsv"),
        Report_2=os.path.join(DATA_DIR,"{project}/PathoFact_report/AMR_MGE_prediction_{sample}_report.tsv")
    log:
        os.path.join(DATA_DIR,"{project}/logs/{sample}/MGE_AMR_prediction_detail_temp.log")
    params:
        runtime=config["pathofact"]["runtime"]["medium"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    conda:
        "../../envs/R.yaml"
    message: "Merge AMR and MGE predictions for the following sample: {wildcards.project} - {wildcards.sample}"
    script:
        "../../scripts/AMR_MGE.R"

