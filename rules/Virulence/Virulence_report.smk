#Virulence

import glob
import os

####################################
#         Final Report             #
####################################

# create final report by combining all files

rule Virulence_merge_final:
    input:
        translation="{datadir}/{project}/renamed/{sample}_translation.tsv",
        prediction="{datadir}/{project}/VIRULENCE/Virulence_prediction/{sample}_Confidence_Virulence_ensembled.csv"
    output:
        temp("{datadir}/{project}/VIRULENCE/Virulence_prediction/{sample}_Virulence_translation_ensembled.csv")
    params:
        outdir="{datadir}",
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    shell:
        """
        join -t $'\\t' <(sort {input.translation}) <(sort {input.prediction}) > {output}
        """

rule virulence_report:
    input:
        "{datadir}/{project}/VIRULENCE/Virulence_prediction/{sample}_Virulence_translation_ensembled.csv"
    output:
        "{datadir}/{project}/Virulence_prediction_{sample}_report.csv"
    params:
        outdir="{datadir}",
        runtime=config["pathofact"]["runtime"]["short"],
        mem=config["pathofact"]["mem"]["normal_mem_per_core_gb"]
    shell:
        "sed -i $'1 i\\\ Sequence no.\\tSequence Query\\tHMM prediction\\tclassifier prediction\\tVirulence_prediction\\tSignalP\\tConfidence level' {input};"
        "cp {input} {output}"
