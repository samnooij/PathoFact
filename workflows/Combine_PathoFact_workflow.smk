# include rules for the Analysis part of the workflow
include:
    '../rules/Universal/Preprocessing.smk'
include:
    '../rules/Universal/SignalP.smk'
include:
    '../rules/Virulence/Virulence.smk'
include:
    '../rules/Virulence/Combine_Virulence_SignalP.smk'
include:
    '../rules/Toxin/Toxin.smk'
include:
    '../rules/Toxin/Combine_Toxin_SignalP.smk'
include:
    '../rules/Universal/Preprocessing_contig.smk'
include:
    '../rules/AMR/AMR.smk'
include:
    '../rules/AMR/Plasmid.smk'
include:
    '../rules/AMR/Phage.smk'
include:
    '../rules/AMR/Combine_MGE_AMR.smk'
include:
    '../rules/Universal/Combine_PathoFact.smk'
include: 
    '../rules/Universal/Clean_up.smk'
# master command
rule Analysis:
    input:
        expand(
            [
                os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/AMR/{sample}_AMR_MGE_prediction_detailed.tsv"),
                os.path.join(DATA_DIR,"{project}/PathoFact_report/Toxin_gene_library_{sample}_report.tsv"),
                os.path.join(DATA_DIR,"{project}/PathoFact_report/PathoFact_{sample}_predictions.tsv"),
                os.path.join(DATA_DIR,"{project}/logs/{sample}_compressed.zip")
            ],
            project=config["pathofact"]["project"], sample=config["pathofact"]["sample"]
        )
    output:
        touch('PathoFact_analyis.done')
