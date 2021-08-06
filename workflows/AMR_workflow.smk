# include rules for the Analysis part of the workflow
include:
    '../rules/Universal/Preprocessing.smk'
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
    '../rules/Universal/Clean_up_individual.smk'

# master command
rule AMR_Analysis:
    input: 
        expand(
            [
                os.path.join(DATA_DIR,"{project}/PathoFact_report/AMR_MGE_prediction_{sample}_report.tsv"),        
                os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/AMR/{sample}_AMR_MGE_prediction_detailed.tsv"),
                os.path.join(DATA_DIR,"{project}/logs/AMR_{sample}_compressed.zip")
           ],
            project=config["pathofact"]["project"], sample=config["pathofact"]["sample"]
        )
    output:
        touch('AMR_analyis.done')
