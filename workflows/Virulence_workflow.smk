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
    '../rules/Universal/Clean_up_individual.smk'

# master command
rule Analysis_Virulence:
    input:
        expand(
            [
                os.path.join(DATA_DIR,"{project}/PathoFact_report/Virulence_prediction_{sample}_report.tsv"),
                os.path.join(DATA_DIR,"{project}/logs/VF_{sample}_compressed.zip")
            ],
            project=config["pathofact"]["project"], sample=config["pathofact"]["sample"]
        )
    output:
        touch('virulence_analyis.done')
