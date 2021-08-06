# include rules for the Analysis part of the workflow
include:
    '../rules/Universal/Preprocessing.smk'
include:
    '../rules/Universal/SignalP.smk'
include:
    '../rules/Toxin/Toxin.smk'
include:
    '../rules/Toxin/Combine_Toxin_SignalP.smk'
include:
    '../rules/Universal/Clean_up_individual.smk'
# master command
rule Analysis:
    input:
        expand(
            [
               os.path.join(DATA_DIR, "{project}/PathoFact_report/Toxin_gene_library_{sample}_report.tsv"),
               os.path.join(DATA_DIR,"{project}/PathoFact_report/Toxin_prediction_{sample}_report.tsv"),
               os.path.join(DATA_DIR,"{project}/logs/Tox_{sample}_compressed.zip")
            ],
            DATA_DIR=config["pathofact"]["datadir"], project=config["pathofact"]["project"], sample=config["pathofact"]["sample"]
        )
    output:
        touch('Toxin_analyis.done')
