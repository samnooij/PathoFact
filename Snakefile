#Snakefile

configfile: "config.yaml"
DATA_DIR=config["pathofact"]["datadir"]

if config["pathofact"]["workflow"] == "complete":
    include:
        "workflows/Combine_PathoFact_workflow.smk"
    rule all:
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
elif config["pathofact"]["workflow"] == "Tox":
    include:
        "workflows/Toxin_workflow.smk"
    rule all:
        input:
            expand(
                [
                    os.path.join(DATA_DIR,"{project}/PathoFact_report/Toxin_prediction_{sample}_report.tsv"),
                    os.path.join(DATA_DIR,"{project}/PathoFact_report/Toxin_gene_library_{sample}_report.tsv"),
                    os.path.join(DATA_DIR,"{project}/logs/Tox_{sample}_compressed.zip")
                ],
                project=config["pathofact"]["project"], sample=config["pathofact"]["sample"]
            )
elif config["pathofact"]["workflow"] == "Vir":
    include:
        "workflows/Virulence_workflow.smk"
    rule all:
        input:
            expand(
                [
                    os.path.join(DATA_DIR,"{project}/PathoFact_report/Virulence_prediction_{sample}_report.tsv"),
                    os.path.join(DATA_DIR,"{project}/logs/VF_{sample}_compressed.zip")
                ],
                project=config["pathofact"]["project"], sample=config["pathofact"]["sample"]
            )                    
elif config["pathofact"]["workflow"] == "AMR":
    include:
        "workflows/AMR_workflow.smk"
    rule all:
        input:
            expand(
                [
                    os.path.join(DATA_DIR,"{project}/PathoFact_report/AMR_MGE_prediction_{sample}_report.tsv"),
                    os.path.join(DATA_DIR,"{project}/PathoFact_intermediate/AMR/{sample}_AMR_MGE_prediction_detailed.tsv"),
                    os.path.join(DATA_DIR,"{project}/logs/AMR_{sample}_compressed.zip")
                ],
                project=config["pathofact"]["project"], sample=config["pathofact"]["sample"]
            )
else:
    raise Exception("Unknown workflow option: %s" % config["pathofact"]["workflow"])
