# Changes that have been made to PathoFact

# 1. Update snakemake version

The default version of [`snakemake`](https://snakemake.readthedocs.io/)
used in PathoFact cannot use the `{name}` and `{jobid}` variables to 
communicate with SLURM. I like to use those to create named log files.
Therefore, I updated snakemake to version 5.10.0.
(Any later version would likely work as well.)

# Shell commands used to make updates

```bash
# Create the original PathoFact conda environment
conda env create -f envs/PathoFact.yaml
conda activate PathoFact

# Update snakemake
conda install snakemake=5.10.0 -c bioconda -c conda-forge
# Overwrite the conda environment yaml file
conda env export > envs/PathoFact.yaml
```
