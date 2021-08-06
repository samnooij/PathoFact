# Changes that have been made to PathoFact

# 1. Update snakemake version

The default version of [`snakemake`](https://snakemake.readthedocs.io/)
used in PathoFact cannot use the `{name}` and `{jobid}` variables to 
communicate with SLURM. I like to use those to create named log files.
Therefore, I updated snakemake to version 5.10.0.
(Any later version would likely work as well.)
