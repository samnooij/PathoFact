# Changes that have been made to PathoFact

(Shell commands used to make some of the changes are listed
at [the bottom](#shell-commands-used-to-make-updates).)

# 1. Update snakemake version

The default version of [`snakemake`](https://snakemake.readthedocs.io/)
used in PathoFact cannot use the `{name}` and `{jobid}` variables to 
communicate with SLURM. I like to use those to create named log files.
Therefore, I updated snakemake to version 5.10.0.
(Any later version would likely work as well.)

# 2. Control snakemake through `--profile` and move config to a directory

I prefer to use snakemake using the `--profile` parameter and a text file
that controls all the optional command-line parameters.
Besides, snakemake seems to work best with a separate file for these
command-line arguments.
Therefore, the parameters related to the samples used within the
workflow are stored as a separate file and these two configuration files
are stored under [`config/`](config).

The [profile file](config/config.yaml) reads as follows:

```yaml
jobs: 50
jobname: PathoFact-{name}.{jobid}
use-conda: true
latency-wait: 60
keep-going: true
cluster: "sbatch --parsable -N 1 -n 1 -c {threads} --mem={params.mem} -t {params.runtime} -D . -e log/{name}-{jobid}.err -o log/{name}-{jobid}.out"
```

And this profile can be used with:

```bash
snakemake --profile config
```

(Any command-line parameters can still be appended to this command, so
for example add a `-n` to first try a dry-run.)

Furthermore, since the newly added `cluster` argument was added,
all snakemake rules that are sent as jobs to scheduler (SLURM)
need to have their required memory (RAM) and runtime specified as
`params` within the snakefiles. The rules that did not have this
already are:

 - `rules/Universal/Preprocessing.smk`: Prodigal  
 - `rules/Universal/Preprocessing.smk`: mapping_file  
 - `rules/Universal/Clean_up_individual.smk`: all rules  
 - `rules/Universal/Clean_up.smk`: clean_all

So these lines have been added.

# Shell commands used to make updates

```bash
# Create the original PathoFact conda environment
conda env create -f envs/PathoFact.yaml
conda activate PathoFact

# 1. Update snakemake
conda install snakemake=5.10.0 -c bioconda -c conda-forge
# Overwrite the conda environment yaml file
conda env export > envs/PathoFact.yaml

# 2. Reorganise configuration files
mkdir config
mv config.yaml config/parameters.yaml
touch config/config.yaml
# (I edited both using a text editor program.)

# the 'log' directory is necessary to work with the profile config
mkdir log
```
