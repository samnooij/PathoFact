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

# 3. Use `slurm-cluster-status`

On the HPC cluster of the Leiden University Medical Center,
snakemake does not always communicate well with SLURM and
using a custom script can help.
(It should not hurt to use this, but not everybody may need it.)

The script is [on GitHub](https://github.com/LUMC/slurm-cluster-status)
and was downloaded to the current directory using the command:

```bash
git clone https://github.com/LUMC/slurm-cluster-status.git
```

The following line was added to the [profile](config/config.yaml) to use it:

```yaml
cluster-status: "slurm-cluster-status/slurm-cluster-status.py"
```

# 4. Save split fasta files in one place

By default, PathoFact splits large input (FASTA) files into 
smaller files to speed up processing.
However, these split files are stored in two different locations,
causing more disk space to be used than is necessary.
To fix this, I changed the path to which these split files are
written and from which they are read so that there is only one
place to put them.
The rules that I changed are:

 - `rules/Universal/Preprocessing.smk`: checkpoint splitting  
 - `rules/AMR/AMR.smk`: run_deepARG, run_RGI  
 - `rules/Toxin/Toxin.smk`: run_HMM_tox  
 - `rules/Virulence/Virulence.smk`: run_HMM_vir, AAC, DPC, CTDC, CTDT, CTDD

Now all these files are stored under:
`{project}/PathoFact_intermediate/splitted/`.

## 4.2. Make split output explicit

Related to the split files, the rules 'splitcontig' from
[`Preprocessing_contig.smk`](rules/Universal/Preprocessing_contig.smk)
and 'splitplasmid' from
[`Plasmid.smk`](rules/AMR/Plasmid.smk) have named variables
in their snakefiles under 'output', but did not use these
names in the shell command.
These have been adjusted to make the output more explicit.  
(This was a manual change, I have no code for this.)

# 5. Reserve extra time for very long jobs

The rules 'run_VirFinder' and 'run_VirSorter'
from [`Phage.smk`](rules/AMR/Phage.smk) may take way longer
than the default 2 hours for large datasets.
These two jobs have been assigned a new class of runtimes
from the configuration file: 'extra_long'.
The time has been set to 48 hours, but can be adjusted if
necessary.

(This has also been changed manually.
Edit the line `    extra_long: "48:00:00"` near the bottom
of [`parameters.yaml`](config/parameters.yaml) to
change this to any time you need.)

# 6. Reserve CPU threads for DeepARG

PathoFact uses [DeepARG](https://bitbucket.org/gusphdproj/deeparg-ss)
to predict antibiotic resistance genes, which in turn uses
[DIAMOND](https://github.com/bbuchfink/diamond).
Now DIAMOND by default uses the maximum number of CPU threads available
and DeepARG has no way of setting the threads manually.
(This is a [known issue](https://bitbucket.org/gusphdproj/deeparg-ss/issues/16/add-option-to-set-diamond-threads).)

So now, by default, PathoFact reserves 1 CPU thread on the
HPC node, the node may have 48 threads available,
DeepARG starts using 48 threads on a system that may also
have other jobs running. This is generally not a good idea
and may get several jobs crashed.

As a workaround, in the meantime, I'm reserving 24 threads
for this rule, because most of the nodes on the local HPC
cluster have 24 threads.

(This was also a manual change, with no code.)

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

# 3. (optional) use slurm-cluster-status script
# Download the script from GitHub
git clone https://github.com/LUMC/slurm-cluster-status.git
# Then add the line to the profile configuration file
echo "cluster-status: \"slurm-cluster-status/slurm-cluster-status.py\"" >> config/config.yaml

# 4. save all split files in the same place
# (this command should change all references of 'splitted' to the same
# directory, although I didn't test this!)
sed -i 's/{project}\/splitted/{project}\/PathoFact_intermediate\/splitted/g' rules/*/*.smk
```
