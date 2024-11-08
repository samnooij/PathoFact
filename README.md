⚠️⚠️⚠️⚠️⚠️

This is just a fork of the original PathoFact repository, which is hosted here: 
[https://git-r3lab.uni.lu/laura.denies/PathoFact](https://git-r3lab.uni.lu/laura.denies/PathoFact)

This repository is not maintained and I cannot help with questions and problems posted in the issues section.
Sorry!

⚠️⚠️⚠️⚠️⚠️

# Alternative PathoFact code

This repository contains updated code for PathoFact to make it work
more conveniently with SLURM on High Performance Computing (HPC)
clusters with big datasets.
To keep the repository small, database files and other big files
are _not included_.
Before use, please download the original PathoFact from
[here](https://gitlab.lcsb.uni.lu/laura.denies/PathoFact).

The files that have been left out are:

 - `databases/*`  
 - `localDB/*`  
 - `analytical_code/VF_*_subset_final*faa`  
 - `scripts/Virulence_factor_model.sav`  
 - `scripts/finalized_model_3f.sav`

For a complete list of changes and short rationale, see
[this document](Changes_made.md).

# PathoFact v1.0 

PathoFact is an easy-to-use modular pipeline for the metagenomic analyses of toxins, virulence factors and antimicrobial resistance. 
Additionally, PathoFact combines the prediction of these pathogenic factors with the identification of mobile genetic elements. 
This provides further depth to the analysis by considering the localization of the genes on mobile genetic elements (MGEs), as well as on the chromosome. 
Furthermore, each module (toxins, virulence factors, and antimicrobial resistance) of PathoFact is also a standalone component, making it a flexible and versatile tool. 

For further information regarding usage and generated reports, please see the [Documentation](https://git-r3lab.uni.lu/laura.denies/PathoFact/-/wikis/home)

# Requirements and installation

The main requirements are:
- `gcc/g++`
- `git` and `git lfs`
- `conda`

Most other dependencies are either included as a `git submodule` or will be installed automatically by `snakemake` using the `conda` YAML files in `envs/`.
However, some tools need to be installed and/or configured manually.


## Miniconda (conda)

```bash
# install miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
chmod u+x Miniconda3-latest-Linux-x86_64.sh
./Miniconda3-latest-Linux-x86_64.sh # follow the instructions
```

## PathoFact

Clone the repository and its sub-modules:

```bash
# activate git lfs
git lfs install
# clone branch incl. sub-modules
git clone -b master --recursive https://git-r3lab.uni.lu/laura.denies/PathoFact.git
```

## Pipeline environment

```bash
# create the conda environment
conda env create -f=envs/PathoFact.yaml
```

You can activate and deactivate the environment using `conda activate PathoFact` and `conda deactivate`, respectively. PathoFact requires a `snakemake version >= 5.5.4`.

## Dependency: `SignalP`

Required version: `5.0`

To download the tool you need to submit a request at https://services.healthtech.dtu.dk/:

- Look for "SignalP" and click on the link
- Click on "Downloads"
- Click on the link for your platform for the version "Version 5.0g"
- Fill and submit the form

After the installation, adjust the path for this tool in `config.yaml` and `test/test_config.yaml` (keyword `signalp`).

# Usage

## Input files

Each sample should have one input file:

- `*.fna`: FASTA file containing nucleotide sequences of the contigs
    - no whitespaces in FASTA headers (i.e. to remove whitespaces: sed -i '/^>/ s/ .*//' file.fna)
    - for prediction of mobile genetic elements and input for prodigal

The following files are generated by PathoFact itself:

- `*.faa`: FASTA file conatining translated gene sequences, i.e. amino acid sequences
    - no whitespaces in FASTA headers
    - for prediction of toxins, virulence factors and antimicrobial resistance genes
- `*.contig`: TAB-delimited file containing a mapping from contig ID (1st column) to gene ID (2nd column)
    - no header, one gene ID per line
    - contig and gene IDs should be the same as in the FASTA files

The input file for each sample should be located in the same directory.
For each sample, the corresponding input files should have the same basename, e.g. `SAMPLE_A.fna` for sample `SAMPLE_A`.

**NOTE**: For preprocessing and assembly of metagenomic reads we would suggest using IMP (https://imp.pages.uni.lu/web/)

## Run PathoFact

### Configuration

To run PathoFact you need to adjust some parameters in `config.yaml`.

- `sample`: This is a list of sample names, e.g.

```yaml
  sample:
    - Sample1
    - Sample2
```

(Where 'Sample1' and 'Sample2' match the file names
`Sample1.fna` and `Sample2.fna` in the input directory.)

- `project`: A unique project name which will be used as the name of the output directory in `datapath` path (see below).
- `datadir`: Path to directory containing the sample data; the output directory will be created there.
- `workflow`: Pathofact can run the complete pipeline (default) or a specific step:
    - "complete": complete pipeline = toxin + virulence + AMR + MGE prediction
    - "Tox": toxin prediction
    - "Vir": virulence prediction
    - "AMR": antimicrobial resistance (AMR) & mobile genetic elements (MGEs) prediction

### Execution

Basic command to run the pipeline on SLURM:

```bash
# activate the env
conda activate PathoFact
# run the pipeline
snakemake --profile config
```

**NOTE**: Add parameter `-n` (or `--dry-run`) to the command to see which steps will be executed without running them.

**NOTE**: Add `--configfile <configfile.yaml>` to use a different config file than `config.yaml`. 

For more options, see the [snakemake documentation](https://snakemake.readthedocs.io/en/stable/index.html).  
The default settings are stored in [config/config.yaml](config/config.yaml).


### Test module

To test for the correct installation of the pipeline the testmodule can be run:

Include the required path to the SingalP v5.0 installation to the config file `test/test_config.yaml`

```
# activate env
conda activate PathoFact
# run the pipeline
snakemake -s test/Snakefile --use-conda --reason --cores 1 -p
```
