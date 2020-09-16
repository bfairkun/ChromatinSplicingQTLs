# Snakemake workflow: ChromatinSplicing

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.3.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/ChromatinSplicingQTLs.svg?branch=master)](https://travis-ci.org/snakemake-workflows/ChromatinSplicingQTLs)


## Authors

* Benjamin Fair (@bfairkun)

## Usage

### Step 1: Install workflow

If you simply want to use this workflow, clone the [latest release](https://github.com/bfairkun/ChromatinSplicingQTLs). This workflow is all contained as a snakemake workflow in the `code` directory.
If you intend to modify and further develop this workflow, fork this repository. Please consider providing any generally applicable modifications via a pull request.

### Step 2: Install dependencies

Once in the code directory, install the dependencies as a conda environment using the `./envs/MAIN_ENV.yaml` file:

    conda env create -f envs/MAIN_ENV.yaml

If solving the environment somehow fails, I recommend installing [mamba](https://github.com/mamba-org/mamba) and using `mamba` as a drop-in replacement for `conda` in the above command. Once created, activate the environment.

    conda activate ChromatinSplicingQTLs_env

Because some rules in the snakemake require their own isolated conda environments, these additionally need to be created by snakemake:

    snakemake --use-conda --create-envs-only


### Step 3: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml`. Configure cluster settings in `cluster-config.json`

### Step 4: Execute workflow

Test your configuration by performing a dry-run via

    snakemake -n

Execute the workflow locally via

    snakemake --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --cluster --cluster-config cluster-config.json --cluster "sbatch --partition={cluster.partition} --job-name={cluster.name} --output=/dev/null --job-name={cluster.name} --nodes={cluster.n} --mem={cluster.mem}"

or by executing the included sbatch script to execute the snakemake process from a cluster

    sbatch snakemake.sbatch

See the [Snakemake documentation](https://snakemake.readthedocs.io) for further details.
