# Snakemake workflow: ChromatinSplicing

[![Snakemake](https://img.shields.io/badge/snakemake-≥5.3.0-brightgreen.svg)](https://snakemake.bitbucket.io)
[![Build Status](https://travis-ci.org/snakemake-workflows/ChromatinSplicingQTLs.svg?branch=master)](https://travis-ci.org/snakemake-workflows/ChromatinSplicingQTLs)


## Authors

* Benjamin Fair (@bfairkun)

## Usage

### Step 1: Install workflow

If you simply want to use this workflow, clone the [latest release](https://github.com/bfairkun/ChromatinSplicingQTLs).
If you intend to modify and further develop this workflow, fork this repository. Please consider providing any generally applicable modifications via a pull request.

### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml`. Configure cluster settings in `cluster-config.json`

### Step 3: Execute workflow

Test your configuration by performing a dry-run via

    snakemake -n

Execute the workflow locally via

    snakemake --cores $N

using `$N` cores or run it in a cluster environment via

    snakemake --cluster --cluster-config cluster-config.json --cluster "sbatch --partition={cluster.partition} --job-name={cluster.name} --output=/dev/null --job-name={cluster.name} --nodes={cluster.n} --mem={cluster.mem}"

or by executing the included sbatch script to execute the snakemake process from a cluster

    sbatch snakemake.sbatch

See the [Snakemake documentation](https://snakemake.readthedocs.io) for further details.
