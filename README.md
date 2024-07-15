# 2024-protein-universe

[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/miniconda/en/latest/)

## Purpose

This repository contains code related to the pub "The known protein universe is phylogenetically biased".

## Installation and Setup

This repository uses conda to manage software environments and installations. You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/projects/miniconda/en/latest/). After installing conda and [mamba](https://mamba.readthedocs.io/en/latest/), run the following command to create the pipeline run environment.

```{bash}
mamba env create -n protein_universe --file envs/dev.yml
conda activate protein_universe
```

## Data



## Overview

### Description of the folder structure

The repository is organized into the following top-level directories.

* **code**:
  R scripts used for downloading and cleaning data, performing analysis, and generating figures presented in the pub.
* **data**:
 .RDS files used in analyses.
* **envs**:
  YAML file including the packages and dependencies used for creating the conda environment.

```
─ code
  ├── README.md
  ├── protein-universe-analysis.R
  ├── protein-universe-data.R
  └── protein-universe-utils.R
─ data
  ├── README.md
  ├── afdb_cluster_stats.RDS
  ├── afdb_cluster_taxonomy.RDS
  ├── afdb_genome_size_stats.RDS
  ├── pdb_metadata.RDS
  ├── pdb_taxonomy.RDS
  ├── timetree_phylogeny_cleaned.RDS
  └── timetree_taxonomy.RDS
─ envs
  └── dev.yml
```

### Methods

> 1. Download, clean, and organize data using `protein-universe-data.R`.
> 2. Load supporting functions using `protein-universe-utils.R`
> 3. Run analyses using `protein-universe-analysis.R`

### Compute Specifications

All analyses were done on an Apple MacBook Pro running macOS Montery with 32GB RAM, 10 cores, and 1TB of storage.

## Contributing

See how we recognize [feedback and contributions to our code](https://github.com/Arcadia-Science/arcadia-software-handbook/blob/main/guides-and-standards/guide-credit-for-contributions.md).
