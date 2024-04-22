# TODO: Replace with the name of the repo

[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/miniconda/en/latest/)

Note: Analysis repo names should be prefixed with the year (ie `2024-noveltree-analysis`). This prefix can be changed at time of publication if appropriate.

## Purpose

TODO: Briefly describe the core analyses performed in the repository and the motivation behind them.

## Installation and Setup

This repository uses conda to manage software environments and installations. You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/projects/miniconda/en/latest/). After installing conda and [mamba](https://mamba.readthedocs.io/en/latest/), run the following command to create the pipeline run environment.

```{bash}
TODO: Replace <NAME> with the name of your environment
mamba env create -n <NAME> --file envs/dev.yml
conda activate <NAME>
```

**Tips for Developers**

You can use the following command to export your current conda environment to a `yml` file.  
This command will only export the packages that you have installed directly, not the ones that were installed as dependencies. When you're ready to share, please delete this section.

```{bash}
conda env export --from-history --no-builds > envs/dev.yml
```

## Data

TODO: Add details about the description of input / output data and links to Zenodo depositions, if applicable.

## Overview

### Description of the folder structure

### Methods

TODO: Include a brief, step-wise overview of analyses performed.

> Example:
>
> 1. Download scripts using `download.sh`.
> 2. Preprocess using `./preprocessing.sh -a data/`
> 3. Run analysis script using `analysis.Rscript`
> 4. Generate figures using `pub/make_figures.R`.

### Compute Specifications

TODO: Describe what compute resources were used to develop and run the analysis. For example, you could list the operating system, number of cores, RAM, and storage space. You should log any major changes to the compute specifications here as they happen.

## Contributing

See how we recognize [feedback and contributions to our code](https://github.com/Arcadia-Science/arcadia-software-handbook/blob/main/guides-and-standards/guide-credit-for-contributions.md).

---
## For Developers

This section contains information for developers who are working off of this template. Please delete this section when you're ready to share your repository.

### GitHub templates
This template uses GitHub templates to provide checklists when making new pull requests. These templates are stored in the [.github/](./.github/) directory.

### `.gitignore`
This template uses a `.gitignore` file to prevent certain files from being committed to the repository.

### Linting
This template automates linting using GitHub Actions and the [`lintr` linter](https://cran.r-project.org/web/packages/lintr/vignettes/lintr.html). When you push changes to your repository, GitHub will automatically run the linter and report any errors, blocking merges until they are resolved. 
