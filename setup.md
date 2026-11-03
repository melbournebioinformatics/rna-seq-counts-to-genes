---
title: Setup
---

Please follow the steps below and install the required software **before** the scheduled workshop.


## Data Sets

**Data:** [Obers et al. "Retinoic acid and TGF-β orchestrate organ-specific programs of tissue residency", *Immunity* (2024)](https://www.cell.com/immunity/abstract/S1074-7613(24)00459-X)

Download the [GSE232852_CountsTable.txt data file](episodes/data/GSE232852_CountsTable.txt) and place it in the workshop folder.

## RStudio Setup

We use RStudio for coding in R.

[Click here and follow the instructions](https://posit.co/download/rstudio-desktop/) to install RStudio Desktop in your system.

::::::::::::::::: discussion

### R packages

Most workshops using R will require the installation of specific packages. Make sure to check in advance with the workshop organisers what packages need to be installed. 

You can install packages from CRAN using:

```r
install.packages("package_name")
```

If your package is in a different R repository, such as Bioconductor or GitHub, you may need the [BiocManager](https://www.bioconductor.org/install/) or [devtools](https://devtools.r-lib.org/) packages to install them. To install BiocManager:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
```

For devtools, you can simply do:

```r
install.packages("devtools")
```

You can then install packages directly from GitHub with:

```r
devtools::install_github("username/reponame")
```

::::::::::::::::::::::::::::


<!--
READ HERE FOR OS-SPECIFIC INSTRUCTIONS

Setup for different systems can be presented in dropdown menus via a `spoiler`
tag. They will join to this discussion block, so you can give a general overview
of the software used in this lesson here and fill out the individual operating
systems (and potentially add more, e.g. online setup) in the solutions blocks.
-->

<!--
:::::::::::::::: spoiler

### Windows

Use PuTTY

::::::::::::::::::::::::

:::::::::::::::: spoiler

### MacOS

Use Terminal.app

::::::::::::::::::::::::


:::::::::::::::: spoiler

### Linux

Use Terminal

::::::::::::::::::::::::
-->

