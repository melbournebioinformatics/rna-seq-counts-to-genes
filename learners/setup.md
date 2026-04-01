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

For this workshop, you can use the following code to install the necessary packages:

```r
# CRAN packages
install.packages(c("ggplot2","ggrepel","ggfortify","scales","pheatmap","matrixStats","openxlsx"))

# Bioconductor packages
install.packages("BiocManager")
BiocManager::install(c("limma","edgeR"))

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

