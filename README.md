# CSmiR
Exploring cell-specific miRNA regulation with single-cell miRNA-mRNA co-sequencing data

## Background
Single-cell miRNA-mRNA co-sequencing opens up an insight to investigate miRNA regulation at single-cell level. However, existing methods are only applicable to study miRNA regulation from a group of cells rather than a single cell. Consequently, the unique miRNA regulation in each single cell will be annihilated. To understand the heterogeneity of miRNA regulation in each single cell, there is a strong need of novel methods to explore cell-specific miRNA regulation at a single-cell resolution level. In this work, we propose a novel method to explore Cell-Specific miRNA regulation (CSmiR) for each single cell by integrating single-cell miRNA-mRNA co-sequencing data and putative miRNA-mRNA target binding information. To the best of our knowledge, CSmiR is the first method to explore miRNA regulation at a single-cell resolution level, and we believe that it can be applied to uncover miRNA regulation in other single cells.

## Description of each file
Exp_K562_19_single_cells.RData: Matched miRNA and mRNA expression data across 19 half K562 cells, Putative miRNA-target binding information from TargetScan, Validated miRNA-target binding information from miRTarBase v8.0 and TarBase v8.0, the list of CML-related miRNAs and mRNAs from HMDD v3.0 and DisGeNET v5.0.

CSmiR.R: Utility functions for exploring cell-specific miRNA regulation.

CSmiRsyn.R: Utility functions for exploring cell-specific miRNA synergism.

Case_study.R: Running scripts for exploring cell-specific miRNA regulation.

## The usage of CSmiR
Paste all files into a single folder (set the folder as the directory of R environment), the workflow of CSmiR is implemented in CSmiR.R. The users can simply run the scripts as follows.

```{r echo=FALSE, results='hide', message=FALSE}
source("Case_study.R")
```
## Quick example to use CSmiR
For identifying cell-specific miRNA regulation, users should prepare matched miRNA and mRNA single-cell co-expression data. Paste the datasets and our source file (CSmiR.R) into a single folder (set the folder as the directory of R environment), users can use the following scripts to identify cell-specific miRNA regulation. For convenience, the datasets prepared for users are from our datasets (Exp_K562_19_single_cells.RData).

```{r echo=FALSE, results='hide', message=FALSE}
## Load required R packages, please firstly install the following R packages before running scripts
library(pracma)
library(WGCNA)
library(igraph)
library(miRspongeR)
library(biclique)
library(corrplot)
library(dendextend)
library(ggplot2)
library(ggsignif)
library(cowplot)
library(clusterProfiler)
library(msigdbr)

## Load utility functions
source("CSmiR.R")

## Load prepared datasets
load("Exp_K562_19_single_cells.RData")

## Preprocess the single-cell sequencing data including log2(x+1), compute the average expression values of duplicate genes
## and remove genes with constant expression values in all cells

    # Transformation using log2(x+1)
    miRNA_scRNA_norm <- log2(miRNA_scRNA_raw+1)
    mRNA_scRNA_norm <- log2(mRNA_scRNA_raw+1)

    # Compute the average expression values of duplicate genes
    miRNA_scRNA_norm_average <- Averg_Duplicate(miRNA_scRNA_norm)
    mRNA_scRNA_norm_average <- Averg_Duplicate(mRNA_scRNA_norm)

    # Remove genes with constant expression values in all cells
    miRNA_scRNA_norm_sd <- unlist(lapply(seq(dim(miRNA_scRNA_norm_average)[2]), function(i) sd(miRNA_scRNA_norm_average[, i])))
    miRNA_scRNA_norm_filter <- miRNA_scRNA_norm_average[, which(miRNA_scRNA_norm_sd > 0)]
    mRNA_scRNA_norm_sd <- unlist(lapply(seq(dim(mRNA_scRNA_norm_average)[2]), function(i) sd(mRNA_scRNA_norm_average[, i])))
    mRNA_scRNA_norm_filter <- mRNA_scRNA_norm_average[, which(mRNA_scRNA_norm_sd > 0)]

## Discovering cell-specific miRNA-mRNA regulatory network. The running time depends on the number of bootstrapping (100 in default).
## Users can set a small number of bootstrapping (i.e. bootstrap_num = 3) to run the following script.
    CSmiR_network_bootstrap <- CSmiR_net_bootstrap(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, 
                                                   boxsize = 0.1, bootstrap_betw_point = 5, 
						   bootstrap_num = 100, p.value.cutoff = 0.01)
```




