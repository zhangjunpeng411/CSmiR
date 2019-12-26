# CSmiR
Exploring cell-specific miRNA regulation from single-cell miRNA-mRNA co-sequencing data

## Background
Single-cell miRNA-mRNA co-sequencing opens up an insight to investigate miRNA regulation at single-cell level. However, existing methods are only applicable to study miRNA regulation from a group of cells rather than a single cell. Consequently, the unique miRNA regulation in each single cell will be annihilated. To understand the heterogeneity of miRNA regulation in each single cell, there is a strong need of novel methods to explore cell-specific miRNA regulation at a single-cell resolution level. In this work, we propose a novel method to explore Cell-Specific miRNA regulation (CSmiR) for each single cell from single-cell miRNA-mRNA co-sequencing data. To the best of our knowledge, CSmiR is the first method to explore miRNA regulation at a single-cell resolution level, and we believe that it can be applied to uncover miRNA regulation in other single cells.

## Description of each file
Exp_K562_19_single_cells.RData: Matched miRNA and mRNA expression data across 19 half K562 cells, Putative miRNA-target binding information from miRTarBase v8.0 and TarBase v8.0, the list of CML-related miRNAs and mRNAs from HMDD v3.0 and DisGeNET v5.0.

Conserved and rewired results.md: The results of conserved and rewired miRNA-mRNA regulatory networks and hub miRNAs across 19 half K562 cells.

CSmiR.R: Scripts for exploring cell-specific miRNA regulation.

## The usage of CSmiR
Paste all files into a single folder (set the folder as the directory of R environment), the workflow of CSmiR is implemented in CSmiR.R. The users can simply run the scripts as follows.

```{r echo=FALSE, results='hide', message=FALSE}
source("CSmiR.R")
```
