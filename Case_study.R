##########################################################################################
############## Running scripts for exploring cell-specific miRNA regulation ##############
##########################################################################################

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
library(pheatmap)

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

## Discovering cell-specific miRNA-mRNA regulatory network   
    CSmiR_network_constant <- CSmiR_net(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, 
                                        boxsize = 0.1, interp_betw_point = 5, 
			                interp_type = "constant", p.value.cutoff = 0.01)

    CSmiR_network_linear <- CSmiR_net(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, 
                                        boxsize = 0.1, interp_betw_point = 5, 
			                interp_type = "linear", p.value.cutoff = 0.01)

    CSmiR_network_nearest <- CSmiR_net(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, 
                                        boxsize = 0.1, interp_betw_point = 5, 
			                interp_type = "nearest", p.value.cutoff = 0.01)

    CSmiR_network_spline <- CSmiR_net(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, 
                           boxsize = 0.1, interp_betw_point = 5, 
			   interp_type = "spline", p.value.cutoff = 0.01)

    CSmiR_network_cubic <- CSmiR_net(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, 
                                        boxsize = 0.1, interp_betw_point = 5, 
			                interp_type = "cubic", p.value.cutoff = 0.01)
  
## Experimentally validated cell-specific miRNA-mRNA interactions, the ground-truth (miRTarget variable) is from miRTarBase v8.0 and TarBase v8.0
    miRTarget_graph <- make_graph(c(t(miRTarget[, 1:2])), directed = TRUE)
    
    CSmiR_network_constant_graph <- lapply(seq(CSmiR_network_constant), function(i) make_graph(c(t(CSmiR_network_constant[[i]][, 1:2])), directed = TRUE))
    CSmiR_network_constant_validated <- lapply(seq(CSmiR_network_constant), function(i) as_data_frame(CSmiR_network_constant_graph[[i]] %s% miRTarget_graph))

    CSmiR_network_linear_graph <- lapply(seq(CSmiR_network_linear), function(i) make_graph(c(t(CSmiR_network_linear[[i]][, 1:2])), directed = TRUE))
    CSmiR_network_linear_validated <- lapply(seq(CSmiR_network_linear), function(i) as_data_frame(CSmiR_network_linear_graph[[i]] %s% miRTarget_graph))

    CSmiR_network_nearest_graph <- lapply(seq(CSmiR_network_nearest), function(i) make_graph(c(t(CSmiR_network_nearest[[i]][, 1:2])), directed = TRUE))
    CSmiR_network_nearest_validated <- lapply(seq(CSmiR_network_nearest), function(i) as_data_frame(CSmiR_network_nearest_graph[[i]] %s% miRTarget_graph))

    CSmiR_network_spline_graph <- lapply(seq(CSmiR_network_spline), function(i) make_graph(c(t(CSmiR_network_spline[[i]][, 1:2])), directed = TRUE))    
    CSmiR_network_spline_validated <- lapply(seq(CSmiR_network_spline), function(i) as_data_frame(CSmiR_network_spline_graph[[i]] %s% miRTarget_graph))

    CSmiR_network_cubic_graph <- lapply(seq(CSmiR_network_cubic), function(i) make_graph(c(t(CSmiR_network_cubic[[i]][, 1:2])), directed = TRUE))
    CSmiR_network_cubic_validated <- lapply(seq(CSmiR_network_cubic), function(i) as_data_frame(CSmiR_network_cubic_graph[[i]] %s% miRTarget_graph))

## CML-related cell-specific miRNA-mRNA interactions. the list of CML-related miRNAs and mRNAs (CML variable) is from HMDD v3.0 and DisGeNET v5.0, respectively    
    CSmiR_network_constant_CML <- lapply(seq(CSmiR_network_constant), function(i) CSmiR_network_constant[[i]][intersect(which(CSmiR_network_constant[[i]][, 1] %in% as.matrix(CML)), which(CSmiR_network_constant[[i]][, 2] %in% as.matrix(CML))), ])

    CSmiR_network_linear_CML <- lapply(seq(CSmiR_network_linear), function(i) CSmiR_network_linear[[i]][intersect(which(CSmiR_network_linear[[i]][, 1] %in% as.matrix(CML)), which(CSmiR_network_linear[[i]][, 2] %in% as.matrix(CML))), ])
    
    CSmiR_network_nearest_CML <- lapply(seq(CSmiR_network_nearest), function(i) CSmiR_network_nearest[[i]][intersect(which(CSmiR_network_nearest[[i]][, 1] %in% as.matrix(CML)), which(CSmiR_network_nearest[[i]][, 2] %in% as.matrix(CML))), ])

    CSmiR_network_spline_CML <- lapply(seq(CSmiR_network_spline), function(i) CSmiR_network_spline[[i]][intersect(which(CSmiR_network_spline[[i]][, 1] %in% as.matrix(CML)), which(CSmiR_network_spline[[i]][, 2] %in% as.matrix(CML))), ])
   
    CSmiR_network_cubic_CML <- lapply(seq(CSmiR_network_cubic), function(i) CSmiR_network_cubic[[i]][intersect(which(CSmiR_network_cubic[[i]][, 1] %in% as.matrix(CML)), which(CSmiR_network_cubic[[i]][, 2] %in% as.matrix(CML))), ])

## Overlap of cell-specific miRNA-mRNA interactions across cells    
    Overlap_network_constant <- Overlap_net(CSmiR_network_constant, Intersect_num = round(length(CSmiR_network_constant)*0.9))
    Overlap_network_constant_union <- Overlap_net(CSmiR_network_constant, Intersect_num = 1)
    Overlap_network_constant_two <- Overlap_net(CSmiR_network_constant, Intersect_num = 2)
    Overlap_network_constant_union_graph <- make_graph(c(t(Overlap_network_constant_union[, 1:2])), directed = TRUE)
    Overlap_network_constant_two_graph <- make_graph(c(t(Overlap_network_constant_two[, 1:2])), directed = TRUE)
    Overlap_network_constant_rewired <- as_data_frame(Overlap_network_constant_union_graph %m% Overlap_network_constant_two_graph)
    Overlap_network_constant_graph <- make_graph(c(t(Overlap_network_constant[, 1:2])), directed = TRUE)
    Overlap_network_constant_conserved_validated <- as_data_frame(Overlap_network_constant_graph %s% miRTarget_graph)
    Overlap_network_constant_rewired_graph <- make_graph(c(t(Overlap_network_constant_rewired[, 1:2])), directed = TRUE)
    Overlap_network_constant_rewired_validated <- as_data_frame(Overlap_network_constant_rewired_graph %s% miRTarget_graph)

    Overlap_network_linear <- Overlap_net(CSmiR_network_linear, Intersect_num = round(length(CSmiR_network_linear)*0.9))
    Overlap_network_linear_union <- Overlap_net(CSmiR_network_linear, Intersect_num = 1)
    Overlap_network_linear_two <- Overlap_net(CSmiR_network_linear, Intersect_num = 2)
    Overlap_network_linear_union_graph <- make_graph(c(t(Overlap_network_linear_union[, 1:2])), directed = TRUE)
    Overlap_network_linear_two_graph <- make_graph(c(t(Overlap_network_linear_two[, 1:2])), directed = TRUE)
    Overlap_network_linear_rewired <- as_data_frame(Overlap_network_linear_union_graph %m% Overlap_network_linear_two_graph)
    Overlap_network_linear_graph <- make_graph(c(t(Overlap_network_linear[, 1:2])), directed = TRUE)
    Overlap_network_linear_conserved_validated <- as_data_frame(Overlap_network_linear_graph %s% miRTarget_graph)
    Overlap_network_linear_rewired_graph <- make_graph(c(t(Overlap_network_linear_rewired[, 1:2])), directed = TRUE)
    Overlap_network_linear_rewired_validated <- as_data_frame(Overlap_network_linear_rewired_graph %s% miRTarget_graph)

    Overlap_network_nearest <- Overlap_net(CSmiR_network_nearest, Intersect_num = round(length(CSmiR_network_nearest)*0.9))
    Overlap_network_nearest_union <- Overlap_net(CSmiR_network_nearest, Intersect_num = 1)
    Overlap_network_nearest_two <- Overlap_net(CSmiR_network_nearest, Intersect_num = 2)
    Overlap_network_nearest_union_graph <- make_graph(c(t(Overlap_network_nearest_union[, 1:2])), directed = TRUE)
    Overlap_network_nearest_two_graph <- make_graph(c(t(Overlap_network_nearest_two[, 1:2])), directed = TRUE)
    Overlap_network_nearest_rewired <- as_data_frame(Overlap_network_nearest_union_graph %m% Overlap_network_nearest_two_graph)
    Overlap_network_nearest_graph <- make_graph(c(t(Overlap_network_nearest[, 1:2])), directed = TRUE)
    Overlap_network_nearest_conserved_validated <- as_data_frame(Overlap_network_nearest_graph %s% miRTarget_graph)
    Overlap_network_nearest_rewired_graph <- make_graph(c(t(Overlap_network_nearest_rewired[, 1:2])), directed = TRUE)
    Overlap_network_nearest_rewired_validated <- as_data_frame(Overlap_network_nearest_rewired_graph %s% miRTarget_graph)

    Overlap_network_spline <- Overlap_net(CSmiR_network_spline, Intersect_num = round(length(CSmiR_network_spline)*0.9))
    Overlap_network_spline_union <- Overlap_net(CSmiR_network_spline, Intersect_num = 1)
    Overlap_network_spline_two <- Overlap_net(CSmiR_network_spline, Intersect_num = 2)
    Overlap_network_spline_union_graph <- make_graph(c(t(Overlap_network_spline_union[, 1:2])), directed = TRUE)
    Overlap_network_spline_two_graph <- make_graph(c(t(Overlap_network_spline_two[, 1:2])), directed = TRUE)
    Overlap_network_spline_rewired <- as_data_frame(Overlap_network_spline_union_graph %m% Overlap_network_spline_two_graph)
    Overlap_network_spline_graph <- make_graph(c(t(Overlap_network_spline[, 1:2])), directed = TRUE)
    Overlap_network_spline_conserved_validated <- as_data_frame(Overlap_network_spline_graph %s% miRTarget_graph)
    Overlap_network_spline_rewired_graph <- make_graph(c(t(Overlap_network_spline_rewired[, 1:2])), directed = TRUE)
    Overlap_network_spline_rewired_validated <- as_data_frame(Overlap_network_spline_rewired_graph %s% miRTarget_graph)

    Overlap_network_cubic <- Overlap_net(CSmiR_network_cubic, Intersect_num = round(length(CSmiR_network_cubic)*0.9))
    Overlap_network_cubic_union <- Overlap_net(CSmiR_network_cubic, Intersect_num = 1)
    Overlap_network_cubic_two <- Overlap_net(CSmiR_network_cubic, Intersect_num = 2)
    Overlap_network_cubic_union_graph <- make_graph(c(t(Overlap_network_cubic_union[, 1:2])), directed = TRUE)
    Overlap_network_cubic_two_graph <- make_graph(c(t(Overlap_network_cubic_two[, 1:2])), directed = TRUE)
    Overlap_network_cubic_rewired <- as_data_frame(Overlap_network_cubic_union_graph %m% Overlap_network_cubic_two_graph)
    Overlap_network_cubic_graph <- make_graph(c(t(Overlap_network_cubic[, 1:2])), directed = TRUE)
    Overlap_network_cubic_conserved_validated <- as_data_frame(Overlap_network_cubic_graph %s% miRTarget_graph)
    Overlap_network_cubic_rewired_graph <- make_graph(c(t(Overlap_network_cubic_rewired[, 1:2])), directed = TRUE)
    Overlap_network_cubic_rewired_validated <- as_data_frame(Overlap_network_cubic_rewired_graph %s% miRTarget_graph)

## CML-related conserved and rewired miRNA-mRNA regulatory network    
    Overlap_network_constant_CML <- Overlap_network_constant[intersect(which(Overlap_network_constant[, 1] %in% as.matrix(CML)), which(Overlap_network_constant[, 2] %in% as.matrix(CML))), ]
    Overlap_network_constant_rewired_CML <- Overlap_network_constant_rewired[intersect(which(Overlap_network_constant_rewired[, 1] %in% as.matrix(CML)), which(Overlap_network_constant_rewired[, 2] %in% as.matrix(CML))), ]

    Overlap_network_linear_CML <- Overlap_network_linear[intersect(which(Overlap_network_linear[, 1] %in% as.matrix(CML)), which(Overlap_network_linear[, 2] %in% as.matrix(CML))), ]
    Overlap_network_linear_rewired_CML <- Overlap_network_linear_rewired[intersect(which(Overlap_network_linear_rewired[, 1] %in% as.matrix(CML)), which(Overlap_network_linear_rewired[, 2] %in% as.matrix(CML))), ]

    Overlap_network_nearest_CML <- Overlap_network_nearest[intersect(which(Overlap_network_nearest[, 1] %in% as.matrix(CML)), which(Overlap_network_nearest[, 2] %in% as.matrix(CML))), ]
    Overlap_network_nearest_rewired_CML <- Overlap_network_nearest_rewired[intersect(which(Overlap_network_nearest_rewired[, 1] %in% as.matrix(CML)), which(Overlap_network_nearest_rewired[, 2] %in% as.matrix(CML))), ]

    Overlap_network_spline_CML <- Overlap_network_spline[intersect(which(Overlap_network_spline[, 1] %in% as.matrix(CML)), which(Overlap_network_spline[, 2] %in% as.matrix(CML))), ]
    Overlap_network_spline_rewired_CML <- Overlap_network_spline_rewired[intersect(which(Overlap_network_spline_rewired[, 1] %in% as.matrix(CML)), which(Overlap_network_spline_rewired[, 2] %in% as.matrix(CML))), ]

    Overlap_network_cubic_CML <- Overlap_network_cubic[intersect(which(Overlap_network_cubic[, 1] %in% as.matrix(CML)), which(Overlap_network_cubic[, 2] %in% as.matrix(CML))), ]
    Overlap_network_cubic_rewired_CML <- Overlap_network_cubic_rewired[intersect(which(Overlap_network_cubic_rewired[, 1] %in% as.matrix(CML)), which(Overlap_network_cubic_rewired[, 2] %in% as.matrix(CML))), ]

## Identifying cell-specific hub miRNAs    
    CSmiR_network_constant_outdegree <- lapply(seq(CSmiR_network_constant), function(i) degree(CSmiR_network_constant_graph[[i]], mode="out"))
    hub_miRNAs_constant <- lapply(seq(CSmiR_network_constant), function(i) names(sort(CSmiR_network_constant_outdegree[[i]][which(CSmiR_network_constant_outdegree[[i]]!=0)], decreasing=TRUE))[1:ceiling(0.2*length(which(CSmiR_network_constant_outdegree[[i]]!=0)))])
    
    CSmiR_network_linear_outdegree <- lapply(seq(CSmiR_network_linear), function(i) degree(CSmiR_network_linear_graph[[i]], mode="out"))
    hub_miRNAs_linear <- lapply(seq(CSmiR_network_linear), function(i) names(sort(CSmiR_network_linear_outdegree[[i]][which(CSmiR_network_linear_outdegree[[i]]!=0)], decreasing=TRUE))[1:ceiling(0.2*length(which(CSmiR_network_linear_outdegree[[i]]!=0)))])
    
    CSmiR_network_nearest_outdegree <- lapply(seq(CSmiR_network_nearest), function(i) degree(CSmiR_network_nearest_graph[[i]], mode="out"))
    hub_miRNAs_nearest <- lapply(seq(CSmiR_network_nearest), function(i) names(sort(CSmiR_network_nearest_outdegree[[i]][which(CSmiR_network_nearest_outdegree[[i]]!=0)], decreasing=TRUE))[1:ceiling(0.2*length(which(CSmiR_network_nearest_outdegree[[i]]!=0)))])

    CSmiR_network_spline_outdegree <- lapply(seq(CSmiR_network_spline), function(i) degree(CSmiR_network_spline_graph[[i]], mode="out"))
    hub_miRNAs_spline <- lapply(seq(CSmiR_network_spline), function(i) names(sort(CSmiR_network_spline_outdegree[[i]][which(CSmiR_network_spline_outdegree[[i]]!=0)], decreasing=TRUE))[1:ceiling(0.2*length(which(CSmiR_network_spline_outdegree[[i]]!=0)))])

    CSmiR_network_cubic_outdegree <- lapply(seq(CSmiR_network_cubic), function(i) degree(CSmiR_network_cubic_graph[[i]], mode="out"))
    hub_miRNAs_cubic <- lapply(seq(CSmiR_network_cubic), function(i) names(sort(CSmiR_network_cubic_outdegree[[i]][which(CSmiR_network_cubic_outdegree[[i]]!=0)], decreasing=TRUE))[1:ceiling(0.2*length(which(CSmiR_network_cubic_outdegree[[i]]!=0)))])

## Overlap of cell-specific hub miRNAs across cells    
    Overlap_hub_miRNAs_constant_conserved <- Overlap_hub(hub_miRNAs_constant, Intersect_num = round(length(hub_miRNAs_constant)*0.9))
    Overlap_hub_miRNAs_constant_union <- Overlap_hub(hub_miRNAs_constant, Intersect_num = 1)
    Overlap_hub_miRNAs_constant_two <- Overlap_hub(hub_miRNAs_constant, Intersect_num = 2)
    Overlap_hub_miRNAs_constant_rewired <- setdiff(Overlap_hub_miRNAs_constant_union, Overlap_hub_miRNAs_constant_two)

    Overlap_hub_miRNAs_linear_conserved <- Overlap_hub(hub_miRNAs_linear, Intersect_num = round(length(hub_miRNAs_linear)*0.9))
    Overlap_hub_miRNAs_linear_union <- Overlap_hub(hub_miRNAs_linear, Intersect_num = 1)
    Overlap_hub_miRNAs_linear_two <- Overlap_hub(hub_miRNAs_linear, Intersect_num = 2)
    Overlap_hub_miRNAs_linear_rewired <- setdiff(Overlap_hub_miRNAs_linear_union, Overlap_hub_miRNAs_linear_two)

    Overlap_hub_miRNAs_nearest_conserved <- Overlap_hub(hub_miRNAs_nearest, Intersect_num = round(length(hub_miRNAs_nearest)*0.9))
    Overlap_hub_miRNAs_nearest_union <- Overlap_hub(hub_miRNAs_nearest, Intersect_num = 1)
    Overlap_hub_miRNAs_nearest_two <- Overlap_hub(hub_miRNAs_nearest, Intersect_num = 2)
    Overlap_hub_miRNAs_nearest_rewired <- setdiff(Overlap_hub_miRNAs_nearest_union, Overlap_hub_miRNAs_nearest_two)

    Overlap_hub_miRNAs_spline_conserved <- Overlap_hub(hub_miRNAs_spline, Intersect_num = round(length(hub_miRNAs_spline)*0.9))
    Overlap_hub_miRNAs_spline_union <- Overlap_hub(hub_miRNAs_spline, Intersect_num = 1)
    Overlap_hub_miRNAs_spline_two <- Overlap_hub(hub_miRNAs_spline, Intersect_num = 2)
    Overlap_hub_miRNAs_spline_rewired <- setdiff(Overlap_hub_miRNAs_spline_union, Overlap_hub_miRNAs_spline_two)

    Overlap_hub_miRNAs_cubic_conserved <- Overlap_hub(hub_miRNAs_cubic, Intersect_num = round(length(hub_miRNAs_cubic)*0.9))
    Overlap_hub_miRNAs_cubic_union <- Overlap_hub(hub_miRNAs_cubic, Intersect_num = 1)
    Overlap_hub_miRNAs_cubic_two <- Overlap_hub(hub_miRNAs_cubic, Intersect_num = 2)
    Overlap_hub_miRNAs_cubic_rewired <- setdiff(Overlap_hub_miRNAs_cubic_union, Overlap_hub_miRNAs_cubic_two)

## CML-related cell-specific hub miRNAs    
    hub_miRNAs_constant_CML <- lapply(seq(hub_miRNAs_constant), function(i) hub_miRNAs_constant[[i]][which(hub_miRNAs_constant[[i]] %in% as.matrix(CML))])

    hub_miRNAs_linear_CML <- lapply(seq(hub_miRNAs_linear), function(i) hub_miRNAs_linear[[i]][which(hub_miRNAs_linear[[i]] %in% as.matrix(CML))])

    hub_miRNAs_nearest_CML <- lapply(seq(hub_miRNAs_nearest), function(i) hub_miRNAs_nearest[[i]][which(hub_miRNAs_nearest[[i]] %in% as.matrix(CML))])

    hub_miRNAs_spline_CML <- lapply(seq(hub_miRNAs_spline), function(i) hub_miRNAs_spline[[i]][which(hub_miRNAs_spline[[i]] %in% as.matrix(CML))])
    
    hub_miRNAs_cubic_CML <- lapply(seq(hub_miRNAs_cubic), function(i) hub_miRNAs_cubic[[i]][which(hub_miRNAs_cubic[[i]] %in% as.matrix(CML))])

## CML-related conserved and rewired hub miRNAs    
    Overlap_hub_miRNAs_constant_conserved_CML <- Overlap_hub_miRNAs_constant_conserved[which(Overlap_hub_miRNAs_constant_conserved %in% as.matrix(CML))]    
    Overlap_hub_miRNAs_constant_rewired_CML <- Overlap_hub_miRNAs_constant_rewired[which(Overlap_hub_miRNAs_constant_rewired %in% as.matrix(CML))]

    Overlap_hub_miRNAs_linear_conserved_CML <- Overlap_hub_miRNAs_linear_conserved[which(Overlap_hub_miRNAs_linear_conserved %in% as.matrix(CML))]    
    Overlap_hub_miRNAs_linear_rewired_CML <- Overlap_hub_miRNAs_linear_rewired[which(Overlap_hub_miRNAs_linear_rewired %in% as.matrix(CML))]

    Overlap_hub_miRNAs_nearest_conserved_CML <- Overlap_hub_miRNAs_nearest_conserved[which(Overlap_hub_miRNAs_nearest_conserved %in% as.matrix(CML))]    
    Overlap_hub_miRNAs_nearest_rewired_CML <- Overlap_hub_miRNAs_nearest_rewired[which(Overlap_hub_miRNAs_nearest_rewired %in% as.matrix(CML))]

    Overlap_hub_miRNAs_spline_conserved_CML <- Overlap_hub_miRNAs_spline_conserved[which(Overlap_hub_miRNAs_spline_conserved %in% as.matrix(CML))]    
    Overlap_hub_miRNAs_spline_rewired_CML <- Overlap_hub_miRNAs_spline_rewired[which(Overlap_hub_miRNAs_spline_rewired %in% as.matrix(CML))]

    Overlap_hub_miRNAs_cubic_conserved_CML <- Overlap_hub_miRNAs_cubic_conserved[which(Overlap_hub_miRNAs_cubic_conserved %in% as.matrix(CML))]    
    Overlap_hub_miRNAs_cubic_rewired_CML <- Overlap_hub_miRNAs_cubic_rewired[which(Overlap_hub_miRNAs_cubic_rewired %in% as.matrix(CML))]

## Calculating similarity matrix of cell-specific miRNA-mRNA regulatory network across cells    
    CSmiR_network_constant_Sim <- Sim.network(CSmiR_network_constant, CSmiR_network_constant, directed = TRUE)

    CSmiR_network_linear_Sim <- Sim.network(CSmiR_network_linear, CSmiR_network_linear, directed = TRUE)

    CSmiR_network_nearest_Sim <- Sim.network(CSmiR_network_nearest, CSmiR_network_nearest, directed = TRUE)

    CSmiR_network_spline_Sim <- Sim.network(CSmiR_network_spline, CSmiR_network_spline, directed = TRUE)

    CSmiR_network_cubic_Sim <- Sim.network(CSmiR_network_cubic, CSmiR_network_cubic, directed = TRUE)

## Calculating similarity matrix of cell-specific hub miRNAs across cells    
    CSmiR_hub_constant_Sim <- Sim.hub(hub_miRNAs_constant, hub_miRNAs_constant)

    CSmiR_hub_linear_Sim <- Sim.hub(hub_miRNAs_linear, hub_miRNAs_linear)

    CSmiR_hub_nearest_Sim <- Sim.hub(hub_miRNAs_nearest, hub_miRNAs_nearest)

    CSmiR_hub_spline_Sim <- Sim.hub(hub_miRNAs_spline, hub_miRNAs_spline)

    CSmiR_hub_cubic_Sim <- Sim.hub(hub_miRNAs_cubic, hub_miRNAs_cubic)

## Identifying cell-cell crosstalk network in terms of network similarity matrix    
    CSmiR_network_constant_adjacency_matrix <- ifelse(CSmiR_network_constant_Sim > median(CSmiR_network_constant_Sim[lower.tri(CSmiR_network_constant_Sim)]), 1, 0)
    diag(CSmiR_network_constant_adjacency_matrix) <- 0
    CSmiR_network_constant_adjacency_matrix_graph <- graph_from_adjacency_matrix(CSmiR_network_constant_adjacency_matrix, mode = "undirected")

    CSmiR_network_linear_adjacency_matrix <- ifelse(CSmiR_network_linear_Sim > median(CSmiR_network_linear_Sim[lower.tri(CSmiR_network_linear_Sim)]), 1, 0)
    diag(CSmiR_network_linear_adjacency_matrix) <- 0
    CSmiR_network_linear_adjacency_matrix_graph <- graph_from_adjacency_matrix(CSmiR_network_linear_adjacency_matrix, mode = "undirected")

    CSmiR_network_nearest_adjacency_matrix <- ifelse(CSmiR_network_nearest_Sim > median(CSmiR_network_nearest_Sim[lower.tri(CSmiR_network_nearest_Sim)]), 1, 0)
    diag(CSmiR_network_nearest_adjacency_matrix) <- 0
    CSmiR_network_nearest_adjacency_matrix_graph <- graph_from_adjacency_matrix(CSmiR_network_nearest_adjacency_matrix, mode = "undirected")

    CSmiR_network_spline_adjacency_matrix <- ifelse(CSmiR_network_spline_Sim > median(CSmiR_network_spline_Sim[lower.tri(CSmiR_network_spline_Sim)]), 1, 0)    
    diag(CSmiR_network_spline_adjacency_matrix) <- 0
    CSmiR_network_spline_adjacency_matrix_graph <- graph_from_adjacency_matrix(CSmiR_network_spline_adjacency_matrix, mode = "undirected")

    CSmiR_network_cubic_adjacency_matrix <- ifelse(CSmiR_network_cubic_Sim > median(CSmiR_network_cubic_Sim[lower.tri(CSmiR_network_cubic_Sim)]), 1, 0)
    diag(CSmiR_network_cubic_adjacency_matrix) <- 0
    CSmiR_network_cubic_adjacency_matrix_graph <- graph_from_adjacency_matrix(CSmiR_network_cubic_adjacency_matrix, mode = "undirected")

## Identifying cell-cell crosstalk network in terms of hub miRNA similarity matrix    
    CSmiR_hub_constant_adjacency_matrix <- ifelse(CSmiR_hub_constant_Sim > median(CSmiR_hub_constant_Sim[lower.tri(CSmiR_hub_constant_Sim)]), 1, 0)
    diag(CSmiR_hub_constant_adjacency_matrix) <- 0
    CSmiR_hub_constant_adjacency_matrix_graph <- graph_from_adjacency_matrix(CSmiR_hub_constant_adjacency_matrix, mode = "undirected")

    CSmiR_hub_linear_adjacency_matrix <- ifelse(CSmiR_hub_linear_Sim > median(CSmiR_hub_linear_Sim[lower.tri(CSmiR_hub_linear_Sim)]), 1, 0)
    diag(CSmiR_hub_linear_adjacency_matrix) <- 0
    CSmiR_hub_linear_adjacency_matrix_graph <- graph_from_adjacency_matrix(CSmiR_hub_linear_adjacency_matrix, mode = "undirected")

    CSmiR_hub_nearest_adjacency_matrix <- ifelse(CSmiR_hub_nearest_Sim > median(CSmiR_hub_nearest_Sim[lower.tri(CSmiR_hub_nearest_Sim)]), 1, 0)
    diag(CSmiR_hub_nearest_adjacency_matrix) <- 0
    CSmiR_hub_nearest_adjacency_matrix_graph <- graph_from_adjacency_matrix(CSmiR_hub_nearest_adjacency_matrix, mode = "undirected")

    CSmiR_hub_spline_adjacency_matrix <- ifelse(CSmiR_hub_spline_Sim > median(CSmiR_hub_spline_Sim[lower.tri(CSmiR_hub_spline_Sim)]), 1, 0)    
    diag(CSmiR_hub_spline_adjacency_matrix) <- 0
    CSmiR_hub_spline_adjacency_matrix_graph <- graph_from_adjacency_matrix(CSmiR_hub_spline_adjacency_matrix, mode = "undirected")

    CSmiR_hub_cubic_adjacency_matrix <- ifelse(CSmiR_hub_cubic_Sim > median(CSmiR_hub_cubic_Sim[lower.tri(CSmiR_hub_cubic_Sim)]), 1, 0)
    diag(CSmiR_hub_cubic_adjacency_matrix) <- 0
    CSmiR_hub_cubic_adjacency_matrix_graph <- graph_from_adjacency_matrix(CSmiR_hub_cubic_adjacency_matrix, mode = "undirected")

## Identifying hub cells in terms of network similarity    
    CSmiR_network_constant_cell_degree <- degree(CSmiR_network_constant_adjacency_matrix_graph)
    CSmiR_network_constant_hub_cells <- names(sort(CSmiR_network_constant_cell_degree[which(CSmiR_network_constant_cell_degree!=0)], decreasing=TRUE))[1:ceiling(0.2*length(which(CSmiR_network_constant_cell_degree!=0)))]

    CSmiR_network_linear_cell_degree <- degree(CSmiR_network_linear_adjacency_matrix_graph)
    CSmiR_network_linear_hub_cells <- names(sort(CSmiR_network_linear_cell_degree[which(CSmiR_network_linear_cell_degree!=0)], decreasing=TRUE))[1:ceiling(0.2*length(which(CSmiR_network_linear_cell_degree!=0)))]

    CSmiR_network_nearest_cell_degree <- degree(CSmiR_network_nearest_adjacency_matrix_graph)
    CSmiR_network_nearest_hub_cells <- names(sort(CSmiR_network_nearest_cell_degree[which(CSmiR_network_nearest_cell_degree!=0)], decreasing=TRUE))[1:ceiling(0.2*length(which(CSmiR_network_nearest_cell_degree!=0)))]

    CSmiR_network_spline_cell_degree <- degree(CSmiR_network_spline_adjacency_matrix_graph)
    CSmiR_network_spline_hub_cells <- names(sort(CSmiR_network_spline_cell_degree[which(CSmiR_network_spline_cell_degree!=0)], decreasing=TRUE))[1:ceiling(0.2*length(which(CSmiR_network_spline_cell_degree!=0)))]

    CSmiR_network_cubic_cell_degree <- degree(CSmiR_network_cubic_adjacency_matrix_graph)
    CSmiR_network_cubic_hub_cells <- names(sort(CSmiR_network_cubic_cell_degree[which(CSmiR_network_cubic_cell_degree!=0)], decreasing=TRUE))[1:ceiling(0.2*length(which(CSmiR_network_cubic_cell_degree!=0)))]

## Identifying hub cells in terms of hub miRNA similarity    
    CSmiR_hub_constant_cell_degree <- degree(CSmiR_hub_constant_adjacency_matrix_graph)
    CSmiR_hub_constant_hub_cells <- names(sort(CSmiR_hub_constant_cell_degree[which(CSmiR_hub_constant_cell_degree!=0)], decreasing=TRUE))[1:ceiling(0.2*length(which(CSmiR_hub_constant_cell_degree!=0)))]

    CSmiR_hub_linear_cell_degree <- degree(CSmiR_hub_linear_adjacency_matrix_graph)
    CSmiR_hub_linear_hub_cells <- names(sort(CSmiR_hub_linear_cell_degree[which(CSmiR_hub_linear_cell_degree!=0)], decreasing=TRUE))[1:ceiling(0.2*length(which(CSmiR_hub_linear_cell_degree!=0)))]

    CSmiR_hub_nearest_cell_degree <- degree(CSmiR_hub_nearest_adjacency_matrix_graph)
    CSmiR_hub_nearest_hub_cells <- names(sort(CSmiR_hub_nearest_cell_degree[which(CSmiR_hub_nearest_cell_degree!=0)], decreasing=TRUE))[1:ceiling(0.2*length(which(CSmiR_hub_nearest_cell_degree!=0)))]

    CSmiR_hub_spline_cell_degree <- degree(CSmiR_hub_spline_adjacency_matrix_graph)
    CSmiR_hub_spline_hub_cells <- names(sort(CSmiR_hub_spline_cell_degree[which(CSmiR_hub_spline_cell_degree!=0)], decreasing=TRUE))[1:ceiling(0.2*length(which(CSmiR_hub_spline_cell_degree!=0)))]

    CSmiR_hub_cubic_cell_degree <- degree(CSmiR_hub_cubic_adjacency_matrix_graph)
    CSmiR_hub_cubic_hub_cells <- names(sort(CSmiR_hub_cubic_cell_degree[which(CSmiR_hub_cubic_cell_degree!=0)], decreasing=TRUE))[1:ceiling(0.2*length(which(CSmiR_hub_cubic_cell_degree!=0)))]

## Identifying cell-cell crosstalk modules in terms of network similarity matrix   
    CSmiR_network_constant_cell_module <- netModule(CSmiR_network_constant_adjacency_matrix_graph %>% as_data_frame)

    CSmiR_network_linear_cell_module <- netModule(CSmiR_network_linear_adjacency_matrix_graph %>% as_data_frame)

    CSmiR_network_nearest_cell_module <- netModule(CSmiR_network_nearest_adjacency_matrix_graph %>% as_data_frame)

    CSmiR_network_spline_cell_module <- netModule(CSmiR_network_spline_adjacency_matrix_graph %>% as_data_frame)

    CSmiR_network_cubic_cell_module <- netModule(CSmiR_network_cubic_adjacency_matrix_graph %>% as_data_frame)

## Identifying cell-cell crosstalk modules in terms of hub miRNA similarity matrix    
    CSmiR_hub_constant_cell_module <- netModule(CSmiR_hub_constant_adjacency_matrix_graph %>% as_data_frame)

    CSmiR_hub_linear_cell_module <- netModule(CSmiR_hub_linear_adjacency_matrix_graph %>% as_data_frame)

    CSmiR_hub_nearest_cell_module <- netModule(CSmiR_hub_nearest_adjacency_matrix_graph %>% as_data_frame)

    CSmiR_hub_spline_cell_module <- netModule(CSmiR_hub_spline_adjacency_matrix_graph %>% as_data_frame)

    CSmiR_hub_cubic_cell_module <- netModule(CSmiR_hub_cubic_adjacency_matrix_graph %>% as_data_frame)

## The miR-17/92 family analysis
    miR_17_92_family <- c("hsa-miR-17-3p", "hsa-miR-17-5p","hsa-miR-18a-3p", "hsa-miR-18a-5p", "hsa-miR-19a-3p", "hsa-miR-19a-5p",  
                          "hsa-miR-19b-3p", "hsa-miR-19b-1-5p", "hsa-miR-20a-3p", "hsa-miR-20a-5p", "hsa-miR-92a-3p", "hsa-miR-92a-1-5p")
    
    # Heatmap of miR-17/92 family expression
    miRfamily_exp <- miRNA_scRNA_norm_average[, which(colnames(miRNA_scRNA_norm_average) %in% miR_17_92_family)]
    miRfamily_exp <- t(miRfamily_exp)
    colnames(miRfamily_exp) <- c(1:19)    
    pheatmap(miRfamily_exp, cluster_rows = FALSE, cluster_cols = FALSE)
    
    # Extracting miR-17/92 family related cell-specific miRNA-mRNA regulatory networks
    CSmiR_network_constant_miRfamily <- lapply(seq(CSmiR_network_constant), function(i) CSmiR_network_constant[[i]][which(CSmiR_network_constant[[i]][, 1] %in% as.matrix(miR_17_92_family)), ])
    
    CSmiR_network_linear_miRfamily <- lapply(seq(CSmiR_network_linear), function(i) CSmiR_network_linear[[i]][which(CSmiR_network_linear[[i]][, 1] %in% as.matrix(miR_17_92_family)), ])
    
    CSmiR_network_nearest_miRfamily <- lapply(seq(CSmiR_network_nearest), function(i) CSmiR_network_nearest[[i]][which(CSmiR_network_nearest[[i]][, 1] %in% as.matrix(miR_17_92_family)), ])
 
    CSmiR_network_spline_miRfamily <- lapply(seq(CSmiR_network_spline), function(i) CSmiR_network_spline[[i]][which(CSmiR_network_spline[[i]][, 1] %in% as.matrix(miR_17_92_family)), ])
    
    CSmiR_network_cubic_miRfamily <- lapply(seq(CSmiR_network_cubic), function(i) CSmiR_network_cubic[[i]][which(CSmiR_network_cubic[[i]][, 1] %in% as.matrix(miR_17_92_family)), ])
 
    # Extracting conserved and rewired miRNA-mRNA regulatory networks associated with the miR-17/92 family 
    Overlap_network_constant_miRfamily <- Overlap_network_constant[which(Overlap_network_constant[, 1] %in% as.matrix(miR_17_92_family)), ]
    Overlap_network_constant_rewired_miRfamily <- Overlap_network_constant_rewired[which(Overlap_network_constant_rewired[, 1] %in% as.matrix(miR_17_92_family)), ]

    Overlap_network_linear_miRfamily <- Overlap_network_linear[which(Overlap_network_linear[, 1] %in% as.matrix(miR_17_92_family)), ]
    Overlap_network_linear_rewired_miRfamily <- Overlap_network_linear_rewired[which(Overlap_network_linear_rewired[, 1] %in% as.matrix(miR_17_92_family)), ]

    Overlap_network_nearest_miRfamily <- Overlap_network_nearest[which(Overlap_network_nearest[, 1] %in% as.matrix(miR_17_92_family)), ]
    Overlap_network_nearest_rewired_miRfamily <- Overlap_network_nearest_rewired[which(Overlap_network_nearest_rewired[, 1] %in% as.matrix(miR_17_92_family)), ]

    Overlap_network_spline_miRfamily <- Overlap_network_spline[which(Overlap_network_spline[, 1] %in% as.matrix(miR_17_92_family)), ]
    Overlap_network_spline_rewired_miRfamily <- Overlap_network_spline_rewired[which(Overlap_network_spline_rewired[, 1] %in% as.matrix(miR_17_92_family)), ]

    Overlap_network_cubic_miRfamily <- Overlap_network_cubic[which(Overlap_network_cubic[, 1] %in% as.matrix(miR_17_92_family)), ]
    Overlap_network_cubic_rewired_miRfamily <- Overlap_network_cubic_rewired[which(Overlap_network_cubic_rewired[, 1] %in% as.matrix(miR_17_92_family)), ]
    
    # Extracting CML-related conserved and rewired miRNA-mRNA regulatory networks associated with the miR-17/92 family
    CSmiR_network_constant_miRfamily_CML <- lapply(seq(CSmiR_network_constant_miRfamily), function(i) CSmiR_network_constant_miRfamily[[i]][intersect(which(CSmiR_network_constant_miRfamily[[i]][, 1] %in% as.matrix(CML)), which(CSmiR_network_constant_miRfamily[[i]][, 2] %in% as.matrix(CML))), ])

    CSmiR_network_linear_miRfamily_CML <- lapply(seq(CSmiR_network_linear_miRfamily), function(i) CSmiR_network_linear_miRfamily[[i]][intersect(which(CSmiR_network_linear_miRfamily[[i]][, 1] %in% as.matrix(CML)), which(CSmiR_network_linear_miRfamily[[i]][, 2] %in% as.matrix(CML))), ])
    
    CSmiR_network_nearest_miRfamily_CML <- lapply(seq(CSmiR_network_nearest_miRfamily), function(i) CSmiR_network_nearest_miRfamily[[i]][intersect(which(CSmiR_network_nearest_miRfamily[[i]][, 1] %in% as.matrix(CML)), which(CSmiR_network_nearest_miRfamily[[i]][, 2] %in% as.matrix(CML))), ])

    CSmiR_network_spline_miRfamily_CML <- lapply(seq(CSmiR_network_spline_miRfamily), function(i) CSmiR_network_spline_miRfamily[[i]][intersect(which(CSmiR_network_spline_miRfamily[[i]][, 1] %in% as.matrix(CML)), which(CSmiR_network_spline_miRfamily[[i]][, 2] %in% as.matrix(CML))), ])
    
    CSmiR_network_cubic_miRfamily_CML <- lapply(seq(CSmiR_network_cubic_miRfamily), function(i) CSmiR_network_cubic_miRfamily[[i]][intersect(which(CSmiR_network_cubic_miRfamily[[i]][, 1] %in% as.matrix(CML)), which(CSmiR_network_cubic_miRfamily[[i]][, 2] %in% as.matrix(CML))), ])
    
    # Experimentally validated conserved and rewired miRNA-mRNA interactions associated with the miR-17/92 family
    CSmiR_network_constant_miRfamily_graph <- lapply(seq(CSmiR_network_constant_miRfamily), function(i) make_graph(c(t(CSmiR_network_constant_miRfamily[[i]][, 1:2])), directed = TRUE))
    CSmiR_network_constant_miRfamily_validated <- lapply(seq(CSmiR_network_constant_miRfamily), function(i) as_data_frame(CSmiR_network_constant_miRfamily_graph[[i]] %s% miRTarget_graph))

    CSmiR_network_linear_miRfamily_graph <- lapply(seq(CSmiR_network_linear_miRfamily), function(i) make_graph(c(t(CSmiR_network_linear_miRfamily[[i]][, 1:2])), directed = TRUE))
    CSmiR_network_linear_miRfamily_validated <- lapply(seq(CSmiR_network_linear_miRfamily), function(i) as_data_frame(CSmiR_network_linear_miRfamily_graph[[i]] %s% miRTarget_graph))

    CSmiR_network_nearest_miRfamily_graph <- lapply(seq(CSmiR_network_nearest_miRfamily), function(i) make_graph(c(t(CSmiR_network_nearest_miRfamily[[i]][, 1:2])), directed = TRUE))
    CSmiR_network_nearest_miRfamily_validated <- lapply(seq(CSmiR_network_nearest_miRfamily), function(i) as_data_frame(CSmiR_network_nearest_miRfamily_graph[[i]] %s% miRTarget_graph))

    CSmiR_network_spline_miRfamily_graph <- lapply(seq(CSmiR_network_spline_miRfamily), function(i) make_graph(c(t(CSmiR_network_spline_miRfamily[[i]][, 1:2])), directed = TRUE))
    CSmiR_network_spline_miRfamily_validated <- lapply(seq(CSmiR_network_spline_miRfamily), function(i) as_data_frame(CSmiR_network_spline_miRfamily_graph[[i]] %s% miRTarget_graph))

    CSmiR_network_cubic_miRfamily_graph <- lapply(seq(CSmiR_network_cubic_miRfamily), function(i) make_graph(c(t(CSmiR_network_cubic_miRfamily[[i]][, 1:2])), directed = TRUE))
    CSmiR_network_cubic_miRfamily_validated <- lapply(seq(CSmiR_network_cubic_miRfamily), function(i) as_data_frame(CSmiR_network_cubic_miRfamily_graph[[i]] %s% miRTarget_graph))

## Discovering maximal bicliques of conserved and rewired the miR-17/92 family regulation   
    CSmiR_network_miRfamily_constant_biclique <- biclique_network(list(Sub_miR(Overlap_network_constant_miRfamily), Sub_miR(Overlap_network_constant_rewired_miRfamily)), lleast = 2, rleast = 3)
   
    CSmiR_network_miRfamily_linear_biclique <- biclique_network(list(Sub_miR(Overlap_network_linear_miRfamily), Sub_miR(Overlap_network_linear_rewired_miRfamily)), lleast = 2, rleast = 3)
    
    CSmiR_network_miRfamily_nearest_biclique <- biclique_network(list(Sub_miR(Overlap_network_nearest_miRfamily), Sub_miR(Overlap_network_nearest_rewired_miRfamily)), lleast = 2, rleast = 3)
    
    CSmiR_network_miRfamily_spline_biclique <- biclique_network(list(Sub_miR(Overlap_network_spline_miRfamily), Sub_miR(Overlap_network_spline_rewired_miRfamily)), lleast = 2, rleast = 3)
    
    #No rewired bicliques found
    CSmiR_network_miRfamily_cubic_biclique <- biclique_network(list(Sub_miR(Overlap_network_cubic_miRfamily)), lleast = 2, rleast = 3)

## Enrichment analysis of maximal bicliques of conserved and rewired the miR-17/92 family regulation   
   miRfamily_constant_biclique_conserved_gene_list <- lapply(seq(CSmiR_network_miRfamily_constant_biclique[[1]]), function(i) CSmiR_network_miRfamily_constant_biclique[[1]][[i]]$right)
   miRfamily_constant_biclique_rewired_gene_list <- lapply(seq(CSmiR_network_miRfamily_constant_biclique[[2]]), function(i) CSmiR_network_miRfamily_constant_biclique[[2]][[i]]$right)

   miRfamily_linear_biclique_conserved_gene_list <- lapply(seq(CSmiR_network_miRfamily_linear_biclique[[1]]), function(i) CSmiR_network_miRfamily_linear_biclique[[1]][[i]]$right)
   miRfamily_linear_biclique_rewired_gene_list <- lapply(seq(CSmiR_network_miRfamily_linear_biclique[[2]]), function(i) CSmiR_network_miRfamily_linear_biclique[[2]][[i]]$right)

   miRfamily_nearest_biclique_conserved_gene_list <- lapply(seq(CSmiR_network_miRfamily_nearest_biclique[[1]]), function(i) CSmiR_network_miRfamily_nearest_biclique[[1]][[i]]$right)
   miRfamily_nearest_biclique_rewired_gene_list <- lapply(seq(CSmiR_network_miRfamily_nearest_biclique[[2]]), function(i) CSmiR_network_miRfamily_nearest_biclique[[2]][[i]]$right)

   miRfamily_spline_biclique_conserved_gene_list <- lapply(seq(CSmiR_network_miRfamily_spline_biclique[[1]]), function(i) CSmiR_network_miRfamily_spline_biclique[[1]][[i]]$right)
   miRfamily_spline_biclique_rewired_gene_list <- lapply(seq(CSmiR_network_miRfamily_spline_biclique[[2]]), function(i) CSmiR_network_miRfamily_spline_biclique[[2]][[i]]$right)

   miRfamily_cubic_biclique_conserved_gene_list <- lapply(seq(CSmiR_network_miRfamily_cubic_biclique[[1]]), function(i) CSmiR_network_miRfamily_cubic_biclique[[1]][[i]]$right)
   
   # GO, KEGG and Reactome enrichment analysis
   miRfamily_constant_biclique_conserved_FEA <- moduleFEA(miRfamily_constant_biclique_conserved_gene_list, padjustvaluecutoff = 0.05)
   miRfamily_constant_biclique_rewired_FEA <- moduleFEA(miRfamily_constant_biclique_rewired_gene_list, padjustvaluecutoff = 0.05)

   miRfamily_linear_biclique_conserved_FEA <- moduleFEA(miRfamily_linear_biclique_conserved_gene_list, padjustvaluecutoff = 0.05)
   miRfamily_linear_biclique_rewired_FEA <- moduleFEA(miRfamily_linear_biclique_rewired_gene_list, padjustvaluecutoff = 0.05)

   miRfamily_nearest_biclique_conserved_FEA <- moduleFEA(miRfamily_nearest_biclique_conserved_gene_list, padjustvaluecutoff = 0.05)
   miRfamily_nearest_biclique_rewired_FEA <- moduleFEA(miRfamily_nearest_biclique_rewired_gene_list, padjustvaluecutoff = 0.05)

   miRfamily_spline_biclique_conserved_FEA <- moduleFEA(miRfamily_spline_biclique_conserved_gene_list, padjustvaluecutoff = 0.05)
   miRfamily_spline_biclique_rewired_FEA <- moduleFEA(miRfamily_spline_biclique_rewired_gene_list, padjustvaluecutoff = 0.05) 
   
   miRfamily_cubic_biclique_conserved_FEA <- moduleFEA(miRfamily_cubic_biclique_conserved_gene_list, padjustvaluecutoff = 0.05)
      
   miRfamily_constant_biclique_conserved_miRmR_list <- lapply(seq(CSmiR_network_miRfamily_constant_biclique[[1]]), function(i) c(CSmiR_network_miRfamily_constant_biclique[[1]][[i]]$left, CSmiR_network_miRfamily_constant_biclique[[1]][[i]]$right))
   miRfamily_constant_biclique_rewired_miRmR_list <- lapply(seq(CSmiR_network_miRfamily_constant_biclique[[2]]), function(i) c(CSmiR_network_miRfamily_constant_biclique[[2]][[i]]$left, CSmiR_network_miRfamily_constant_biclique[[2]][[i]]$right))

   miRfamily_linear_biclique_conserved_miRmR_list <- lapply(seq(CSmiR_network_miRfamily_linear_biclique[[1]]), function(i) c(CSmiR_network_miRfamily_linear_biclique[[1]][[i]]$left, CSmiR_network_miRfamily_linear_biclique[[1]][[i]]$right))
   miRfamily_linear_biclique_rewired_miRmR_list <- lapply(seq(CSmiR_network_miRfamily_linear_biclique[[2]]), function(i) c(CSmiR_network_miRfamily_linear_biclique[[2]][[i]]$left, CSmiR_network_miRfamily_linear_biclique[[2]][[i]]$right))

   miRfamily_nearest_biclique_conserved_miRmR_list <- lapply(seq(CSmiR_network_miRfamily_nearest_biclique[[1]]), function(i) c(CSmiR_network_miRfamily_nearest_biclique[[1]][[i]]$left, CSmiR_network_miRfamily_nearest_biclique[[1]][[i]]$right))
   miRfamily_nearest_biclique_rewired_miRmR_list <- lapply(seq(CSmiR_network_miRfamily_nearest_biclique[[2]]), function(i) c(CSmiR_network_miRfamily_nearest_biclique[[2]][[i]]$left, CSmiR_network_miRfamily_nearest_biclique[[2]][[i]]$right))
 
   miRfamily_spline_biclique_conserved_miRmR_list <- lapply(seq(CSmiR_network_miRfamily_spline_biclique[[1]]), function(i) c(CSmiR_network_miRfamily_spline_biclique[[1]][[i]]$left, CSmiR_network_miRfamily_spline_biclique[[1]][[i]]$right))
   miRfamily_spline_biclique_rewired_miRmR_list <- lapply(seq(CSmiR_network_miRfamily_spline_biclique[[2]]), function(i) c(CSmiR_network_miRfamily_spline_biclique[[2]][[i]]$left, CSmiR_network_miRfamily_spline_biclique[[2]][[i]]$right))

   miRfamily_cubic_biclique_conserved_miRmR_list <- lapply(seq(CSmiR_network_miRfamily_cubic_biclique[[1]]), function(i) c(CSmiR_network_miRfamily_cubic_biclique[[1]][[i]]$left, CSmiR_network_miRfamily_cubic_biclique[[1]][[i]]$right))
  
   # CML enrichment analysis   
   miRfamily_constant_biclique_conserved_CML_EA <- module_CML_EA(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, CML, miRfamily_constant_biclique_conserved_miRmR_list)
   miRfamily_constant_biclique_rewired_CML_EA <- module_CML_EA(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, CML, miRfamily_constant_biclique_rewired_miRmR_list)

   miRfamily_linear_biclique_conserved_CML_EA <- module_CML_EA(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, CML, miRfamily_linear_biclique_conserved_miRmR_list)
   miRfamily_linear_biclique_rewired_CML_EA <- module_CML_EA(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, CML, miRfamily_linear_biclique_rewired_miRmR_list)

   miRfamily_nearest_biclique_conserved_CML_EA <- module_CML_EA(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, CML, miRfamily_nearest_biclique_conserved_miRmR_list)
   miRfamily_nearest_biclique_rewired_CML_EA <- module_CML_EA(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, CML, miRfamily_nearest_biclique_rewired_miRmR_list)
   
   miRfamily_spline_biclique_conserved_CML_EA <- module_CML_EA(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, CML, miRfamily_spline_biclique_conserved_miRmR_list)
   miRfamily_spline_biclique_rewired_CML_EA <- module_CML_EA(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, CML, miRfamily_spline_biclique_rewired_miRmR_list)

   miRfamily_cubic_biclique_conserved_CML_EA <- module_CML_EA(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, CML, miRfamily_cubic_biclique_conserved_miRmR_list)
   
   # Hallmark enrichment analysis
   m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, human_gene_symbol)
   
   miRfamily_constant_biclique_conserved_Hallmark <- lapply(seq(miRfamily_constant_biclique_conserved_gene_list), function(i) enricher(miRfamily_constant_biclique_conserved_gene_list[[i]], TERM2GENE=m_t2g, minGSSize=1) %>% as.data.frame)
   miRfamily_constant_biclique_rewired_Hallmark <- lapply(seq(miRfamily_constant_biclique_rewired_gene_list), function(i) enricher(miRfamily_constant_biclique_rewired_gene_list[[i]], TERM2GENE=m_t2g, minGSSize=1) %>% as.data.frame)

   miRfamily_linear_biclique_conserved_Hallmark <- lapply(seq(miRfamily_linear_biclique_conserved_gene_list), function(i) enricher(miRfamily_linear_biclique_conserved_gene_list[[i]], TERM2GENE=m_t2g, minGSSize=1) %>% as.data.frame)
   miRfamily_linear_biclique_rewired_Hallmark <- lapply(seq(miRfamily_linear_biclique_rewired_gene_list), function(i) enricher(miRfamily_linear_biclique_rewired_gene_list[[i]], TERM2GENE=m_t2g, minGSSize=1) %>% as.data.frame)
   
   miRfamily_nearest_biclique_conserved_Hallmark <- lapply(seq(miRfamily_nearest_biclique_conserved_gene_list), function(i) enricher(miRfamily_nearest_biclique_conserved_gene_list[[i]], TERM2GENE=m_t2g, minGSSize=1) %>% as.data.frame)
   miRfamily_nearest_biclique_rewired_Hallmark <- llapply(seq(miRfamily_nearest_biclique_rewired_gene_list), function(i) enricher(miRfamily_nearest_biclique_rewired_gene_list[[i]], TERM2GENE=m_t2g, minGSSize=1) %>% as.data.frame)
   
   miRfamily_spline_biclique_conserved_Hallmark <- lapply(seq(miRfamily_spline_biclique_conserved_gene_list), function(i) enricher(miRfamil_spliney_biclique_conserved_gene_list[[i]], TERM2GENE=m_t2g, minGSSize=1) %>% as.data.frame)
   miRfamily_spline_biclique_rewired_Hallmark <- lapply(seq(miRfamily_spline_biclique_rewired_gene_list), function(i) enricher(miRfamily_spline_biclique_rewired_gene_list[[i]], TERM2GENE=m_t2g, minGSSize=1) %>% as.data.frame)
  
   miRfamily_cubic_biclique_conserved_Hallmark <- lapply(seq(miRfamily_cubic_biclique_conserved_gene_list), function(i) enricher(miRfamily_cubic_biclique_conserved_gene_list[[i]], TERM2GENE=m_t2g, minGSSize=1) %>% as.data.frame)
  
   # Cell marker enrichment analsyis
   cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>%
   tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
   dplyr::select(cellMarker, geneSymbol)    
   
   miRfamily_constant_biclique_conserved_Cellmarker <- lapply(seq(miRfamily_constant_biclique_conserved_gene_list), function(i) enricher(miRfamily_constant_biclique_conserved_gene_list[[i]], TERM2GENE=cell_markers, minGSSize=1) %>% as.data.frame)
   miRfamily_constant_biclique_rewired_Cellmarker <- lapply(seq(miRfamily_constant_biclique_rewired_gene_list), function(i) enricher(miRfamily_constant_biclique_rewired_gene_list[[i]], TERM2GENE=cell_markers, minGSSize=1) %>% as.data.frame)

   miRfamily_linear_biclique_conserved_Cellmarker <- lapply(seq(miRfamily_linear_biclique_conserved_gene_list), function(i) enricher(miRfamily_linear_biclique_conserved_gene_list[[i]], TERM2GENE=cell_markers, minGSSize=1) %>% as.data.frame)
   miRfamily_linear_biclique_rewired_Cellmarker <- lapply(seq(miRfamily_linear_biclique_rewired_gene_list), function(i) enricher(miRfamily_linear_biclique_rewired_gene_list[[i]], TERM2GENE=cell_markers, minGSSize=1) %>% as.data.frame)
   
   miRfamily_nearest_biclique_conserved_Cellmarker <- lapply(seq(miRfamily_nearest_biclique_conserved_gene_list), function(i) enricher(miRfamily_nearest_biclique_conserved_gene_list[[i]], TERM2GENE=cell_markers, minGSSize=1) %>% as.data.frame)
   miRfamily_nearest_biclique_rewired_Cellmarker <- lapply(seq(miRfamily_nearest_biclique_rewired_gene_list), function(i) enricher(miRfamily_nearest_biclique_rewired_gene_list[[i]], TERM2GENE=cell_markers, minGSSize=1) %>% as.data.frame)

   miRfamily_spline_biclique_conserved_Cellmarker <- lapply(seq(miRfamily_spline_biclique_conserved_gene_list), function(i) enricher(miRfamily_spline_biclique_conserved_gene_list[[i]], TERM2GENE=cell_markers, minGSSize=1) %>% as.data.frame)
   miRfamily_spline_biclique_rewired_Cellmarker <- lapply(seq(miRfamily_spline_biclique_rewired_gene_list), function(i) enricher(miRfamily_spline_biclique_rewired_gene_list[[i]], TERM2GENE=cell_markers, minGSSize=1) %>% as.data.frame)

   miRfamily_cubic_biclique_conserved_Cellmarker <- lapply(seq(miRfamily_cubic_biclique_conserved_gene_list), function(i) enricher(miRfamily_cubic_biclique_conserved_gene_list[[i]], TERM2GENE=cell_markers, minGSSize=1) %>% as.data.frame)
  
## Hierarchical cluster analysis of cell-specific miRNA-mRNA regulatory network    
    rownames(CsmiR_network_constant_Sim) <- colnames(CsmiR_network_constant_Sim) <- rownames(miRNA_scRNA_norm_filter)
    rownames(CsmiR_network_linear_Sim) <- colnames(CsmiR_network_linear_Sim) <- rownames(miRNA_scRNA_norm_filter)
    rownames(CsmiR_network_nearest_Sim) <- colnames(CsmiR_network_nearest_Sim) <- rownames(miRNA_scRNA_norm_filter)
    rownames(CsmiR_network_spline_Sim) <- colnames(CsmiR_network_spline_Sim) <- rownames(miRNA_scRNA_norm_filter)
    rownames(CsmiR_network_cubic_Sim) <- colnames(CsmiR_network_cubic_Sim) <- rownames(miRNA_scRNA_norm_filter)
     
    hclust_res_network_constant <- hclust(as.dist(1-CsmiR_network_constant_Sim), "complete")
    hclust_res_network_linear <- hclust(as.dist(1-CsmiR_network_linear_Sim), "complete")
    hclust_res_network_nearest <- hclust(as.dist(1-CsmiR_network_nearest_Sim), "complete")
    hclust_res_network_spline <- hclust(as.dist(1-CsmiR_network_spline_Sim), "complete")  
    hclust_res_network_cubic <- hclust(as.dist(1-CsmiR_network_cubic_Sim), "complete")    
    
    dend_network_constant <- as.dendrogram(hclust_res_network_constant)
    dend_network_constant %>% set("branches_k_color", value = c("red", "blue"), k=2) %>% 
    set("labels_col", value = c("red", "blue"), k=2) %>% 
    set("labels_cex", value = 0.7) %>% plot(main="k=2 using network dissimilarity")
    
    dend_network_constant %>% set("branches_k_color", value = c("red", "green", "blue"), k=3) %>% 
    set("labels_col", value = c("red", "green", "blue"), k=3) %>% 
    set("labels_cex", value = 0.7) %>% plot(main="k=3 using network dissimilarity")
    
    dend_network_linear <- as.dendrogram(hclust_res_network_linear)
    dend_network_linear %>% set("branches_k_color", value = c("red", "blue"), k=2) %>% 
    set("labels_col", value = c("red", "blue"), k=2) %>% 
    set("labels_cex", value = 0.7) %>% plot(main="k=2 using network dissimilarity")
    
    dend_network_linear %>% set("branches_k_color", value = c("red", "green", "blue"), k=3) %>% 
    set("labels_col", value = c("red", "green", "blue"), k=3) %>% 
    set("labels_cex", value = 0.7) %>% plot(main="k=3 using network dissimilarity")
    
    dend_network_nearest <- as.dendrogram(hclust_res_network_nearest)
    dend_network_nearest %>% set("branches_k_color", value = c("red", "blue"), k=2) %>% 
    set("labels_col", value = c("red", "blue"), k=2) %>% 
    set("labels_cex", value = 0.7) %>% plot(main="k=2 using network dissimilarity")
    
    dend_network_nearest %>% set("branches_k_color", value = c("red", "green", "blue"), k=3) %>% 
    set("labels_col", value = c("red", "green", "blue"), k=3) %>% 
    set("labels_cex", value = 0.7) %>% plot(main="k=3 using network dissimilarity")

    dend_network_spline <- as.dendrogram(hclust_res_network_spline)
    dend_network_spline %>% set("branches_k_color", value = c("red", "blue"), k=2) %>% 
    set("labels_col", value = c("red", "blue"), k=2) %>% 
    set("labels_cex", value = 0.7) %>% plot(main="k=2 using network dissimilarity")
    
    dend_network_spline %>% set("branches_k_color", value = c("red", "green", "blue"), k=3) %>% 
    set("labels_col", value = c("red", "green", "blue"), k=3) %>% 
    set("labels_cex", value = 0.7) %>% plot(main="k=3 using network dissimilarity")
    
    dend_network_cubic <- as.dendrogram(hclust_res_network_cubic)
    dend_network_cubic %>% set("branches_k_color", value = c("red", "blue"), k=2) %>% 
    set("labels_col", value = c("red", "blue"), k=2) %>% 
    set("labels_cex", value = 0.7) %>% plot(main="k=2 using network dissimilarity")
    
    dend_network_cubic %>% set("branches_k_color", value = c("red", "green", "blue"), k=3) %>% 
    set("labels_col", value = c("red", "green", "blue"), k=3) %>% 
    set("labels_cex", value = 0.7) %>% plot(main="k=3 using network dissimilarity")

## Hierarchical cluster analysis of cell-specific hub miRNAs    
    rownames(CsmiR_hub_constant_Sim) <- colnames(CsmiR_hub_constant_Sim) <- rownames(miRNA_scRNA_norm_filter)
    rownames(CsmiR_hub_linear_Sim) <- colnames(CsmiR_hub_linear_Sim) <- rownames(miRNA_scRNA_norm_filter)
    rownames(CsmiR_hub_nearest_Sim) <- colnames(CsmiR_hub_nearest_Sim) <- rownames(miRNA_scRNA_norm_filter)
    rownames(CsmiR_hub_spline_Sim) <- colnames(CsmiR_hub_spline_Sim) <- rownames(miRNA_scRNA_norm_filter)
    rownames(CsmiR_hub_cubic_Sim) <- colnames(CsmiR_hub_cubic_Sim) <- rownames(miRNA_scRNA_norm_filter)
      
    hclust_res_hub_constant <- hclust(as.dist(1-CsmiR_hub_constant_Sim), "complete")
    hclust_res_hub_linear <- hclust(as.dist(1-CsmiR_hub_linear_Sim), "complete")
    hclust_res_hub_nearest <- hclust(as.dist(1-CsmiR_hub_nearest_Sim), "complete")
    hclust_res_hub_spline <- hclust(as.dist(1-CsmiR_hub_spline_Sim), "complete") 
    hclust_res_hub_cubic <- hclust(as.dist(1-CsmiR_hub_cubic_Sim), "complete")    

    dend_hub_constant <- as.dendrogram(hclust_res_hub_constant)
    dend_hub_constant %>% set("branches_k_color", value = c("red", "blue"), k=2) %>% 
    set("labels_col", value = c("red", "blue"), k=2) %>% 
    set("labels_cex", value = 0.7) %>% plot(main="k=2 using hub miRNA dissimilarity")
    
    dend_hub_constant %>% set("branches_k_color", value = c("red", "green", "blue"), k=3) %>% 
    set("labels_col", value = c("red", "green", "blue"), k=3) %>% 
    set("labels_cex", value = 0.7) %>% plot(main="k=3 using hub miRNA dissimilarity")
    
    dend_hub_linear <- as.dendrogram(hclust_res_hub_linear)
    dend_hub_linear %>% set("branches_k_color", value = c("red", "blue"), k=2) %>% 
    set("labels_col", value = c("red", "blue"), k=2) %>% 
    set("labels_cex", value = 0.7) %>% plot(main="k=2 using hub miRNA dissimilarity")
    
    dend_hub_linear %>% set("branches_k_color", value = c("red", "green", "blue"), k=3) %>% 
    set("labels_col", value = c("red", "green", "blue"), k=3) %>% 
    set("labels_cex", value = 0.7) %>% plot(main="k=3 using hub miRNA dissimilarity")
    
    dend_hub_nearest <- as.dendrogram(hclust_res_hub_nearest)
    dend_hub_nearest %>% set("branches_k_color", value = c("red", "blue"), k=2) %>% 
    set("labels_col", value = c("red", "blue"), k=2) %>% 
    set("labels_cex", value = 0.7) %>% plot(main="k=2 using hub miRNA dissimilarity")
    
    dend_hub_nearest %>% set("branches_k_color", value = c("red", "green", "blue"), k=3) %>% 
    set("labels_col", value = c("red", "green", "blue"), k=3) %>% 
    set("labels_cex", value = 0.7) %>% plot(main="k=3 using hub miRNA dissimilarity")

    dend_hub_spline <- as.dendrogram(hclust_res_hub_spline)
    dend_hub_spline %>% set("branches_k_color", value = c("red", "blue"), k=2) %>% 
    set("labels_col", value = c("red", "blue"), k=2) %>% 
    set("labels_cex", value = 0.7) %>% plot(main="k=2 using hub miRNA dissimilarity")
    
    dend_hub_spline %>% set("branches_k_color", value = c("red", "green", "blue"), k=3) %>% 
    set("labels_col", value = c("red", "green", "blue"), k=3) %>% 
    set("labels_cex", value = 0.7) %>% plot(main="k=3 using hub miRNA dissimilarity")
    
    dend_hub_cubic <- as.dendrogram(hclust_res_hub_cubic)
    dend_hub_cubic %>% set("branches_k_color", value = c("red", "blue"), k=2) %>% 
    set("labels_col", value = c("red", "blue"), k=2) %>% 
    set("labels_cex", value = 0.7) %>% plot(main="k=2 using hub miRNA dissimilarity")
    
    dend_hub_cubic %>% set("branches_k_color", value = c("red", "green", "blue"), k=3) %>% 
    set("labels_col", value = c("red", "green", "blue"), k=3) %>% 
    set("labels_cex", value = 0.7) %>% plot(main="k=3 using hub miRNA dissimilarity")
    
## Similarity network plot in terms of cell-specific miRNA-mRNA regulatory netowork
    col3 <- colorRampPalette(c("black", "blue", "red"))

    rownames(CsmiR_network_constant_Sim) <- colnames(CsmiR_network_constant_Sim) <- paste("Cell",c(1:19),sep=" ")
    corrplot.mixed(CsmiR_network_constant_Sim, lower = "number", upper = "pie", lower.col = col3(100), upper.col = col3(100), cl.lim = c(0, 1), number.cex = 0.8, tl.cex = 0.6)

    rownames(CsmiR_network_linear_Sim) <- colnames(CsmiR_network_linear_Sim) <- paste("Cell",c(1:19),sep=" ")
    corrplot.mixed(CsmiR_network_linear_Sim, lower = "number", upper = "pie", lower.col = col3(100), upper.col = col3(100), cl.lim = c(0, 1), number.cex = 0.8, tl.cex = 0.6)

    rownames(CsmiR_network_nearest_Sim) <- colnames(CsmiR_network_nearest_Sim) <- paste("Cell",c(1:19),sep=" ")
    corrplot.mixed(CsmiR_network_nearest_Sim, lower = "number", upper = "pie", lower.col = col3(100), upper.col = col3(100), cl.lim = c(0, 1), number.cex = 0.8, tl.cex = 0.6)

    rownames(CsmiR_network_spline_Sim) <- colnames(CsmiR_network_spline_Sim) <- paste("Cell",c(1:19),sep=" ")
    corrplot.mixed(CsmiR_network_spline_Sim, lower = "number", upper = "pie", lower.col = col3(100), upper.col = col3(100), cl.lim = c(0, 1), number.cex = 0.8, tl.cex = 0.6)

    rownames(CsmiR_network_cubic_Sim) <- colnames(CsmiR_network_cubic_Sim) <- paste("Cell",c(1:19),sep=" ")
    corrplot.mixed(CsmiR_network_cubic_Sim, lower = "number", upper = "pie", lower.col = col3(100), upper.col = col3(100), cl.lim = c(0, 1), number.cex = 0.8, tl.cex = 0.6)

## Similarity network plot in terms of cell-specific hub miRNAs
    rownames(CsmiR_hub_constant_Sim) <- colnames(CsmiR_hub_constant_Sim) <- paste("Cell",c(1:19),sep=" ")
    corrplot.mixed(CsmiR_hub_constant_Sim, lower = "number", upper = "pie", lower.col = col3(100), upper.col = col3(100), cl.lim = c(0, 1), number.cex = 0.8, tl.cex = 0.6)

    rownames(CsmiR_hub_linear_Sim) <- colnames(CsmiR_hub_linear_Sim) <- paste("Cell",c(1:19),sep=" ")
    corrplot.mixed(CsmiR_hub_linear_Sim, lower = "number", upper = "pie", lower.col = col3(100), upper.col = col3(100), cl.lim = c(0, 1), number.cex = 0.8, tl.cex = 0.6)

    rownames(CsmiR_hub_nearest_Sim) <- colnames(CsmiR_hub_nearest_Sim) <- paste("Cell",c(1:19),sep=" ")
    corrplot.mixed(CsmiR_hub_nearest_Sim, lower = "number", upper = "pie", lower.col = col3(100), upper.col = col3(100), cl.lim = c(0, 1), number.cex = 0.8, tl.cex = 0.6)

    rownames(CsmiR_hub_spline_Sim) <- colnames(CsmiR_hub_spline_Sim) <- paste("Cell",c(1:19),sep=" ")
    corrplot.mixed(CsmiR_hub_spline_Sim, lower = "number", upper = "pie", lower.col = col3(100), upper.col = col3(100), cl.lim = c(0, 1), number.cex = 0.8, tl.cex = 0.6)

    rownames(CsmiR_hub_cubic_Sim) <- colnames(CsmiR_hub_cubic_Sim) <- paste("Cell",c(1:19),sep=" ")
    corrplot.mixed(CsmiR_hub_cubic_Sim, lower = "number", upper = "pie", lower.col = col3(100), upper.col = col3(100), cl.lim = c(0, 1), number.cex = 0.8, tl.cex = 0.6)

## Violin plots
    index_net1 <- data.frame(value = unlist(lapply(seq(CSmiR_network_constant), function(i) nrow(CSmiR_network_constant[[i]]))), group = unlist(lapply(seq(CSmiR_network_constant), function(i) c("constant"))))
    index_net2 <- data.frame(value = unlist(lapply(seq(CSmiR_network_linear), function(i) nrow(CSmiR_network_linear[[i]]))), group = unlist(lapply(seq(CSmiR_network_linear), function(i) c("linear"))))
    index_net3 <- data.frame(value = unlist(lapply(seq(CSmiR_network_nearest), function(i) nrow(CSmiR_network_nearest[[i]]))), group = unlist(lapply(seq(CSmiR_network_nearest), function(i) c("nearest neighbor"))))
    index_net4 <- data.frame(value = unlist(lapply(seq(CSmiR_network_spline), function(i) nrow(CSmiR_network_spline[[i]]))), group = unlist(lapply(seq(CSmiR_network_spline), function(i) c("cubic spline"))))
    index_net5 <- data.frame(value = unlist(lapply(seq(CSmiR_network_cubic), function(i) nrow(CSmiR_network_cubic[[i]]))), group = unlist(lapply(seq(CSmiR_network_cubic), function(i) c("cubic Hermite"))))
    index_net <- rbind(index_net1, index_net2, index_net3, index_net4, index_net5)

p1 <- ggplot(index_net, aes(x = group, y = value, fill = group)) +
    geom_violin(trim=FALSE) + scale_y_log10() + 
    xlab("Pseudo-cells interpolation methods") +
    ylab("#Predicted cell-specific miRNA-mRNA interactions") +
    theme(legend.position = "none", 
          axis.text.x = element_text(face = "bold"),
	  axis.text.y = element_text(face = "bold"),
	  axis.title.x = element_text(face = "bold"),
	  axis.title.y = element_text(face = "bold")) +
    geom_boxplot(width = 0.1) +
    geom_signif(comparisons= list(c("constant", "linear"), c("constant", "nearest neighbor"), 
                                  c("constant", "cubic spline"), c("constant", "cubic Hermite"),
				  c("linear", "nearest neighbor"), c("linear", "cubic spline"),
				  c("linear", "cubic Hermite"),c("nearest neighbor", "cubic spline"),
				  c("nearest neighbor", "cubic Hermite"), c("cubic spline", "cubic Hermite")),
				  step_increase = 0.1)

    index_net1_validated <- data.frame(value = unlist(lapply(seq(CSmiR_network_constant_validated), function(i) nrow(CSmiR_network_constant_validated[[i]])/nrow(CSmiR_network_constant[[i]])*100)), group = unlist(lapply(seq(CSmiR_network_constant_validated), function(i) c("constant"))))
    index_net2_validated <- data.frame(value = unlist(lapply(seq(CSmiR_network_linear_validated), function(i) nrow(CSmiR_network_linear_validated[[i]])/nrow(CSmiR_network_linear[[i]])*100)), group = unlist(lapply(seq(CSmiR_network_linear_validated), function(i) c("linear"))))
    index_net3_validated <- data.frame(value = unlist(lapply(seq(CSmiR_network_nearest_validated), function(i) nrow(CSmiR_network_nearest_validated[[i]])/nrow(CSmiR_network_nearest[[i]])*100)), group = unlist(lapply(seq(CSmiR_network_nearest_validated), function(i) c("nearest neighbor"))))
    index_net4_validated <- data.frame(value = unlist(lapply(seq(CSmiR_network_spline_validated), function(i) nrow(CSmiR_network_spline_validated[[i]])/nrow(CSmiR_network_spline[[i]])*100)), group = unlist(lapply(seq(CSmiR_network_spline_validated), function(i) c("cubic spline"))))
    index_net5_validated <- data.frame(value = unlist(lapply(seq(CSmiR_network_cubic_validated), function(i) nrow(CSmiR_network_cubic_validated[[i]])/nrow(CSmiR_network_cubic[[i]])*100)), group = unlist(lapply(seq(CSmiR_network_cubic_validated), function(i) c("cubic Hermite"))))
    index_net_validated <- rbind(index_net1_validated, index_net2_validated, index_net3_validated, index_net4_validated, index_net5_validated)

p2 <- ggplot(index_net_validated, aes(x = group, y = value, fill = group)) +
    geom_violin(trim=FALSE) +  
    xlab("Pseudo-cells interpolation methods") +
    ylab("%Validated cell-specific miRNA-mRNA interactions") +    
    theme(legend.position = "none", 
          axis.text.x = element_text(face = "bold"),
	  axis.text.y = element_text(face = "bold"),
	  axis.title.x = element_text(face = "bold"),
	  axis.title.y = element_text(face = "bold")) +
    geom_boxplot(width = 0.1) +
    geom_signif(comparisons= list(c("constant", "linear"), c("constant", "nearest neighbor"), 
                                  c("constant", "cubic spline"), c("constant", "cubic Hermite"),
				  c("linear", "nearest neighbor"), c("linear", "cubic spline"),
				  c("linear", "cubic Hermite"),c("nearest neighbor", "cubic spline"),
				  c("nearest neighbor", "cubic Hermite"), c("cubic spline", "cubic Hermite")),
				  step_increase = 0.1)

    index_net1_CML <- data.frame(value = unlist(lapply(seq(CSmiR_network_constant_CML), function(i) nrow(CSmiR_network_constant_CML[[i]])/nrow(CSmiR_network_constant[[i]])*100)), group = unlist(lapply(seq(CSmiR_network_constant_CML), function(i) c("constant"))))
    index_net2_CML <- data.frame(value = unlist(lapply(seq(CSmiR_network_linear_CML), function(i) nrow(CSmiR_network_linear_CML[[i]])/nrow(CSmiR_network_linear[[i]])*100)), group = unlist(lapply(seq(CSmiR_network_linear_CML), function(i) c("linear"))))
    index_net3_CML <- data.frame(value = unlist(lapply(seq(CSmiR_network_nearest_CML), function(i) nrow(CSmiR_network_nearest_CML[[i]])/nrow(CSmiR_network_nearest[[i]])*100)), group = unlist(lapply(seq(CSmiR_network_nearest_CML), function(i) c("nearest neighbor"))))
    index_net4_CML <- data.frame(value = unlist(lapply(seq(CSmiR_network_spline_CML), function(i) nrow(CSmiR_network_spline_CML[[i]])/nrow(CSmiR_network_spline[[i]])*100)), group = unlist(lapply(seq(CSmiR_network_spline_CML), function(i) c("cubic spline"))))
    index_net5_CML <- data.frame(value = unlist(lapply(seq(CSmiR_network_cubic_CML), function(i) nrow(CSmiR_network_cubic_CML[[i]])/nrow(CSmiR_network_cubic[[i]])*100)), group = unlist(lapply(seq(CSmiR_network_cubic_CML), function(i) c("cubic Hermite"))))
    index_net_CML <- rbind(index_net1_CML, index_net2_CML, index_net3_CML, index_net4_CML, index_net5_CML)

p3 <- ggplot(index_net_CML, aes(x = group, y = value, fill = group)) +
    geom_violin(trim=FALSE) +  
    xlab("Pseudo-cells interpolation methods") +
    ylab("%CML-related cell-specific miRNA-mRNA interactions") +
    theme(legend.position = "none", 
          axis.text.x = element_text(face = "bold"),
	  axis.text.y = element_text(face = "bold"),
	  axis.title.x = element_text(face = "bold"),
	  axis.title.y = element_text(face = "bold")) +
    geom_boxplot(width = 0.1) +
    geom_signif(comparisons= list(c("constant", "linear"), c("constant", "nearest neighbor"), 
                                  c("constant", "cubic spline"), c("constant", "cubic Hermite"),
				  c("linear", "nearest neighbor"), c("linear", "cubic spline"),
				  c("linear", "cubic Hermite"),c("nearest neighbor", "cubic spline"),
				  c("nearest neighbor", "cubic Hermite"), c("cubic spline", "cubic Hermite")),
				  step_increase = 0.1)

    index_hub1_CML <- data.frame(value = unlist(lapply(seq(hub_miRNAs_constant_CML), function(i) length(hub_miRNAs_constant_CML[[i]])/length(hub_miRNAs_constant[[i]])*100)), group = unlist(lapply(seq(hub_miRNAs_constant_CML), function(i) c("constant"))))
    index_hub2_CML <- data.frame(value = unlist(lapply(seq(hub_miRNAs_linear_CML), function(i) length(hub_miRNAs_linear_CML[[i]])/length(hub_miRNAs_linear[[i]])*100)), group = unlist(lapply(seq(hub_miRNAs_linear_CML), function(i) c("linear"))))
    index_hub3_CML <- data.frame(value = unlist(lapply(seq(hub_miRNAs_nearest_CML), function(i) length(hub_miRNAs_nearest_CML[[i]])/length(hub_miRNAs_nearest[[i]])*100)), group = unlist(lapply(seq(hub_miRNAs_nearest_CML), function(i) c("nearest neighbor"))))
    index_hub4_CML <- data.frame(value = unlist(lapply(seq(hub_miRNAs_spline_CML), function(i) length(hub_miRNAs_spline_CML[[i]])/length(hub_miRNAs_spline[[i]])*100)), group = unlist(lapply(seq(hub_miRNAs_spline_CML), function(i) c("cubic spline"))))
    index_hub5_CML <- data.frame(value = unlist(lapply(seq(hub_miRNAs_cubic_CML), function(i) length(hub_miRNAs_cubic_CML[[i]])/length(hub_miRNAs_cubic[[i]])*100)), group = unlist(lapply(seq(hub_miRNAs_cubic_CML), function(i) c("cubic Hermite"))))
    index_hub_CML <- rbind(index_hub1_CML, index_hub2_CML, index_hub3_CML, index_hub4_CML, index_hub5_CML)

p4 <- ggplot(index_hub_CML, aes(x = group, y = value, fill = group)) +
    geom_violin(trim=FALSE) +  
    xlab("Pseudo-cells interpolation methods") +
    ylab("%CML-related cell-specific hub miRNAs") +
    theme(legend.position = "none", 
          axis.text.x = element_text(face = "bold"),
	  axis.text.y = element_text(face = "bold"),
	  axis.title.x = element_text(face = "bold"),
	  axis.title.y = element_text(face = "bold")) +
    geom_boxplot(width = 0.1) +
    geom_signif(comparisons= list(c("constant", "linear"), c("constant", "nearest neighbor"), 
                                  c("constant", "cubic spline"), c("constant", "cubic Hermite"),
				  c("linear", "nearest neighbor"), c("linear", "cubic spline"),
				  c("linear", "cubic Hermite"),c("nearest neighbor", "cubic spline"),
				  c("nearest neighbor", "cubic Hermite"), c("cubic spline", "cubic Hermite")),
				  step_increase = 0.1)

    index_overlap1_validated <- data.frame(value = c(nrow(Overlap_network_constant_conserved_validated)/nrow(Overlap_network_constant)*100, nrow(Overlap_network_constant_rewired_validated)/nrow(Overlap_network_constant_rewired)*100), group = c("constant", "constant"), type = c("conserved", "rewired"))
    index_overlap2_validated <- data.frame(value = c(nrow(Overlap_network_linear_conserved_validated)/nrow(Overlap_network_linear)*100, nrow(Overlap_network_linear_rewired_validated)/nrow(Overlap_network_linear_rewired)*100), group = c("linear", "linear"), type = c("conserved", "rewired"))
    index_overlap3_validated <- data.frame(value = c(nrow(Overlap_network_nearest_conserved_validated)/nrow(Overlap_network_nearest)*100, nrow(Overlap_network_nearest_rewired_validated)/nrow(Overlap_network_nearest_rewired)*100), group = c("nearest neighbor", "nearest neighbor"), type = c("conserved", "rewired"))
    index_overlap4_validated <- data.frame(value = c(nrow(Overlap_network_spline_conserved_validated)/nrow(Overlap_network_spline)*100, nrow(Overlap_network_spline_rewired_validated)/nrow(Overlap_network_spline_rewired)*100), group = c("cubic spline", "cubic spline"), type = c("conserved", "rewired"))
    index_overlap5_validated <- data.frame(value = c(nrow(Overlap_network_cubic_conserved_validated)/nrow(Overlap_network_cubic)*100, nrow(Overlap_network_cubic_rewired_validated)/nrow(Overlap_network_cubic_rewired)*100), group = c("cubic Hermite", "cubic Hermite"), type = c("conserved", "rewired"))
    index_overlap_validated <- rbind(index_overlap1_validated, index_overlap2_validated, index_overlap3_validated, index_overlap4_validated, index_overlap5_validated)

p5 <- ggplot(index_overlap_validated, aes(x = group, y = value, fill = type)) +
    geom_bar(stat = 'identity', position = 'stack') +  
    xlab("Pseudo-cells interpolation methods") +
    ylab("%Validated miRNA-mRNA interactions") +
    theme(legend.position = c(0.8, 0.9),          
          axis.text.x = element_text(face = "bold"),
	  axis.text.y = element_text(face = "bold"),
	  axis.title.x = element_text(face = "bold"),
	  axis.title.y = element_text(face = "bold")) + 
    guides(fill=guide_legend(title = "miRNA regulation type")) 

    index_overlap1_CML <- data.frame(value = c(nrow(Overlap_network_constant_CML)/nrow(Overlap_network_constant)*100, nrow(Overlap_network_constant_rewired_CML)/nrow(Overlap_network_constant_rewired)*100), group = c("constant", "constant"), type = c("conserved", "rewired"))
    index_overlap2_CML <- data.frame(value = c(nrow(Overlap_network_linear_CML)/nrow(Overlap_network_linear)*100, nrow(Overlap_network_linear_rewired_CML)/nrow(Overlap_network_linear_rewired)*100), group = c("linear", "linear"), type = c("conserved", "rewired"))
    index_overlap3_CML <- data.frame(value = c(nrow(Overlap_network_nearest_CML)/nrow(Overlap_network_nearest)*100, nrow(Overlap_network_nearest_rewired_CML)/nrow(Overlap_network_nearest_rewired)*100), group = c("nearest neighbor", "nearest neighbor"), type = c("conserved", "rewired"))
    index_overlap4_CML <- data.frame(value = c(nrow(Overlap_network_spline_CML)/nrow(Overlap_network_spline)*100, nrow(Overlap_network_spline_rewired_CML)/nrow(Overlap_network_spline_rewired)*100), group = c("cubic spline", "cubic spline"), type = c("conserved", "rewired"))
    index_overlap5_CML <- data.frame(value = c(nrow(Overlap_network_cubic_CML)/nrow(Overlap_network_cubic)*100, nrow(Overlap_network_cubic_rewired_CML)/nrow(Overlap_network_cubic_rewired)*100), group = c("cubic Hermite", "cubic Hermite"), type = c("conserved", "rewired"))
    index_overlap_CML <- rbind(index_overlap1_CML, index_overlap2_CML, index_overlap3_CML, index_overlap4_CML, index_overlap5_CML)

p6 <- ggplot(index_overlap_CML, aes(x = group, y = value, fill = type)) +
    geom_bar(stat = 'identity', position = 'stack') +  
    xlab("Pseudo-cells interpolation methods") +
    ylab("%CML-related miRNA-mRNA interactions") +
    theme(legend.position = c(0.8, 0.9),          
          axis.text.x = element_text(face = "bold"),
	  axis.text.y = element_text(face = "bold"),
	  axis.title.x = element_text(face = "bold"),
	  axis.title.y = element_text(face = "bold")) + 
    guides(fill=guide_legend(title = "miRNA regulation type")) +
    ylim(0, 1.5)

    index_overlap1_hub_CML <- data.frame(value = c(length(Overlap_hub_miRNAs_constant_conserved_CML)/length(Overlap_hub_miRNAs_constant_conserved)*100, length(Overlap_network_constant_rewired_CML)/length(Overlap_hub_miRNAs_constant_rewired)*100), group = c("constant", "constant"), type = c("conserved", "rewired"))
    index_overlap2_hub_CML <- data.frame(value = c(length(Overlap_hub_miRNAs_linear_conserved_CML)/length(Overlap_hub_miRNAs_linear_conserved)*100, length(Overlap_network_linear_rewired_CML)/length(Overlap_hub_miRNAs_linear_rewired)*100), group = c("linear", "linear"), type = c("conserved", "rewired"))
    index_overlap3_hub_CML <- data.frame(value = c(length(Overlap_hub_miRNAs_nearest_conserved_CML)/length(Overlap_hub_miRNAs_nearest_conserved)*100, length(Overlap_network_nearest_rewired_CML)/length(Overlap_hub_miRNAs_nearest_rewired)*100), group = c("nearest neighbor", "nearest neighbor"), type = c("conserved", "rewired"))
    index_overlap4_hub_CML <- data.frame(value = c(length(Overlap_hub_miRNAs_spline_conserved_CML)/length(Overlap_hub_miRNAs_spline_conserved)*100, length(Overlap_network_spline_rewired_CML)/length(Overlap_hub_miRNAs_spline_rewired)*100), group = c("cubic spline", "cubic spline"), type = c("conserved", "rewired"))
    index_overlap5_hub_CML <- data.frame(value = c(length(Overlap_hub_miRNAs_cubic_conserved_CML)/length(Overlap_hub_miRNAs_cubic_conserved)*100, length(Overlap_network_cubic_rewired_CML)/length(Overlap_hub_miRNAs_cubic_rewired)*100), group = c("cubic Hermite", "cubic Hermite"), type = c("conserved", "rewired"))
    index_overlap_hub_CML <- rbind(index_overlap1_hub_CML, index_overlap2_hub_CML, index_overlap3_hub_CML, index_overlap4_hub_CML, index_overlap5_hub_CML)

p7 <- ggplot(index_overlap_hub_CML, aes(x = group, y = value, fill = type)) +
    geom_bar(stat = 'identity', position = 'stack') +  
    xlab("Pseudo-cells interpolation methods") +
    ylab("%CML-related hub miRNAs") +
    theme(legend.position = c(0.2, 0.9),          
          axis.text.x = element_text(face = "bold"),
	  axis.text.y = element_text(face = "bold"),
	  axis.title.x = element_text(face = "bold"),
	  axis.title.y = element_text(face = "bold")) + 
    guides(fill=guide_legend(title = "miRNA regulation type"))
   

    index_miRfamily_net1 <- data.frame(value = unlist(lapply(seq(CSmiR_network_constant_miRfamily), function(i) nrow(CSmiR_network_constant_miRfamily[[i]]))), group = unlist(lapply(seq(CSmiR_network_constant_miRfamily), function(i) c("constant"))), id = seq(19))
    index_miRfamily_net2 <- data.frame(value = unlist(lapply(seq(CSmiR_network_linear_miRfamily), function(i) nrow(CSmiR_network_linear_miRfamily[[i]]))), group = unlist(lapply(seq(CSmiR_network_linear_miRfamily), function(i) c("linear"))), id = seq(19))
    index_miRfamily_net3 <- data.frame(value = unlist(lapply(seq(CSmiR_network_nearest_miRfamily), function(i) nrow(CSmiR_network_nearest_miRfamily[[i]]))), group = unlist(lapply(seq(CSmiR_network_nearest_miRfamily), function(i) c("nearest neighbor"))), id = seq(19))
    index_miRfamily_net4 <- data.frame(value = unlist(lapply(seq(CSmiR_network_spline_miRfamily), function(i) nrow(CSmiR_network_spline_miRfamily[[i]]))), group = unlist(lapply(seq(CSmiR_network_spline_miRfamily), function(i) c("cubic spline"))), id = seq(19))
    index_miRfamily_net5 <- data.frame(value = unlist(lapply(seq(CSmiR_network_cubic_miRfamily), function(i) nrow(CSmiR_network_cubic_miRfamily[[i]]))), group = unlist(lapply(seq(CSmiR_network_cubic_miRfamily), function(i) c("cubic Hermite"))), id = seq(19))
    index_miRfamily_net <- rbind(index_miRfamily_net1, index_miRfamily_net2, index_miRfamily_net3, index_miRfamily_net4, index_miRfamily_net5)

p8 <- ggplot(index_miRfamily_net, aes(x = id, y = value, fill = group)) +
    geom_area() +
    xlab("Single-cell ID") +
    ylab("#Predicted targets of miR-17, miR-20a and miR-92a-1") +
    theme(legend.position = c(0.85, 0.83),          
          axis.text.x = element_text(face = "bold"),
	  axis.text.y = element_text(face = "bold"),
	  axis.title.x = element_text(face = "bold"),
	  axis.title.y = element_text(face = "bold")) + 
    guides(fill=guide_legend(title = "Interpolation method")) +
    scale_x_continuous(breaks = seq(1, 19, 1))

    index_miRfamily_validated_net1 <- data.frame(value = unlist(lapply(seq(CSmiR_network_constant_miRfamily_validated), function(i) nrow(CSmiR_network_constant_miRfamily_validated[[i]])/nrow(CSmiR_network_constant_miRfamily[[i]])*100)), group = unlist(lapply(seq(CSmiR_network_constant_miRfamily_validated), function(i) c("constant"))), id = seq(19))
    index_miRfamily_validated_net2 <- data.frame(value = unlist(lapply(seq(CSmiR_network_linear_miRfamily_validated), function(i) nrow(CSmiR_network_linear_miRfamily_validated[[i]])/nrow(CSmiR_network_linear_miRfamily[[i]])*100)), group = unlist(lapply(seq(CSmiR_network_linear_miRfamily_validated), function(i) c("linear"))), id = seq(19))
    index_miRfamily_validated_net3 <- data.frame(value = unlist(lapply(seq(CSmiR_network_nearest_miRfamily_validated), function(i) nrow(CSmiR_network_nearest_miRfamily_validated[[i]])/nrow(CSmiR_network_nearest_miRfamily[[i]])*100)), group = unlist(lapply(seq(CSmiR_network_nearest_miRfamily_validated), function(i) c("nearest neighbor"))), id = seq(19))
    index_miRfamily_validated_net4 <- data.frame(value = unlist(lapply(seq(CSmiR_network_spline_miRfamily_validated), function(i) nrow(CSmiR_network_spline_miRfamily_validated[[i]])/nrow(CSmiR_network_spline_miRfamily[[i]])*100)), group = unlist(lapply(seq(CSmiR_network_spline_miRfamily_validated), function(i) c("cubic spline"))), id = seq(19))
    index_miRfamily_validated_net5 <- data.frame(value = unlist(lapply(seq(CSmiR_network_cubic_miRfamily_validated), function(i) nrow(CSmiR_network_cubic_miRfamily_validated[[i]])/nrow(CSmiR_network_cubic_miRfamily[[i]])*100)), group = unlist(lapply(seq(CSmiR_network_cubic_miRfamily_validated), function(i) c("cubic Hermite"))), id = seq(19))
    index_miRfamily_validated_net <- rbind(index_miRfamily_validated_net1, index_miRfamily_validated_net2, index_miRfamily_validated_net3, index_miRfamily_validated_net4, index_miRfamily_validated_net5)

p9 <- ggplot(index_miRfamily_validated_net, aes(x = id, y = value, fill = group)) +
    geom_area() +
    xlab("Single-cell ID") +
    ylab("%Validated targets of miR-17, miR-20a and miR-92a-1") +
    theme(legend.position = c(0.75, 0.83),          
          axis.text.x = element_text(face = "bold"),
	  axis.text.y = element_text(face = "bold"),
	  axis.title.x = element_text(face = "bold"),
	  axis.title.y = element_text(face = "bold")) + 
    guides(fill=guide_legend(title = "Interpolation method")) +
    scale_x_continuous(breaks = seq(1, 19, 1))

    index_miRfamily_CML_net1 <- data.frame(value = unlist(lapply(seq(CSmiR_network_constant_miRfamily_CML), function(i) nrow(CSmiR_network_constant_miRfamily_CML[[i]])/nrow(CSmiR_network_constant_miRfamily[[i]])*100)), group = unlist(lapply(seq(CSmiR_network_constant_miRfamily_CML), function(i) c("constant"))), id = seq(19))
    index_miRfamily_CML_net2 <- data.frame(value = unlist(lapply(seq(CSmiR_network_linear_miRfamily_CML), function(i) nrow(CSmiR_network_linear_miRfamily_CML[[i]])/nrow(CSmiR_network_linear_miRfamily[[i]])*100)), group = unlist(lapply(seq(CSmiR_network_linear_miRfamily_CML), function(i) c("linear"))), id = seq(19))
    index_miRfamily_CML_net3 <- data.frame(value = unlist(lapply(seq(CSmiR_network_nearest_miRfamily_CML), function(i) nrow(CSmiR_network_nearest_miRfamily_CML[[i]])/nrow(CSmiR_network_nearest_miRfamily[[i]])*100)), group = unlist(lapply(seq(CSmiR_network_nearest_miRfamily_CML), function(i) c("nearest neighbor"))), id = seq(19))
    index_miRfamily_CML_net4 <- data.frame(value = unlist(lapply(seq(CSmiR_network_spline_miRfamily_CML), function(i) nrow(CSmiR_network_spline_miRfamily_CML[[i]])/nrow(CSmiR_network_spline_miRfamily[[i]])*100)), group = unlist(lapply(seq(CSmiR_network_spline_miRfamily_CML), function(i) c("cubic spline"))), id = seq(19))
    index_miRfamily_CML_net5 <- data.frame(value = unlist(lapply(seq(CSmiR_network_cubic_miRfamily_CML), function(i) nrow(CSmiR_network_cubic_miRfamily_CML[[i]])/nrow(CSmiR_network_cubic_miRfamily[[i]])*100)), group = unlist(lapply(seq(CSmiR_network_cubic_miRfamily_CML), function(i) c("cubic Hermite"))), id = seq(19))
    index_miRfamily_CML_net <- rbind(index_miRfamily_CML_net1, index_miRfamily_CML_net2, index_miRfamily_CML_net3, index_miRfamily_CML_net4, index_miRfamily_CML_net5)

p10 <- ggplot(index_miRfamily_CML_net, aes(x = id, y = value, fill = group)) +
    geom_area() +
    xlab("Single-cell ID") +
    ylab("%CML-related targets of miR-17, miR-20a and miR-92a-1") +
    theme(legend.position = c(0.83, 0.83),          
          axis.text.x = element_text(face = "bold"),
	  axis.text.y = element_text(face = "bold"),
	  axis.title.x = element_text(face = "bold"),
	  axis.title.y = element_text(face = "bold")) + 
    guides(fill=guide_legend(title = "Interpolation method")) +
    scale_x_continuous(breaks = seq(1, 19, 1)) +
    ylim(0, 13)

    index_miRfamily_cr_net1 <- data.frame(value = c(nrow(Overlap_network_constant_miRfamily), nrow(Overlap_network_constant_rewired_miRfamily)), group = c("constant", "constant"), type = c("conserved", "rewired"))
    index_miRfamily_cr_net2 <- data.frame(value = c(nrow(Overlap_network_linear_miRfamily), nrow(Overlap_network_linear_rewired_miRfamily)), group = c("linear", "linear"), type = c("conserved", "rewired"))
    index_miRfamily_cr_net3 <- data.frame(value = c(nrow(Overlap_network_nearest_miRfamily), nrow(Overlap_network_nearest_rewired_miRfamily)), group = c("nearest neighbor", "nearest neighbor"), type = c("conserved", "rewired"))
    index_miRfamily_cr_net4 <- data.frame(value = c(nrow(Overlap_network_spline_miRfamily), nrow(Overlap_network_spline_rewired_miRfamily)), group = c("cubic spline", "cubic spline"), type = c("conserved", "rewired"))
    index_miRfamily_cr_net5 <- data.frame(value = c(nrow(Overlap_network_cubic_miRfamily), nrow(Overlap_network_cubic_rewired_miRfamily)), group = c("cubic Hermite", "cubic Hermite"), type = c("conserved", "rewired"))
    index_miRfamily_cr_net <- rbind(index_miRfamily_cr_net1, index_miRfamily_cr_net2, index_miRfamily_cr_net3, index_miRfamily_cr_net4, index_miRfamily_cr_net5)

p11 <- ggplot(index_miRfamily_cr_net, aes(x = type, y = value, fill = group)) +
    geom_bar(stat = 'identity', position = 'dodge') +
    xlab("miRNA regulation type") +
    ylab("#Predicted targets of miR-17, miR-20a and miR-92a-1") +
    theme(legend.position = c(0.2, 0.8),          
          axis.text.x = element_text(face = "bold"),
	  axis.text.y = element_text(face = "bold"),
	  axis.title.x = element_text(face = "bold"),
	  axis.title.y = element_text(face = "bold")) + 
    guides(fill=guide_legend(title = "Interpolation method")) +
    geom_text(mapping = aes(label = value), colour = 'black', vjust = -.5, hjust = .5, position = position_dodge(0.9)) +
    annotate("text", x = 1.825, y = -5, label="*", colour = '#7CAE00', fontface = "bold", size = 5) + 
    annotate("text", x = 2.365, y = -5, label="*", colour = '#E76BF3', fontface = "bold", size = 5)

save.image("CSmiR.RData")

