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
    CSmiR_network_bootstrap <- CSmiR_net_bootstrap(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, 
                                                   boxsize = 0.1, bootstrap_betw_point = 5, 
						   bootstrap_num = 100, p.value.cutoff = 0.01)
  
## Experimentally validated cell-specific miRNA-mRNA interactions, the ground-truth (miRTarget variable) is from miRTarBase v8.0 and TarBase v8.0
    miRTarget_graph <- make_graph(c(t(miRTarget[, 1:2])), directed = TRUE)
    
    CSmiR_network_bootstrap_graph <- lapply(seq(CSmiR_network_bootstrap), function(i) make_graph(c(t(CSmiR_network_bootstrap[[i]][, 1:2])), directed = TRUE))
    CSmiR_network_bootstrap_validated <- lapply(seq(CSmiR_network_bootstrap), function(i) as_data_frame(CSmiR_network_bootstrap_graph[[i]] %s% miRTarget_graph))

## CML-related cell-specific miRNA-mRNA interactions. the list of CML-related miRNAs and mRNAs (CML variable) is from HMDD v3.0 and DisGeNET v5.0, respectively    
    CSmiR_network_bootstrap_CML <- lapply(seq(CSmiR_network_bootstrap), function(i) CSmiR_network_bootstrap[[i]][intersect(which(CSmiR_network_bootstrap[[i]][, 1] %in% as.matrix(CML)), which(CSmiR_network_bootstrap[[i]][, 2] %in% as.matrix(CML))), ])

## Overlap of cell-specific miRNA-mRNA interactions across cells    
    Overlap_network_bootstrap <- Overlap_net(CSmiR_network_bootstrap, Intersect_num = round(length(CSmiR_network_bootstrap)*0.9))
    Overlap_network_bootstrap_union <- Overlap_net(CSmiR_network_bootstrap, Intersect_num = 1)
    Overlap_network_bootstrap_two <- Overlap_net(CSmiR_network_bootstrap, Intersect_num = 2)
    Overlap_network_bootstrap_union_graph <- make_graph(c(t(Overlap_network_bootstrap_union[, 1:2])), directed = TRUE)
    Overlap_network_bootstrap_two_graph <- make_graph(c(t(Overlap_network_bootstrap_two[, 1:2])), directed = TRUE)
    Overlap_network_bootstrap_rewired <- as_data_frame(Overlap_network_bootstrap_union_graph %m% Overlap_network_bootstrap_two_graph)
    Overlap_network_bootstrap_graph <- make_graph(c(t(Overlap_network_bootstrap[, 1:2])), directed = TRUE)
    Overlap_network_bootstrap_conserved_validated <- as_data_frame(Overlap_network_bootstrap_graph %s% miRTarget_graph)
    Overlap_network_bootstrap_rewired_graph <- make_graph(c(t(Overlap_network_bootstrap_rewired[, 1:2])), directed = TRUE)
    Overlap_network_bootstrap_rewired_validated <- as_data_frame(Overlap_network_bootstrap_rewired_graph %s% miRTarget_graph)

## CML-related conserved and rewired miRNA-mRNA regulatory network    
    Overlap_network_bootstrap_CML <- Overlap_network_bootstrap[intersect(which(Overlap_network_bootstrap[, 1] %in% as.matrix(CML)), which(Overlap_network_bootstrap[, 2] %in% as.matrix(CML))), ]
    Overlap_network_bootstrap_rewired_CML <- Overlap_network_bootstrap_rewired[intersect(which(Overlap_network_bootstrap_rewired[, 1] %in% as.matrix(CML)), which(Overlap_network_bootstrap_rewired[, 2] %in% as.matrix(CML))), ]

## Identifying cell-specific hub miRNAs    
    CSmiR_network_bootstrap_outdegree <- lapply(seq(CSmiR_network_bootstrap), function(i) degree(CSmiR_network_bootstrap_graph[[i]], mode="out"))
    hub_miRNAs_bootstrap <- lapply(seq(CSmiR_network_bootstrap), function(i) names(sort(CSmiR_network_bootstrap_outdegree[[i]][which(CSmiR_network_bootstrap_outdegree[[i]]!=0)], decreasing=TRUE))[1:ceiling(0.2*length(which(CSmiR_network_bootstrap_outdegree[[i]]!=0)))])
 
## Overlap of cell-specific hub miRNAs across cells    
    Overlap_hub_miRNAs_bootstrap_conserved <- Overlap_hub(hub_miRNAs_bootstrap, Intersect_num = round(length(hub_miRNAs_bootstrap)*0.9))
    Overlap_hub_miRNAs_bootstrap_union <- Overlap_hub(hub_miRNAs_bootstrap, Intersect_num = 1)
    Overlap_hub_miRNAs_bootstrap_two <- Overlap_hub(hub_miRNAs_bootstrap, Intersect_num = 2)
    Overlap_hub_miRNAs_bootstrap_rewired <- setdiff(Overlap_hub_miRNAs_bootstrap_union, Overlap_hub_miRNAs_bootstrap_two)

## CML-related cell-specific hub miRNAs    
    hub_miRNAs_bootstrap_CML <- lapply(seq(hub_miRNAs_bootstrap), function(i) hub_miRNAs_bootstrap[[i]][which(hub_miRNAs_bootstrap[[i]] %in% as.matrix(CML))])

## CML-related conserved and rewired hub miRNAs    
    Overlap_hub_miRNAs_bootstrap_conserved_CML <- Overlap_hub_miRNAs_bootstrap_conserved[which(Overlap_hub_miRNAs_bootstrap_conserved %in% as.matrix(CML))]    
    Overlap_hub_miRNAs_bootstrap_rewired_CML <- Overlap_hub_miRNAs_bootstrap_rewired[which(Overlap_hub_miRNAs_bootstrap_rewired %in% as.matrix(CML))]

## Calculating similarity matrix of cell-specific miRNA-mRNA regulatory network across cells    
    CSmiR_network_bootstrap_Sim <- Sim.network(CSmiR_network_bootstrap, CSmiR_network_bootstrap, directed = TRUE)

## Calculating similarity matrix of cell-specific hub miRNAs across cells    
    CSmiR_hub_bootstrap_Sim <- Sim.hub(hub_miRNAs_bootstrap, hub_miRNAs_bootstrap)

## Identifying cell-cell crosstalk network in terms of network similarity matrix    
    CSmiR_network_bootstrap_adjacency_matrix <- ifelse(CSmiR_network_bootstrap_Sim > median(CSmiR_network_bootstrap_Sim[lower.tri(CSmiR_network_bootstrap_Sim)]), 1, 0)
    diag(CSmiR_network_bootstrap_adjacency_matrix) <- 0
    CSmiR_network_bootstrap_adjacency_matrix_graph <- graph_from_adjacency_matrix(CSmiR_network_bootstrap_adjacency_matrix, mode = "undirected")

## Identifying cell-cell crosstalk network in terms of hub miRNA similarity matrix    
    CSmiR_hub_bootstrap_adjacency_matrix <- ifelse(CSmiR_hub_bootstrap_Sim > median(CSmiR_hub_bootstrap_Sim[lower.tri(CSmiR_hub_bootstrap_Sim)]), 1, 0)
    diag(CSmiR_hub_bootstrap_adjacency_matrix) <- 0
    CSmiR_hub_bootstrap_adjacency_matrix_graph <- graph_from_adjacency_matrix(CSmiR_hub_bootstrap_adjacency_matrix, mode = "undirected")

## Identifying hub cells in terms of network similarity    
    CSmiR_network_bootstrap_cell_degree <- degree(CSmiR_network_bootstrap_adjacency_matrix_graph)
    CSmiR_network_bootstrap_hub_cells <- names(sort(CSmiR_network_bootstrap_cell_degree[which(CSmiR_network_bootstrap_cell_degree!=0)], decreasing=TRUE))[1:ceiling(0.2*length(which(CSmiR_network_bootstrap_cell_degree!=0)))]

## Identifying hub cells in terms of hub miRNA similarity    
    CSmiR_hub_bootstrap_cell_degree <- degree(CSmiR_hub_bootstrap_adjacency_matrix_graph)
    CSmiR_hub_bootstrap_hub_cells <- names(sort(CSmiR_hub_bootstrap_cell_degree[which(CSmiR_hub_bootstrap_cell_degree!=0)], decreasing=TRUE))[1:ceiling(0.2*length(which(CSmiR_hub_bootstrap_cell_degree!=0)))]

## Identifying cell-cell crosstalk modules in terms of network similarity matrix   
    CSmiR_network_bootstrap_cell_module <- netModule(CSmiR_network_bootstrap_adjacency_matrix_graph %>% as_data_frame)

## Identifying cell-cell crosstalk modules in terms of hub miRNA similarity matrix    
    CSmiR_hub_bootstrap_cell_module <- netModule(CSmiR_hub_bootstrap_adjacency_matrix_graph %>% as_data_frame)

## The miR-17/92 family analysis
    miR_17_92_family <- c("hsa-miR-17-3p", "hsa-miR-17-5p","hsa-miR-18a-3p", "hsa-miR-18a-5p", "hsa-miR-19a-3p", "hsa-miR-19a-5p",  
                          "hsa-miR-19b-3p", "hsa-miR-19b-1-5p", "hsa-miR-20a-3p", "hsa-miR-20a-5p", "hsa-miR-92a-3p", "hsa-miR-92a-1-5p")
    
    # Extracting miR-17/92 family related cell-specific miRNA-mRNA regulatory networks
    CSmiR_network_bootstrap_miRfamily <- lapply(seq(CSmiR_network_bootstrap), function(i) CSmiR_network_bootstrap[[i]][which(CSmiR_network_bootstrap[[i]][, 1] %in% as.matrix(miR_17_92_family)), ])
 
    # Extracting conserved and rewired miRNA-mRNA regulatory networks associated with the miR-17/92 family 
    Overlap_network_bootstrap_miRfamily <- Overlap_network_bootstrap[which(Overlap_network_bootstrap[, 1] %in% as.matrix(miR_17_92_family)), ]
    Overlap_network_bootstrap_rewired_miRfamily <- Overlap_network_bootstrap_rewired[which(Overlap_network_bootstrap_rewired[, 1] %in% as.matrix(miR_17_92_family)), ]
    
    # Extracting CML-related conserved and rewired miRNA-mRNA regulatory networks associated with the miR-17/92 family
    CSmiR_network_bootstrap_miRfamily_CML <- lapply(seq(CSmiR_network_bootstrap_miRfamily), function(i) CSmiR_network_bootstrap_miRfamily[[i]][intersect(which(CSmiR_network_bootstrap_miRfamily[[i]][, 1] %in% as.matrix(CML)), which(CSmiR_network_bootstrap_miRfamily[[i]][, 2] %in% as.matrix(CML))), ])
    
    # Experimentally validated conserved and rewired miRNA-mRNA interactions associated with the miR-17/92 family
    CSmiR_network_bootstrap_miRfamily_graph <- lapply(seq(CSmiR_network_bootstrap_miRfamily), function(i) make_graph(c(t(CSmiR_network_bootstrap_miRfamily[[i]][, 1:2])), directed = TRUE))
    CSmiR_network_bootstrap_miRfamily_validated <- lapply(seq(CSmiR_network_bootstrap_miRfamily), function(i) as_data_frame(CSmiR_network_bootstrap_miRfamily_graph[[i]] %s% miRTarget_graph))

## Discovering maximal bicliques of conserved and rewired the miR-17/92 family regulation   
    CSmiR_network_miRfamily_bootstrap_biclique <- biclique_network(list(Sub_miR(Overlap_network_bootstrap_miRfamily), Sub_miR(Overlap_network_bootstrap_rewired_miRfamily)), lleast = 2, rleast = 3)

## Enrichment analysis of maximal bicliques of conserved and rewired the miR-17/92 family regulation   
   miRfamily_bootstrap_biclique_conserved_gene_list <- lapply(seq(CSmiR_network_miRfamily_bootstrap_biclique[[1]]), function(i) CSmiR_network_miRfamily_bootstrap_biclique[[1]][[i]]$right)
   miRfamily_bootstrap_biclique_rewired_gene_list <- lapply(seq(CSmiR_network_miRfamily_bootstrap_biclique[[2]]), function(i) CSmiR_network_miRfamily_bootstrap_biclique[[2]][[i]]$right)

   # GO, KEGG and Reactome enrichment analysis
   miRfamily_bootstrap_biclique_conserved_FEA <- moduleFEA(miRfamily_bootstrap_biclique_conserved_gene_list, padjustvaluecutoff = 0.05)
   miRfamily_bootstrap_biclique_rewired_FEA <- moduleFEA(miRfamily_bootstrap_biclique_rewired_gene_list, padjustvaluecutoff = 0.05)

   miRfamily_bootstrap_biclique_conserved_miRmR_list <- lapply(seq(CSmiR_network_miRfamily_bootstrap_biclique[[1]]), function(i) c(CSmiR_network_miRfamily_bootstrap_biclique[[1]][[i]]$left, CSmiR_network_miRfamily_bootstrap_biclique[[1]][[i]]$right))
   miRfamily_bootstrap_biclique_rewired_miRmR_list <- lapply(seq(CSmiR_network_miRfamily_bootstrap_biclique[[2]]), function(i) c(CSmiR_network_miRfamily_bootstrap_biclique[[2]][[i]]$left, CSmiR_network_miRfamily_bootstrap_biclique[[2]][[i]]$right))

   # CML enrichment analysis   
   miRfamily_bootstrap_biclique_conserved_CML_EA <- module_CML_EA(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, CML, miRfamily_bootstrap_biclique_conserved_miRmR_list)
   miRfamily_bootstrap_biclique_rewired_CML_EA <- module_CML_EA(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, CML, miRfamily_bootstrap_biclique_rewired_miRmR_list)

   # Hallmark enrichment analysis
   m_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, human_gene_symbol)
   
   miRfamily_bootstrap_biclique_conserved_Hallmark <- lapply(seq(miRfamily_bootstrap_biclique_conserved_gene_list), function(i) enricher(miRfamily_bootstrap_biclique_conserved_gene_list[[i]], TERM2GENE=m_t2g, minGSSize=1) %>% as.data.frame)
   miRfamily_bootstrap_biclique_rewired_Hallmark <- lapply(seq(miRfamily_bootstrap_biclique_rewired_gene_list), function(i) enricher(miRfamily_bootstrap_biclique_rewired_gene_list[[i]], TERM2GENE=m_t2g, minGSSize=1) %>% as.data.frame)

   # Cell marker enrichment analsyis
   cell_markers <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt') %>%
   tidyr::unite("cellMarker", tissueType, cancerType, cellName, sep=", ") %>% 
   dplyr::select(cellMarker, geneSymbol)    
   
   miRfamily_bootstrap_biclique_conserved_Cellmarker <- lapply(seq(miRfamily_bootstrap_biclique_conserved_gene_list), function(i) enricher(miRfamily_bootstrap_biclique_conserved_gene_list[[i]], TERM2GENE=cell_markers, minGSSize=1) %>% as.data.frame)
   miRfamily_bootstrap_biclique_rewired_Cellmarker <- lapply(seq(miRfamily_bootstrap_biclique_rewired_gene_list), function(i) enricher(miRfamily_bootstrap_biclique_rewired_gene_list[[i]], TERM2GENE=cell_markers, minGSSize=1) %>% as.data.frame)

## Hierarchical cluster analysis of cell-specific miRNA-mRNA regulatory network    
    rownames(CsmiR_network_bootstrap_Sim) <- colnames(CsmiR_network_bootstrap_Sim) <- rownames(miRNA_scRNA_norm_filter)

    hclust_res_network_bootstrap <- hclust(as.dist(1-CsmiR_network_bootstrap_Sim), "complete")
 
    dend_network_bootstrap <- as.dendrogram(hclust_res_network_bootstrap)
    dend_network_bootstrap %>% set("branches_k_color", value = c("red", "blue"), k=2) %>% 
    set("labels_col", value = c("red", "blue"), k=2) %>% 
    set("labels_cex", value = 0.7) %>% plot(main="k=2 using network dissimilarity")
    
    dend_network_bootstrap %>% set("branches_k_color", value = c("red", "green", "blue"), k=3) %>% 
    set("labels_col", value = c("red", "green", "blue"), k=3) %>% 
    set("labels_cex", value = 0.7) %>% plot(main="k=3 using network dissimilarity")
    
## Hierarchical cluster analysis of cell-specific hub miRNAs    
    rownames(CsmiR_hub_bootstrap_Sim) <- colnames(CsmiR_hub_bootstrap_Sim) <- rownames(miRNA_scRNA_norm_filter)
 
    hclust_res_hub_bootstrap <- hclust(as.dist(1-CsmiR_hub_bootstrap_Sim), "complete")
 
    dend_hub_bootstrap <- as.dendrogram(hclust_res_hub_bootstrap)
    dend_hub_bootstrap %>% set("branches_k_color", value = c("red", "blue"), k=2) %>% 
    set("labels_col", value = c("red", "blue"), k=2) %>% 
    set("labels_cex", value = 0.7) %>% plot(main="k=2 using hub miRNA dissimilarity")
    
    dend_hub_bootstrap %>% set("branches_k_color", value = c("red", "green", "blue"), k=3) %>% 
    set("labels_col", value = c("red", "green", "blue"), k=3) %>% 
    set("labels_cex", value = 0.7) %>% plot(main="k=3 using hub miRNA dissimilarity")
    
## Similarity network plot in terms of cell-specific miRNA-mRNA regulatory netowork    
    rownames(CsmiR_network_bootstrap_Sim) <- colnames(CsmiR_network_bootstrap_Sim) <- paste("Cell",c(1:19),sep=" ")
    corrplot(CsmiR_network_bootstrap_Sim, method = "pie", type = "upper", diag = FALSE, cl.lim = c(0, 1), tl.cex = 1)

## Similarity network plot in terms of cell-specific hub miRNAs
    rownames(CsmiR_hub_bootstrap_Sim) <- colnames(CsmiR_hub_bootstrap_Sim) <- paste("Cell",c(1:19),sep=" ")
    corrplot(CsmiR_hub_bootstrap_Sim, method = "pie", type = "upper", diag = FALSE, cl.lim = c(0, 1), tl.cex = 1)

## Stem plots
index_net <- data.frame(value = unlist(lapply(seq(CSmiR_network_bootstrap), function(i) nrow(CSmiR_network_bootstrap[[i]]))), id = seq(19))

col1 <- rep("#FF9999", 19)
p1 <- ggplot(index_net, aes(x = id, y = value)) +
    geom_point(aes(color = col1), size = 5) +
    geom_bar(aes(fill = col1), stat = "identity", width = 0.2) +
    #theme_bw(base_family = "Times") +
    xlab("Single-cell ID") +
    ylab("#Predicted cell-specific miRNA-mRNA interactions") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position="none", 
	  panel.border = element_blank(),
          axis.text.x = element_text(face = "bold"),
	  axis.text.y = element_text(face = "bold"),
	  axis.title.x = element_text(face = "bold"),
	  axis.title.y = element_text(face = "bold")) +
    scale_x_continuous(breaks = seq(1, 19, 1)) 

index_validated_net <- data.frame(value = unlist(lapply(seq(CSmiR_network_bootstrap_validated), function(i) nrow(CSmiR_network_bootstrap_validated[[i]])/nrow(CSmiR_network_bootstrap[[i]])*100)), id = seq(19))

col2 <- rep("plum4", 19)
p2 <- ggplot(index_validated_net, aes(x = id, y = value)) +
    geom_point(aes(color = col2), size = 5) +
    geom_bar(aes(fill = col2), stat = "identity", width = 0.2) +
    scale_fill_manual(values=c("plum4"), aesthetics = "fill") +
    scale_colour_manual(values=c("plum4"), aesthetics = "colour") +
    xlab("Single-cell ID") +
    ylab("%Validated cell-specific miRNA-mRNA interactions") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position="none", 
	  panel.border = element_blank(),
          axis.text.x = element_text(face = "bold"),
	  axis.text.y = element_text(face = "bold"),
	  axis.title.x = element_text(face = "bold"),
	  axis.title.y = element_text(face = "bold")) + 
    scale_x_continuous(breaks = seq(1, 19, 1)) 

index_CML_net <- data.frame(value = unlist(lapply(seq(CSmiR_network_bootstrap_CML), function(i) nrow(CSmiR_network_bootstrap_CML[[i]])/nrow(CSmiR_network_bootstrap[[i]])*100)), id = seq(19))

col3 <- rep("blue", 19)
p3 <- ggplot(index_CML_net, aes(x = id, y = value)) +
    geom_point(aes(color = col3), size = 5) +
    geom_bar(aes(fill = col3), stat = "identity", width = 0.2) +
    scale_fill_manual(values=c("blue"), aesthetics = "fill") +
    scale_colour_manual(values=c("blue"), aesthetics = "colour") +
    xlab("Single-cell ID") +
    ylab("%CML-related cell-specific miRNA-mRNA interactions") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position="none", 
	  panel.border = element_blank(),
          axis.text.x = element_text(face = "bold"),
	  axis.text.y = element_text(face = "bold"),
	  axis.title.x = element_text(face = "bold"),
	  axis.title.y = element_text(face = "bold")) +
    scale_x_continuous(breaks = seq(1, 19, 1))

index_CML_hub <- data.frame(value = unlist(lapply(seq(hub_miRNAs_bootstrap_CML), function(i) length(hub_miRNAs_bootstrap_CML[[i]])/length(hub_miRNAs_bootstrap[[i]])*100)), id = seq(19))

col4 <- rep("green", 19)
p4 <- ggplot(index_CML_hub, aes(x = id, y = value)) +
    geom_point(aes(color = col4), size = 5) +
    geom_bar(aes(fill = col4), stat = "identity", width = 0.2) +
    scale_fill_manual(values=c("green"), aesthetics = "fill") +
    scale_colour_manual(values=c("green"), aesthetics = "colour") +
    xlab("Single-cell ID") +
    ylab("%CML-related cell-specific hub miRNAs") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position="none", 
	  panel.border = element_blank(),
          axis.text.x = element_text(face = "bold"),
	  axis.text.y = element_text(face = "bold"),
	  axis.title.x = element_text(face = "bold"),
	  axis.title.y = element_text(face = "bold")) +
    scale_x_continuous(breaks = seq(1, 19, 1))

library(patchwork)
(p1+p2)/(p3+p4) + plot_annotation(tag_levels = 'A')

save.image("CSmiR.RData")

