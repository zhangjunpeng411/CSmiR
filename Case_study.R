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
library(vroom)
library(doParallel)

## Load utility functions
source("CSmiR.R")
source("CSmiRsyn.R")

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
    CSmiR_network_bootstrap_null <- CSmiR_net_bootstrap(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, 
                                                   boxsize = 0.1, bootstrap_betw_point = 5, 
						   bootstrap_num = 100, p.value.cutoff = 0.01)
    prior_graph <- make_graph(c(t(TargetScan[, 1:2])), directed = TRUE)
    CSmiR_network_bootstrap_null_graph <- lapply(seq(CSmiR_network_bootstrap_null), function(i) make_graph(c(t(CSmiR_network_bootstrap_null[[i]][, 1:2])), directed = TRUE))
    CSmiR_network_bootstrap <- lapply(seq(CSmiR_network_bootstrap_null), function(i) as_data_frame(CSmiR_network_bootstrap_null_graph[[i]] %s% prior_graph))
  
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
    colnames(CSmiR_network_bootstrap_adjacency_matrix) <- rownames(CSmiR_network_bootstrap_adjacency_matrix) <- rownames(miRNA_scRNA_norm_filter) 
    CSmiR_network_bootstrap_adjacency_matrix_graph <- graph_from_adjacency_matrix(CSmiR_network_bootstrap_adjacency_matrix, mode = "undirected")
    
## Identifying cell-cell crosstalk network in terms of hub miRNA similarity matrix    
    CSmiR_hub_bootstrap_adjacency_matrix <- ifelse(CSmiR_hub_bootstrap_Sim > median(CSmiR_hub_bootstrap_Sim[lower.tri(CSmiR_hub_bootstrap_Sim)]), 1, 0)
    diag(CSmiR_hub_bootstrap_adjacency_matrix) <- 0
    colnames(CSmiR_hub_bootstrap_adjacency_matrix) <- rownames(CSmiR_hub_bootstrap_adjacency_matrix) <- rownames(miRNA_scRNA_norm_filter)
    CSmiR_hub_bootstrap_adjacency_matrix_graph <- graph_from_adjacency_matrix(CSmiR_hub_bootstrap_adjacency_matrix, mode = "undirected")

## Identifying hub cells in terms of network similarity    
    CSmiR_network_bootstrap_cell_degree <- degree(CSmiR_network_bootstrap_adjacency_matrix_graph)
    names(CSmiR_network_bootstrap_cell_degree) <- rownames(miRNA_scRNA_norm_filter)
    CSmiR_network_bootstrap_hub_cells <- names(sort(CSmiR_network_bootstrap_cell_degree[which(CSmiR_network_bootstrap_cell_degree!=0)], decreasing=TRUE))[1:ceiling(0.2*length(which(CSmiR_network_bootstrap_cell_degree!=0)))]

## Identifying hub cells in terms of hub miRNA similarity    
    CSmiR_hub_bootstrap_cell_degree <- degree(CSmiR_hub_bootstrap_adjacency_matrix_graph)
    names(CSmiR_hub_bootstrap_cell_degree) <- rownames(miRNA_scRNA_norm_filter)
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
    
    # Experimentally validated miRNA-mRNA interactions associated with the miR-17/92 family
    CSmiR_network_bootstrap_miRfamily_graph <- lapply(seq(CSmiR_network_bootstrap_miRfamily), function(i) make_graph(c(t(CSmiR_network_bootstrap_miRfamily[[i]][, 1:2])), directed = TRUE))
    CSmiR_network_bootstrap_miRfamily_validated <- lapply(seq(CSmiR_network_bootstrap_miRfamily), function(i) as_data_frame(CSmiR_network_bootstrap_miRfamily_graph[[i]] %s% miRTarget_graph))
   
## Discovering maximal bicliques of conserved and rewired the miR-17/92 family regulation   
    CSmiR_network_miRfamily_bootstrap_biclique_conserved <- biclique_network(list(Sub_miR(Overlap_network_bootstrap_miRfamily)), lleast = 2, rleast = 3)
    CSmiR_network_miRfamily_bootstrap_biclique_rewired <- biclique_network(list(Sub_miR(Overlap_network_bootstrap_rewired_miRfamily)), lleast = 2, rleast = 3)

## Enrichment analysis of maximal bicliques of conserved and rewired the miR-17/92 family regulation   
   miRfamily_bootstrap_biclique_conserved_gene_list <- lapply(seq(CSmiR_network_miRfamily_bootstrap_biclique_conserved[[1]]), function(i) CSmiR_network_miRfamily_bootstrap_biclique_conserved[[1]][[i]]$right)
   miRfamily_bootstrap_biclique_rewired_gene_list <- lapply(seq(CSmiR_network_miRfamily_bootstrap_biclique_rewired[[1]]), function(i) CSmiR_network_miRfamily_bootstrap_biclique_rewired[[1]][[i]]$right)

   # GO, KEGG and Reactome enrichment analysis
   miRfamily_bootstrap_biclique_conserved_FEA <- moduleFEA(miRfamily_bootstrap_biclique_conserved_gene_list, padjustvaluecutoff = 0.05)
   miRfamily_bootstrap_biclique_rewired_FEA <- moduleFEA(miRfamily_bootstrap_biclique_rewired_gene_list, padjustvaluecutoff = 0.05)

   miRfamily_bootstrap_biclique_conserved_miRmR_list <- lapply(seq(CSmiR_network_miRfamily_bootstrap_biclique_conserved[[1]]), function(i) c(CSmiR_network_miRfamily_bootstrap_biclique_conserved[[1]][[i]]$left, CSmiR_network_miRfamily_bootstrap_biclique[[1]][[i]]$right))
   miRfamily_bootstrap_biclique_rewired_miRmR_list <- lapply(seq(CSmiR_network_miRfamily_bootstrap_biclique_rewired[[1]]), function(i) c(CSmiR_network_miRfamily_bootstrap_biclique_rewired[[1]][[i]]$left, CSmiR_network_miRfamily_bootstrap_biclique[[1]][[i]]$right))

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
    rownames(CSmiR_network_bootstrap_Sim) <- colnames(CSmiR_network_bootstrap_Sim) <- rownames(miRNA_scRNA_norm_filter)

    hclust_res_network_bootstrap <- hclust(as.dist(1-CSmiR_network_bootstrap_Sim), "complete")
 
    dend_network_bootstrap <- as.dendrogram(hclust_res_network_bootstrap)
    dend_network_bootstrap %>% set("branches_k_color", value = c("red", "blue"), k=2) %>% 
    set("labels_col", value = c("red", "blue"), k=2) %>% 
    set("labels_cex", value = 0.7) %>% plot(main="k=2 using network dissimilarity")
    
    dend_network_bootstrap %>% set("branches_k_color", value = c("red", "green", "blue"), k=3) %>% 
    set("labels_col", value = c("red", "green", "blue"), k=3) %>% 
    set("labels_cex", value = 0.7) %>% plot(main="k=3 using network dissimilarity")
    
## Hierarchical cluster analysis of cell-specific hub miRNAs    
    rownames(CSmiR_hub_bootstrap_Sim) <- colnames(CSmiR_hub_bootstrap_Sim) <- rownames(miRNA_scRNA_norm_filter)
 
    hclust_res_hub_bootstrap <- hclust(as.dist(1-CSmiR_hub_bootstrap_Sim), "complete")
 
    dend_hub_bootstrap <- as.dendrogram(hclust_res_hub_bootstrap)
    dend_hub_bootstrap %>% set("branches_k_color", value = c("red", "blue"), k=2) %>% 
    set("labels_col", value = c("red", "blue"), k=2) %>% 
    set("labels_cex", value = 0.7) %>% plot(main="k=2 using hub miRNA dissimilarity")
    
    dend_hub_bootstrap %>% set("branches_k_color", value = c("red", "green", "blue"), k=3) %>% 
    set("labels_col", value = c("red", "green", "blue"), k=3) %>% 
    set("labels_cex", value = 0.7) %>% plot(main="k=3 using hub miRNA dissimilarity")
    
## Similarity network plot in terms of cell-specific miRNA-mRNA regulatory netowork    
    rownames(CSmiR_network_bootstrap_Sim) <- colnames(CSmiR_network_bootstrap_Sim) <- paste("Cell",c(1:19),sep=" ")
    corrplot(CSmiR_network_bootstrap_Sim, method = "pie", type = "upper", diag = FALSE, cl.lim = c(0, 1), tl.cex = 1)

## Similarity network plot in terms of cell-specific hub miRNAs
    rownames(CSmiR_hub_bootstrap_Sim) <- colnames(CSmiR_hub_bootstrap_Sim) <- paste("Cell",c(1:19),sep=" ")
    corrplot(CSmiR_hub_bootstrap_Sim, method = "pie", type = "upper", diag = FALSE, cl.lim = c(0, 1), tl.cex = 1)

## Stem plots
index_net <- data.frame(value = unlist(lapply(seq(CSmiR_network_bootstrap), function(i) nrow(CSmiR_network_bootstrap[[i]]))), id = seq(19))

col1 <- rep("#FF9999", 19)
p1 <- ggplot(index_net, aes(x = id, y = value)) +
    geom_point(aes(color = col1), size = 5) +
    geom_bar(aes(fill = col1), stat = "identity", width = 0.2) +
    #theme_bw(base_family = "Times") +
    xlab("Single-cell ID") +
    ylab("#Predicted miRNA-mRNA interactions") +
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
    ylab("%Validated miRNA-mRNA interactions") +
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
    ylab("%CML-related miRNA-mRNA interactions") +
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
    ylab("%CML-related hub miRNAs") +
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

## Identifying cell-specific miRNA synergistic networks
num.cores <- 6
# get number of cores to run
cl <- makeCluster(num.cores)
registerDoParallel(cl)  
        
CSmiRsyn <- foreach(i = seq(CSmiR_network_bootstrap_validated), .packages = c("igraph", "pracma", "WGCNA"), 
                   .export = c("CSmiRsyn_edge_bootstrap", "CSmiRsyn_net", "csn_edge")) %dopar% {
                   CSmiRsyn_net(CSmiR_network_bootstrap_validated[[i]], i, 
                   miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter, 
		   minSharedmR = 1, p.value.cutoff = 0.05)
}     
    
# shut down the workers
stopCluster(cl)
stopImplicitCluster()

# get number of cores to run
cl <- makeCluster(num.cores)
registerDoParallel(cl)  

## Identifying cell-specific miRNA synergistic modules
CSmiRsyn_module <- foreach(i = seq(CSmiR_network_bootstrap_validated), .packages = c("igraph", "pracma", "WGCNA", "miRspongeR"), 
                   .export = c("CSmiRsyn_edge_bootstrap", "CSmiRsyn_net", "csn_edge")) %dopar% {
                   netModule(CSmiRsyn[[i]][, 1:2], modulesize = 3)
}

# shut down the workers
stopCluster(cl)
stopImplicitCluster()

## Identifying cell-specific hub synergistic miRNAs
CSmiRsyn_graph <- lapply(seq(CSmiRsyn), function(i) make_graph(c(t(CSmiRsyn[[i]][, 1:2])), directed = FALSE))
CSmiRsyn_degree <- lapply(seq(CSmiRsyn), function(i) degree(CSmiRsyn_graph[[i]]))
CSmiRsyn_hub_miRNAs <- lapply(seq(CSmiRsyn), function(i) 
                              names(sort(CSmiRsyn_degree[[i]][which(CSmiRsyn_degree[[i]]!=0)], 
			      decreasing = TRUE))[1:ceiling(0.2*length(which(CSmiRsyn_degree[[i]]!=0)))])

## CML-related cell-specific miRNA-miRNA interactions
CML <- as.matrix(read.csv("CML.csv", header = FALSE, sep=","))    
CSmiRsyn_CML <- lapply(seq(CSmiRsyn), function(i) CSmiRsyn[[i]][intersect(which(CSmiRsyn[[i]][, 1] %in% as.matrix(CML)), which(CSmiRsyn[[i]][, 2] %in% as.matrix(CML))), ])

## CML-related cell-specific hub synergistic miRNAs
CSmiRsyn_hub_miRNAs_CML <- lapply(seq(CSmiRsyn_hub_miRNAs), function(i) CSmiRsyn_hub_miRNAs[[i]][which(CSmiRsyn_hub_miRNAs[[i]] %in% as.matrix(CML))])

## Stem plots
index_CSmiRsyn <- data.frame(value = unlist(lapply(seq(CSmiRsyn), function(i) nrow(CSmiRsyn[[i]]))), id = seq(19))
col5 <- rep("#FF9999", 19)
p5 <- ggplot(index_CSmiRsyn, aes(x = id, y = value)) +
    geom_point(aes(color = col5), size = 5) +
    geom_bar(aes(fill = col5), stat = "identity", width = 0.2) +
    #theme_bw(base_family = "Times") +
    xlab("Single-cell ID") +
    ylab("#Interactions") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position="none", 
	  panel.border = element_blank(),
          axis.text.x = element_text(face = "bold"),
	  axis.text.y = element_text(face = "bold"),
	  axis.title.x = element_text(face = "bold"),
	  axis.title.y = element_text(face = "bold")) +
    scale_x_continuous(breaks = seq(1, 19, 1)) 

index_CSmiRsyn_hub <- data.frame(value = unlist(lapply(seq(CSmiRsyn_hub_miRNAs), function(i) length(CSmiRsyn_hub_miRNAs[[i]]))), id = seq(19))
col6 <- rep("#990066", 19)
p6 <- ggplot(index_CSmiRsyn_hub, aes(x = id, y = value)) +
    geom_point(aes(color = col6), size = 5) +
    geom_bar(aes(fill = col6), stat = "identity", width = 0.2) +
    scale_fill_manual(values=c("#990066"), aesthetics = "fill") +
    scale_colour_manual(values=c("#990066"), aesthetics = "colour") +
    xlab("Single-cell ID") +
    ylab("#Hubs") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position="none", 
	  panel.border = element_blank(),
          axis.text.x = element_text(face = "bold"),
	  axis.text.y = element_text(face = "bold"),
	  axis.title.x = element_text(face = "bold"),
	  axis.title.y = element_text(face = "bold")) + 
    scale_x_continuous(breaks = seq(1, 19, 1)) 

index_CSmiRsyn_CML <- data.frame(value = unlist(lapply(seq(CSmiRsyn_CML), function(i) nrow(CSmiRsyn_CML[[i]])/nrow(CSmiRsyn[[i]])*100)), id = seq(19))
col7 <- rep("blue", 19)
p7 <- ggplot(index_CSmiRsyn_CML, aes(x = id, y = value)) +
    geom_point(aes(color = col7), size = 5) +
    geom_bar(aes(fill = col7), stat = "identity", width = 0.2) +
    scale_fill_manual(values=c("blue"), aesthetics = "fill") +
    scale_colour_manual(values=c("blue"), aesthetics = "colour") +
    xlab("Single-cell ID") +
    ylab("%CML-related interactions") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position="none", 
	  panel.border = element_blank(),
          axis.text.x = element_text(face = "bold"),
	  axis.text.y = element_text(face = "bold"),
	  axis.title.x = element_text(face = "bold"),
	  axis.title.y = element_text(face = "bold")) +
    scale_x_continuous(breaks = seq(1, 19, 1))

index_CSmiRsyn_hub_CML <- data.frame(value = unlist(lapply(seq(CSmiRsyn_hub_miRNAs_CML), function(i) length(CSmiRsyn_hub_miRNAs_CML[[i]])/length(CSmiRsyn_hub_miRNAs[[i]])*100)), id = seq(19))
col8 <- rep("green", 19)
p8 <- ggplot(index_CSmiRsyn_hub_CML, aes(x = id, y = value)) +
    geom_point(aes(color = col8), size = 5) +
    geom_bar(aes(fill = col8), stat = "identity", width = 0.2) +
    scale_fill_manual(values=c("green"), aesthetics = "fill") +
    scale_colour_manual(values=c("green"), aesthetics = "colour") +
    xlab("Single-cell ID") +
    ylab("%CML-related hubs") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position="none", 
	  panel.border = element_blank(),
          axis.text.x = element_text(face = "bold"),
	  axis.text.y = element_text(face = "bold"),
	  axis.title.x = element_text(face = "bold"),
	  axis.title.y = element_text(face = "bold")) +
    scale_x_continuous(breaks = seq(1, 19, 1))

index_CSmiRsyn_module <- data.frame(value = unlist(lapply(seq(CSmiRsyn_module), function(i) length(CSmiRsyn_module[[i]]))), id = seq(19))
col9 <- rep("#FF99FF", 19)
p9 <- ggplot(index_CSmiRsyn_module, aes(x = id, y = value)) +
    geom_point(aes(color = col9), size = 5) +
    geom_bar(aes(fill = col9), stat = "identity", width = 0.2) +
    scale_fill_manual(values=c("#FF99FF"), aesthetics = "fill") +
    scale_colour_manual(values=c("#FF99FF"), aesthetics = "colour") +
    xlab("Single-cell ID") +
    ylab("#Modules") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position="none", 
	  panel.border = element_blank(),
          axis.text.x = element_text(face = "bold"),
	  axis.text.y = element_text(face = "bold"),
	  axis.title.x = element_text(face = "bold"),
	  axis.title.y = element_text(face = "bold")) + 
    scale_x_continuous(breaks = seq(1, 19, 1)) 

index_CSmiRsyn_module_CML <- data.frame(value = unlist(lapply(seq(CSmiRsyn_module), function(i) length(which(module_CML_EA(miRNA_scRNA_norm_filter, CML, CSmiRsyn_module[[i]]) < 0.05))/length(CSmiRsyn_module[[i]])*100)), id = seq(19))
col10 <- rep("#696969", 19)
p10 <- ggplot(index_CSmiRsyn_module_CML, aes(x = id, y = value)) +
    geom_point(aes(color = col10), size = 5) +
    geom_bar(aes(fill = col10), stat = "identity", width = 0.2) +
    scale_fill_manual(values=c("#696969"), aesthetics = "fill") +
    scale_colour_manual(values=c("#696969"), aesthetics = "colour") +
    xlab("Single-cell ID") +
    ylab("%CML-related modules") +
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
(p5+p7)/(p6+p8)/(p9+p10) + plot_annotation(tag_levels = 'A')

## Similarity matrix in terms of networks, hubs and modules
CSmiRsyn_net_Sim <- Sim.network(CSmiRsyn, CSmiRsyn, directed = FALSE)
CSmiRsyn_hub_Sim <- Sim.hub(CSmiRsyn_hub_miRNAs, CSmiRsyn_hub_miRNAs)
CSmiRsyn_module_Sim <- matrix(NA, 19, 19)
for (i in seq(19)) {
    for (j in seq(19)) {
        CSmiRsyn_module_Sim[i, j] <- Sim.module.group(CSmiRsyn_module[[i]], CSmiRsyn_module[[j]])
    }
}

col <- colorRampPalette(c("black", "blue", "red"))
rownames(CSmiRsyn_net_Sim) <- colnames(CSmiRsyn_net_Sim) <- paste("Cell",c(1:19),sep=" ")
corrplot(CSmiRsyn_net_Sim, method = "pie", type = "upper", diag = FALSE, cl.lim = c(0, 1), tl.cex = 1)

rownames(CSmiRsyn_hub_Sim) <- colnames(CSmiRsyn_hub_Sim) <- paste("Cell",c(1:19),sep=" ")
corrplot(CSmiRsyn_hub_Sim, method = "pie", type = "upper", diag = FALSE, cl.lim = c(0, 1), tl.cex = 1)

rownames(CSmiRsyn_module_Sim) <- colnames(CSmiRsyn_module_Sim) <- paste("Cell",c(1:19),sep=" ")
corrplot(CSmiRsyn_module_Sim, method = "pie", type = "upper", diag = FALSE, cl.lim = c(0, 1), tl.cex = 1)

## Overlap of cell-specific miRNA-miRNA interactions across cells
CSmiRsyn_input <- lapply(seq(CSmiRsyn), function(i) CSmiRsyn[[i]][, 1:2])
CSmiRsyn_Overlap_network_conserved <- Overlap_net(CSmiRsyn_input, Intersect_num = round(length(CSmiRsyn)*0.9))
CSmiRsyn_Overlap_network_union <- Overlap_net(CSmiRsyn_input, Intersect_num = 1)
CSmiRsyn_Overlap_network_two <- Overlap_net(CSmiRsyn_input, Intersect_num = 2)
CSmiRsyn_Overlap_network_union_graph <- make_graph(c(t(CSmiRsyn_Overlap_network_union[, 1:2])), directed = FALSE)
CSmiRsyn_Overlap_network_two_graph <- make_graph(c(t(CSmiRsyn_Overlap_network_two[, 1:2])), directed = FALSE)
CSmiRsyn_Overlap_network_rewired <- as_data_frame(CSmiRsyn_Overlap_network_union_graph %m% CSmiRsyn_Overlap_network_two_graph)
CSmiRsyn_Overlap_network_conserved_graph <- make_graph(c(t(CSmiRsyn_Overlap_network_conserved[, 1:2])), directed = FALSE)
CSmiRsyn_Overlap_network_rewired_graph <- make_graph(c(t(CSmiRsyn_Overlap_network_rewired[, 1:2])), directed = FALSE)
    
## CML-related conserved and rewired miRNA-miRNA synergistic network
CSmiRsyn_Overlap_network_conserved_CML <- CSmiRsyn_Overlap_network_conserved[intersect(which(CSmiRsyn_Overlap_network_conserved[, 1] %in% as.matrix(CML)), which(CSmiRsyn_Overlap_network_conserved[, 2] %in% as.matrix(CML))), ]
CSmiRsyn_Overlap_network_rewired_CML <- CSmiRsyn_Overlap_network_rewired[intersect(which(CSmiRsyn_Overlap_network_rewired[, 1] %in% as.matrix(CML)), which(CSmiRsyn_Overlap_network_rewired[, 2] %in% as.matrix(CML))), ]

## Overlap of cell-specific hub miRNAs across cells
CSmiRsyn_Overlap_hub_miRNAs_conserved <- Overlap_hub(CSmiRsyn_hub_miRNAs, Intersect_num = round(length(CSmiRsyn_hub_miRNAs)*0.9))
CSmiRsyn_Overlap_hub_miRNAs_union <- Overlap_hub(CSmiRsyn_hub_miRNAs, Intersect_num = 1)
CSmiRsyn_Overlap_hub_miRNAs_two <- Overlap_hub(CSmiRsyn_hub_miRNAs, Intersect_num = 2)
CSmiRsyn_Overlap_hub_miRNAs_rewired <- setdiff(CSmiRsyn_Overlap_hub_miRNAs_union, CSmiRsyn_Overlap_hub_miRNAs_two)

## CML-related conserved and rewired hub miRNAs
CSmiRsyn_Overlap_hub_miRNAs_conserved_CML <- CSmiRsyn_Overlap_hub_miRNAs_conserved[which(CSmiRsyn_Overlap_hub_miRNAs_conserved %in% as.matrix(CML))]    
CSmiRsyn_Overlap_hub_miRNAs_rewired_CML <- CSmiRsyn_Overlap_hub_miRNAs_rewired[which(CSmiRsyn_Overlap_hub_miRNAs_rewired %in% as.matrix(CML))]

## Overlap of cell-specific miRNA synergistic modules across cells
CSmiRsyn_module_collapse <- list()
for(j in seq(19)) {
    interin <- CSmiRsyn_module[[j]]
    CSmiRsyn_module_collapse[[j]] <- unlist(lapply(seq(interin), function(i) paste(sort(interin[[i]]), collapse = ", ")))
}

Overlap_module_conserved <- Overlap_module(CSmiRsyn_module_collapse, Intersect_num = round(length(CSmiRsyn_module_collapse)*0.9))
Overlap_module_union <- Overlap_module(CSmiRsyn_module_collapse, Intersect_num = 1)
Overlap_module_two <- Overlap_module(CSmiRsyn_module_collapse, Intersect_num = 2)
Overlap_module_rewired <- setdiff(Overlap_module_union, Overlap_module_two)

## Hierarchical cluster analysis
library(dendextend)
rownames(CSmiRsyn_net_Sim) <- colnames(CSmiRsyn_net_Sim) <- colnames(CSmiRsyn_module_Sim) <- rownames(miRNA_scRNA_norm_filter)
hclust_res_CSmiRsyn_net <- hclust(as.dist(1-CSmiRsyn_net_Sim), "complete")
dend_network_bootstrap <- as.dendrogram(hclust_res_CSmiRsyn_net)
dend_network_bootstrap %>% set("branches_k_color", value = c("red", "green", "blue"), k=3) %>% 
set("labels_col", value = c("red", "green", "blue"), k=3) %>% 
set("labels_cex", value = 0.7) %>% plot(main="k=3 using network dissimilarity")

rownames(CSmiRsyn_hub_Sim) <- colnames(CSmiRsyn_hub_Sim) <- rownames(miRNA_scRNA_norm_filter)
hclust_res_CSmiRsyn_hub <- hclust(as.dist(1-CSmiRsyn_hub_Sim), "complete")
dend_hub_bootstrap <- as.dendrogram(hclust_res_CSmiRsyn_hub)    
dend_hub_bootstrap %>% set("branches_k_color", value = c("red", "green", "blue"), k=3) %>% 
set("labels_col", value = c("red", "green", "blue"), k=3) %>% 
set("labels_cex", value = 0.7) %>% plot(main="k=3 using hub miRNA dissimilarity")

rownames(CSmiRsyn_module_Sim) <- colnames(CSmiRsyn_module_Sim) <- rownames(miRNA_scRNA_norm_filter)
hclust_res_CSmiRsyn_module <- hclust(as.dist(1-CSmiRsyn_module_Sim), "complete")
dend_module_bootstrap <- as.dendrogram(hclust_res_CSmiRsyn_module)    
dend_module_bootstrap %>% set("branches_k_color", value = c("red", "green", "blue"), k=3) %>% 
set("labels_col", value = c("red", "green", "blue"), k=3) %>% 
set("labels_cex", value = 0.7) %>% plot(main="k=3 using module dissimilarity")

d <- dist(cbind(miRNA_scRNA_norm_filter, mRNA_scRNA_norm_filter))
d <- (d-min(d))/(max(d)-min(d))
fit <- hclust(d, method = "complete")

dend_fit <- as.dendrogram(fit)   
dend_fit %>% set("branches_k_color", value = c("red", "green", "blue"), k=3) %>% 
set("labels_col", value = c("red", "green", "blue"), k=3) %>% 
set("labels_cex", value = 0.7) %>% plot(main="k=3 using single-cell transcriptomics data")

## Identifying cell-cell crosstalk network in terms of miRNA synergistic network similarity matrix    
CSmiRsyn_net_Sim_adjacency_matrix <- ifelse(CSmiRsyn_net_Sim > median(CSmiRsyn_net_Sim[lower.tri(CSmiRsyn_net_Sim)]), 1, 0)
diag(CSmiRsyn_net_Sim_adjacency_matrix) <- 0
rownames(CSmiRsyn_net_Sim_adjacency_matrix) <- rownames(CSmiRsyn_net_Sim)
colnames(CSmiRsyn_net_Sim_adjacency_matrix) <- colnames(CSmiRsyn_net_Sim)
CSmiRsyn_net_Sim_adjacency_matrix_graph <- graph_from_adjacency_matrix(CSmiRsyn_net_Sim_adjacency_matrix, mode = "undirected")

## Identifying cell-cell crosstalk network in terms of hub miRNA similarity matrix
CSmiRsyn_hub_Sim_adjacency_matrix <- ifelse(CSmiRsyn_hub_Sim > median(CSmiRsyn_hub_Sim[lower.tri(CSmiRsyn_hub_Sim)]), 1, 0)
diag(CSmiRsyn_hub_Sim_adjacency_matrix) <- 0
rownames(CSmiRsyn_hub_Sim_adjacency_matrix) <- rownames(CSmiRsyn_hub_Sim)
colnames(CSmiRsyn_hub_Sim_adjacency_matrix) <- colnames(CSmiRsyn_hub_Sim)
CSmiRsyn_hub_Sim_adjacency_matrix_graph <- graph_from_adjacency_matrix(CSmiRsyn_hub_Sim_adjacency_matrix, mode = "undirected")

## Identifying cell-cell crosstalk network in terms of miRNA synergistic similarity matrix
CSmiRsyn_module_Sim_adjacency_matrix <- ifelse(CSmiRsyn_module_Sim > median(CSmiRsyn_module_Sim[lower.tri(CSmiRsyn_module_Sim)]), 1, 0)
diag(CSmiRsyn_module_Sim_adjacency_matrix) <- 0
rownames(CSmiRsyn_module_Sim_adjacency_matrix) <- rownames(CSmiRsyn_module_Sim)
colnames(CSmiRsyn_module_Sim_adjacency_matrix) <- colnames(CSmiRsyn_module_Sim)
CSmiRsyn_module_Sim_adjacency_matrix_graph <- graph_from_adjacency_matrix(CSmiRsyn_module_Sim_adjacency_matrix, mode = "undirected")

## Identifying hub cells in terms of miRNA synergistic network similarity    
CSmiRsyn_net_Sim_cell_degree <- degree(CSmiRsyn_net_Sim_adjacency_matrix_graph)
CSmiRsyn_net_Sim_hub_cells <- names(sort(CSmiRsyn_net_Sim_cell_degree[which(CSmiRsyn_net_Sim_cell_degree!=0)], decreasing=TRUE))[1:ceiling(0.2*length(which(CSmiRsyn_net_Sim_cell_degree!=0)))]

## Identifying hub cells in terms of hub miRNA similarity
CSmiRsyn_hub_Sim_cell_degree <- degree(CSmiRsyn_hub_Sim_adjacency_matrix_graph)
CSmiRsyn_hub_Sim_hub_cells <- names(sort(CSmiRsyn_hub_Sim_cell_degree[which(CSmiRsyn_hub_Sim_cell_degree!=0)], decreasing=TRUE))[1:ceiling(0.2*length(which(CSmiRsyn_hub_Sim_cell_degree!=0)))]

## Identifying hub cells in terms of miRNA synergistic module similarity
CSmiRsyn_module_Sim_cell_degree <- degree(CSmiRsyn_module_Sim_adjacency_matrix_graph)
CSmiRsyn_module_Sim_hub_cells <- names(sort(CSmiRsyn_module_Sim_cell_degree[which(CSmiRsyn_module_Sim_cell_degree!=0)], decreasing=TRUE))[1:ceiling(0.2*length(which(CSmiRsyn_module_Sim_cell_degree!=0)))]

## Identifying cell-cell crosstalk modules in terms of miRNA synergistic network similarity matrix
CSmiRsyn_net_Sim_cell_module <- netModule(CSmiRsyn_net_Sim_adjacency_matrix_graph %>% as_data_frame)

## Identifying cell-cell crosstalk modules in terms of hub miRNA similarity matrix
CSmiRsyn_hub_Sim_cell_module <- netModule(CSmiRsyn_hub_Sim_adjacency_matrix_graph %>% as_data_frame)

## Identifying cell-cell crosstalk modules in terms of miRNA synergistic module similarity matrix
CSmiRsyn_module_Sim_cell_module <- netModule(CSmiRsyn_module_Sim_adjacency_matrix_graph %>% as_data_frame)

save.image("CSmiR.RData")

