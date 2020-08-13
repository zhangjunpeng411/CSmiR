############################################################################################
############## Utility functions for exploring cell-specific miRNA regulation ##############
############################################################################################

## The original version of the function is written in Matlab at https://github.com/wys8c764/CSN 
## (Dai H, Li L, Zeng T, Chen L. Cell-specific network constructed by single-cell RNA sequencing data. Nucleic Acids Res. 2019, doi: 10.1093/nar/gkz172.), 
## We reimplement the function in R for single cell RNA sequencing data.  
# gx and gy: Gene expression values of gene x (a vector) and gene y (a vector) in n cells 
# boxsize: Size of neighborhood (0.1 in default)
# Output: res is a vector, the normalized statistic of edge gx-gy in n cells
csn_edge <- function(gx, gy, boxsize = 0.1){

# Define the neighborhood of each plot
    n <- length(gx)
    upper <- zeros(1, n)
    lower <- zeros(1, n)
    a <- zeros(2, n)
    B <- list()

    for (i in seq_len(2)){
        g <- gx*(i==1)+gy*(i==2)
        s1 <- sort(g, index.return = TRUE)[[1]]
        s2 <- sort(g, index.return = TRUE)[[2]]
	n0 <- n - sum(sign(s1))
        h <- round(boxsize/2*sum(sign(s1))+eps(1))
	k <- 1
	while (k <= n){
	    s <- 0
	    while ( (n >= k+s+1) && (s1[k+s+1] == s1[k]) ) {
	    s <- s+1
	    }
	    if (s >= h){
	        upper[s2[k:(k+s)]] <- g[s2[k]]
                lower[s2[k:(k+s)]] <- g[s2[k]]
	    } else {	    	    
                upper[s2[k:(k+s)]] = g[s2[min(n,k+s+h)]]
                lower[s2[k:(k+s)]] = g[s2[max(n0*(n0>h)+1,k-h)]]
	    }	
	    k <- k+s+1
	}

	B[[i]] <- (do.call(cbind, lapply(seq_len(n), function(i) g <= upper[i]))) & 
	          (do.call(cbind, lapply(seq_len(n), function(i) g >= lower[i])))
	a[i,] <- colSums(B[[i]])
}

# Calculate the normalized statistic of edge gx-gy
    res <- (colSums(B[[1]] & B[[2]])*n-a[1,]*a[2,])/sqrt(a[1,]*a[2,]*(n-a[1,])*(n-a[2,])/(n-1)+eps(1))

    return(res)

 }

## Function for computing the average expression values of duplicate genes
# Exp_scRNA: Gene expression values of miRNAs or mRNAs in single cells, rows are cells and columns are miRNAs or mRNAs
# Output: temp is single cell expression data without duplicate genes
Averg_Duplicate <- function(Exp_scRNA){
    
    uniqueNameList <- unique(colnames(Exp_scRNA))
    noOfgenes <- length(uniqueNameList)
    temp <- matrix(0, nrow = nrow(Exp_scRNA), ncol = noOfgenes)
    colnames(temp) <- uniqueNameList
    rownames(temp) <- rownames(Exp_scRNA)
    for(c in 1:noOfgenes){
        GeneList <- which(colnames(Exp_scRNA) == colnames(temp)[c])
    for(r in 1:nrow(temp)) {
        temp[r, c] <- mean(as.numeric(Exp_scRNA[r, GeneList]))  
  }
}
    return(temp)
}

## Function for discovering cell-specific miRNA-mRNA regulatory network
# miR and mR: Gene expression values of miRNAs and mRNAs in single cells, rows are cells and columns are miRNAs or mRNAs
# boxsize: Size of neighborhood (0.1 in default)
# interp_betw_point: The number of interpolation points between each cell (5 in default), interp_betw_point = 0 is used to compute
# the normalized statistic of edge gx-gy in large number of cells (more than 100)
# bootstrap_num: The number of bootstrapping for interpolating pseudo-cells 
# p.value.cutoff: Significance p-value for identifying cell-specific miRNA-mRNA regulatory network
# Output: res_list is a list of cell-specific miRNA-mRNA regulatory network
CSmiR_net_bootstrap <- function(miR, mR, boxsize = 0.1, bootstrap_betw_point = 5, bootstrap_num = 100, p.value.cutoff = 0.01) {

    miRs_num <- ncol(miR)
    mRs_num <- ncol(mR)
    cell_num <- nrow(miR)
    bootstrap_sample <- lapply(seq(bootstrap_num), function(i) sample(seq(cell_num), bootstrap_betw_point * (cell_num - 1), replace = TRUE))
    miR_bootstrap <- lapply(seq(bootstrap_num), function(i) rbind(miR, miR[bootstrap_sample[[i]], ]))
    mR_bootstrap <- lapply(seq(bootstrap_num), function(i) rbind(mR, mR[bootstrap_sample[[i]], ]))
    res <- matrix(NA, nrow = miRs_num*mRs_num, ncol = cell_num + 2)
    for (i in seq(miRs_num)){
        for (j in seq(mRs_num)){
	   res[(i-1)*mRs_num+j, 1] <- colnames(miR)[i]
	   res[(i-1)*mRs_num+j, 2] <- colnames(mR)[j]
	   res[(i-1)*mRs_num+j, 3:(cell_num + 2)] <- do.call(pmedian, lapply(seq(bootstrap_num), 
	                                                     function(k) csn_edge(miR_bootstrap[[k]][, i], 
	                                                     mR_bootstrap[[k]][, j], 
	                                                     boxsize = boxsize)[seq(cell_num)]))
        }
    }
    
    q <- -qnorm(p.value.cutoff)

    res_list <- lapply(seq(cell_num), function(i) res[which(as.numeric(res[, i+2]) > q), seq(2)])

    return(res_list)
}

## Function for identifying the overlap of cell-specific miRNA-mRNA regulatory network
# Netlist: list object, a list of cell-specific miRNA-mRNA regulatory network
# Intersect_num£ºThe least number of different cells intersected for overlap.
# The value of 1 means the union of cell-specific miRNA-mRNA interactions from different cells.
# Output: Overlap_res is the overlap of cell-specific miRNA-mRNA regulatory network in all cells
Overlap_net <- function(Netlist, Intersect_num) {

    if (length(Netlist) >= 2 & length(Netlist) >= Intersect_num) {
        Combcase <- t(combn(length(Netlist), Intersect_num))
        Combnum <- dim(Combcase)[1]
        Overlap_Net <- list()

        for (i in seq_len(Combnum)) {
            Interin <- do.call(rbind, lapply(Combcase[i, ], function(i) Netlist[[i]]))
            Interin_paste <- paste(Interin[, 1], Interin[, 2], sep = "-")
            Interin_table <- table(Interin_paste)
            Interin_names <- names(Interin_table)[which(Interin_table == Intersect_num)]
            Overlap_Net[[i]] <- Interin[which(Interin_paste %in% Interin_names), ]
        }

        Overlap_res <- unique(do.call(rbind, Overlap_Net))
        return(Overlap_res)
    } else {
        stop("Please check your input!\n")
    }

}

## Function for identifying the overlap of cell-specific hub miRNAs
# Netlist: list object, a list of cell-specific hub miRNAs
# Intersect_num£ºThe least number of different cells intersected for overlap.
# The value of 1 means the union of cell-specific hub miRNAs from different cells.
# Output: Overlap_res is the overlap of cell-specific hub miRNAs in all cells
Overlap_hub <- function(hublist, Intersect_num) {

    if (length(hublist) >= 2 & length(hublist) >= Intersect_num) {
        Combcase <- t(combn(length(hublist), Intersect_num))
        Combnum <- dim(Combcase)[1]
        Overlap_Hub <- list()

        for (i in seq_len(Combnum)) {
            Interin <- do.call(rbind, lapply(Combcase[i, ], function(i) hublist[[i]]))
            Interin_table <- table(Interin)
            Interin_names <- names(Interin_table)[which(Interin_table == Intersect_num)]
            Overlap_Hub[[i]] <- Interin_names
        }

        Overlap_res <- unique(unlist(Overlap_Hub))
        return(Overlap_res)
    } else {
        stop("Please check your input!\n")
    }

}

## Function for calculating similarity matrix between two list of networks
# net1: List object, the first list of network
# net2: List object, the second list of network
# directed: Logical value, network directed (TRUE) or undirected (FALSE)
# Output: Sim is a similarity matrix between two list of networks
Sim.network <- function(net1, net2, directed = TRUE){

    if(class(net1)!="list" | class(net2)!="list") {
    stop("Please check your input network! The input network should be list object! \n")
    }

    m <- length(net1)
    n <- length(net2)
    Sim <- matrix(NA, m, n)
    for (i in seq(m)){
        for (j in seq(n)){
	    net1_graph_interin <- make_graph(c(t(net1[[i]][, 1:2])), directed = directed)
            net2_graph_interin <- make_graph(c(t(net2[[j]][, 1:2])), directed = directed)
	    overlap_interin <- nrow(as_data_frame(net1_graph_interin %s% net2_graph_interin))
	    Sim[i, j] <- overlap_interin/min(nrow(net1[[i]]), nrow(net2[[j]]))
	}
    }

    return(Sim)
}

## Function for calculating similarity matrix between two list of hubs
# hub1: List object, the first list of hub
# hub2: List object, the second list of hub
# Output: Sim is a similarity matrix between two list of hubs
Sim.hub <- function(hub1, hub2){

    if(class(hub1)!="list" | class(hub2)!="list") {
    stop("Please check your input hub! The input hub should be list object! \n")
    }

    m <- length(hub1)
    n <- length(hub2)
    Sim <- matrix(NA, m, n)
    for (i in seq(m)){
        for (j in seq(n)){	    
	    overlap_interin <- length(intersect(hub1[[i]], hub2[[j]]))
	    Sim[i, j] <- overlap_interin/min(length(hub1[[i]]), length(hub2[[j]]))
	}
    }

    return(Sim)
}

## Function for identifying maximal bicliques from miRNA-mRNA regulatory network
# net: List object, a list of miRNA-mRNA regulatory network
# lleast: Least number of left partite
# rleast: Least number of right partite
# Output: biclique_res is a list of maximal bicliques
biclique_network <- function(net, lleast = 2, rleast = 3){
    
    biclique_res <- list()
    for (i in seq(net)){
        write.table(net[[i]], file = "Interin.csv", 
                    append = FALSE, quote = FALSE, 
		    sep = "\t", eol = "\n", na = "NA", 
	            dec = ".", row.names = FALSE, col.names = FALSE)
        bi.format("Interin.csv")
        # Compute the bicliques
        biclique_res[[i]] <- bi.clique("Interin.csv", left_least = lleast, right_least = rleast)
    }
    return(biclique_res)
}

## CML enrichment analysis using hypergeometric distribution test
# miRExp and mRExp: Gene expression values of miRNAs and mRNAs in single cells, rows are cells and columns are miRNAs or mRNAs 
# CMLgenes: CML-related genes (miRNAs and mRNAs)
# Modulelist: List object, a list of miRNA-mRNA bicliques or networks
# Output: A list of significance p-values enriched in CML
module_CML_EA <- function(miRExp, mRExp, CMLgenes, Modulelist) {

    ExpData <- cbind(miRExp, mRExp)      

    B <- ncol(ExpData)
    N <- length(intersect(colnames(ExpData), as.matrix(CMLgenes)))
    M <- unlist(lapply(seq_along(Modulelist), function(i) length(Modulelist[[i]])))
    x <- unlist(lapply(seq_along(Modulelist), function(i) length(intersect(Modulelist[[i]], as.matrix(CMLgenes)))))    
    p.value <- 1 - phyper(x - 1, N, B - N, M)
    
    names(p.value) <- names(Modulelist)
    return(p.value)
}

## Function for replacing miR-17/92 family names
# Int: miRNA-target interactions
# Output: unique miRNA-target interactions after replacing miRNA names
Sub_miR <- function(Int){
    tmp <- Int[, 1]
    tmp1 <- gsub("hsa-miR-17-5p", "miR-17", tmp)
    tmp2 <- gsub("hsa-miR-17-3p", "miR-17", tmp1)
    tmp3 <- gsub("hsa-miR-18a-5p", "miR-18a", tmp2)
    tmp4 <- gsub("hsa-miR-18a-3p", "miR-18a", tmp3)
    tmp5 <- gsub("hsa-miR-19a-5p", "miR-19a", tmp4)
    tmp6 <- gsub("hsa-miR-19a-3p", "miR-19a", tmp5)
    tmp7 <- gsub("hsa-miR-19b-1-5p", "miR-19b-1", tmp6)
    tmp8 <- gsub("hsa-miR-19b-3p", "miR-19b-1", tmp7)
    tmp9 <- gsub("hsa-miR-20a-5p", "miR-20a", tmp8)
    tmp10 <- gsub("hsa-miR-20a-3p", "miR-20a", tmp9)
    tmp11 <- gsub("hsa-miR-92a-1-5p", "miR-92a-1", tmp10)
    tmp12 <- gsub("hsa-miR-92a-3p", "miR-92a-1", tmp11)

    res <- unique(cbind(tmp12, Int[, 2]))

    colnames(res) <- c("miR", "target")
    return(res)
}