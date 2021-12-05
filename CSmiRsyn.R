## Function for discovering cell-specific miRNA-miRNA positive correlation
# miR1 and miR2: Gene expression values of two miRNAs in single cells
# cell_index: Index of single cells
# boxsize: Size of neighborhood (0.1 in default)
# interp_betw_point: The number of interpolation points between each cell (5 in default), interp_betw_point = 0 is used to compute
# the normalized statistic of edge gx-gy in large number of cells (more than 100)
# bootstrap_num: The number of bootstrapping for interpolating pseudo-cells 
# Output: res is a list of cell-specific miRNA-miRNA positive correlation
CSmiRsyn_edge_bootstrap <- function(miR1, miR2, cell_index, boxsize = 0.1, bootstrap_betw_point = 5, bootstrap_num = 100) {

    cell_num <- length(miR1)
    set.seed(123)
    bootstrap_sample <- lapply(seq(bootstrap_num), function(i) sample(seq(cell_num), bootstrap_betw_point * (cell_num - 1), replace = TRUE))
    miR1_bootstrap <- lapply(seq(bootstrap_num), function(i) c(miR1, miR1[bootstrap_sample[[i]]]))
    miR2_bootstrap <- lapply(seq(bootstrap_num), function(i) c(miR2, miR2[bootstrap_sample[[i]]]))
    res <- do.call(pmedian, lapply(seq(bootstrap_num), 
	                           function(k) csn_edge(miR1_bootstrap[[k]], 
	                           miR2_bootstrap[[k]], 
	                           boxsize = boxsize)[seq(cell_num)]))   
    
    return(res[cell_index])
}

## Identifying cell-specific miRNA synergistic network
CSmiRsyn_net <- function(miRTarget, cell_index, miRExp, mRExp, minSharedmR = 1, p.value.cutoff = 0.05) {
    
    miRTarget <- as.matrix(miRTarget)    
    miRExpNames <- as.matrix(colnames(miRExp))
    
    miR <- miRTarget[, 1]
    mR <- miRTarget[, 2]

    miRSym <- unique(miR)
    mRSym <- unique(mR)

    m2 <- length(miRSym)

    # Initialize variables
    miRInt <- matrix(NA, m2 * (m2 - 1)/2, 2)
    C <- matrix(NA, m2 * (m2 - 1)/2, 3)

    for (i in seq_len(m2 - 1)) {
        for (j in seq(i + 1, m2)) {

            Interin1 <- miRTarget[which(miRTarget[, 1] %in% miRSym[i]), 2]
            Interin2 <- miRTarget[which(miRTarget[, 1] %in% miRSym[j]), 2]

            M1 <- length(Interin1)
            M2 <- length(Interin2)
            M3 <- length(intersect(Interin1, Interin2))
            M4 <- length(mRSym)
            M5 <- 1 - phyper(M3 - 1, M2, M4 - M2, M1)

            if (M3 >= minSharedmR & M5 < p.value.cutoff) {

                miRInt[(i - 1) * m2 + j - sum(seq_len(i)), 1] <- miRSym[i]
                miRInt[(i - 1) * m2 + j - sum(seq_len(i)), 2] <- miRSym[j]

                miRExpIdx1 <- which(miRExpNames %in% miRSym[i])
                miRExpIdx2 <- which(miRExpNames %in% miRSym[j])
		
                # Calculate cell-specific correlation of each miRNA-miRNA pair		
                M6 <- CSmiRsyn_edge_bootstrap(miRExp[, miRExpIdx1], miRExp[, miRExpIdx2], cell_index)
		M7 <- pnorm(-M6)
                
                C[(i - 1) * m2 + j - sum(seq_len(i)), 1] <- M5
                C[(i - 1) * m2 + j - sum(seq_len(i)), 2] <- M6
		C[(i - 1) * m2 + j - sum(seq_len(i)), 3] <- M7
                
            }
        }
    }

    # Extract miRNA-miRNA pairs
    miRInt <- miRInt[which((C[, 1] < p.value.cutoff & C[, 3] < p.value.cutoff) == "TRUE"), ]

    C <- C[which((C[, 1] < p.value.cutoff & C[, 3] < p.value.cutoff) == "TRUE"), ]

    if (is.vector(C)) {
        res_miRInt <- c(miRInt, C)
        names(res_miRInt) <- c("miRNA_1", "miRNA_2", "p_value of shared mRNAs", "correlation", "p_value of correlation")
            
    } else {
        res_miRInt <- cbind(miRInt, C)
        colnames(res_miRInt) <- c("miRNA_1", "miRNA_2", "p_value of shared mRNAs", "correlation", "p_value of correlation")
    }

    return(res_miRInt)
}

## Function for calculating similarity matrix between two list of module groups
# Module.group1: List object, the first list of module group
# Module.group2: List object, the second list of module group
# Output: Sim is a similarity matrix between two list of module groups
Sim.module.group <- function(Module.group1, Module.group2){

    if(class(Module.group1)!="list" | class(Module.group2)!="list") {
    stop("Please check your input module group! The input module group should be list object! \n")
    }

    m <- length(Module.group1)
    n <- length(Module.group2)
    Sim <- matrix(NA, m, n)
    
    for (i in seq(m)){
        for (j in seq(n)){	    
	    overlap_interin <- length(intersect(Module.group1[[i]], Module.group2[[j]]))
	    Sim[i, j] <- overlap_interin/min(length(Module.group1[[i]]), length(Module.group2[[j]]))
	}    
    }

    if (m < n) {        
	GS <- mean(unlist(lapply(seq(m), function(i) Sim[i, max.col(Sim)[i]])))*m/n
    } else if (m == n) {
        GS <- mean(c(unlist(lapply(seq(m), function(i) Sim[i, max.col(Sim)[i]])), 
	           unlist(lapply(seq(n), function(i) Sim[max.col(t(Sim))[i], i]))))
    } else if (m > n) {
        GS <- mean(unlist(lapply(seq(n), function(i) Sim[max.col(t(Sim))[i], i])))*n/m
    }

    return(GS)
}

## Function for calculating cluster coefficients in random networks
# nodes.num: The number of nodes
# edges.num: The number of edges
# perm: The number of permutation
# directed: Logical value, false or true
# Output: Mean and std of cluster coefficients in random networks
Random_net_clusterCoeff <- function(nodes.num, edges.num, perm = 10000, directed = FALSE) {
    set.seed(123)
    res <- c()
    for (i in seq(perm)) {
    g <- sample_pa(n = nodes.num, m = edges.num, directed = directed)
    g <- delete_edges(g, sample(1:gsize(g), size = gsize(g) - edges.num))
        res[i] <- transitivity(g, type="average")
    }

    return(list(mean(res), sd(res)))
}

