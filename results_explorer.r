library(data.table)
library(igraph)

load_results <- function(){
    
    path <- "/home/PERSONALE/federico.magnani9/gepoh/GSEA_results/gsea_significative_pathways.csv"
    df <- as.data.frame(fread(path))
    df$V1 <- NULL
    return(df)

}

diseases <- c('ALL','AML','BL','CLL','CML','DLBCL','FL','HL','MCL','MDS','MM','MZLs','PTCL')
tfs <- as.data.frame(fread("/home/PERSONALE/federico.magnani9/gepoh/diseaseNetworks/tfs.csv"))$TF
genes <- as.data.frame(fread("/home/PERSONALE/federico.magnani9/gepoh/diseaseNetworks/genes.csv"))$genes
cat("Loaded diseases vector with name: diseases\n")
cat("Loaded tfs vector with name: tfs\n")
cat("Loaded genes vector with name: genes\n")

binary_networks <- list()
path_to_binary_nets <- "/home/PERSONALE/federico.magnani9/gepoh/diseaseNetworks/"
extension <- "_PANDAregnet_BINNED__motif.csv"
for(dis in diseases){

    tmp_df <- as.data.frame(fread( paste(path_to_binary_nets, dis, extension, sep='') ))
    names(tmp_df) <- genes
    rownames(tmp_df) <- tfs
    binary_networks[[dis]] <- tmp_df

}
cat("Loaded binary networks in list with name: binary_networks\n")

triclusters_df <- as.data.frame(fread("/home/PERSONALE/daniele.dallolio3/GEPOncoHematology/EXTRA_WORK/Tricluster_20230721_allRun_noOverlap_0.0.csv"))
triclusters_df$V1 <- NULL
cat("Loaded triclusters dataframe with name: triclusters_df\n")

splits <- as.data.frame(fread("/home/PERSONALE/federico.magnani9/gepoh/main_results/splits.csv", colClasses="character"))
rownames(splits) <- splits$V1
splits$V1 <- NULL
cat("Loaded split-codes in dataframe with name: splits\n")

fisher_test <- fread(
    "/home/PERSONALE/federico.magnani9/gepoh/main_results/fisher_wrt_splits.csv", 
    col.names = c("idx","diseases", "clusterIn", "notClusterIn", "clusterOut", "notClusterOut", "pval"), 
    colClasses = c("character", "character", "numeric", "numeric", "numeric", "numeric", "numeric")
)
cat("Loaded fisher test results in dataframe with name: fisher_tests\n")

# Bonferroni corrected Fisher test
# For selecting clusters correlated to the original information
fisherBonferroni <- fisher_test$pval<0.01/nrow(fisher_test)

# Filter by significative clusterIn
tot_entries <- fisher_test[2,]$clusterIn + fisher_test[2,]$clusterOut + fisher_test[2,]$notClusterIn + fisher_test[2,]$notClusterOut
expected_entries <- as.vector(fisher_test$clusterIn+fisher_test$clusterOut)*as.vector(fisher_test$clusterIn+fisher_test$notClusterIn)/tot_entries
significative_clusterIn <- (fisher_test$clusterIn>expected_entries)

# Filter triclusters
filtered_triclusters <- triclusters_df[fisherBonferroni & significative_clusterIn,]
filtered_triclusters$index <- rownames(filtered_triclusters)
filtered_triclusters <- filtered_triclusters[filtered_triclusters$N_diseases!=13,]
cat("Triclusters have been filtered and saved in: filtered_triclusters\n")

diseases2code <- function(str_cluster){

    code = ''
    for(dis in diseases){
        if(dis %in% unlist(strsplit(str_cluster, split=','))){
            code = paste(code, '1', sep='')
        }else{
            code = paste(code, '0', sep='')
        }
    }    

    return(code)

}

digitsum <- function(str_code){
 
    s <- sum(as.numeric(unlist(strsplit(str_code, split=''))))
    return(s)

}

complement_code <- function(str_code){

    num_code <- as.numeric(unlist(strsplit(str_code, split='')))
    complement <- 1-num_code
    str_code <- paste(as.character(complement), collapse='')
    
    return(str_code)

}

edgelist_from_adj <- function(adj, row_genes){
    # From a NAMED matrix (adjacency matrix) it returns the edgelist ready for igraph
    # of the entries with ONE
    # For the entries with zero, just give to the function (1-adj)

    if(sum(adj)==0){
        ones_edgelist <- c()
    }else{

        oneEdges <- which(adj==1, arr.ind=TRUE)
    
        tfsNodes <- rownames(oneEdges)
        geneNodes <- as.numeric(oneEdges[,"col"])
        for(i in unique(as.numeric(oneEdges[,"col"]))){

            geneNodes[geneNodes==i] <- row_genes[i]

        }

        # The igraph edgelist 
        # In odd positions (==1) the GENES
        # In even positions (==2) the TFS
        ones_edgelist <- rep(c(1,2), time=length(tfsNodes))
        ones_edgelist[ones_edgelist==1] <- geneNodes
        ones_edgelist[ones_edgelist==2] <- tfsNodes
    }

    return(ones_edgelist)

}

compute_row_edgelist <- function(row_tfs, row_genes, binary_net, splits, cluster_code){
        
    # If only one column is selected we don't have a (named) matrix but what they call an 'atomic vector'
    # That is, an unnamed unstructured list of numbers
    if(length(row_genes)==1){

        single_gene <- row_genes

        adj_vector <- binary_net[row_tfs, row_genes]

        splits_vector <- splits[row_tfs, row_genes]
        filter_vector <- (splits_vector==cluster_code)
        numericFilterVector <- as.numeric(filter_vector)

        ones_vector <- adj_vector*numericFilterVector
        zeros_vector <- (1-adj_vector)*numericFilterVector

        # In odd positions (==1) the GENES
        # In even positions (==2) the TFS
        row_ones_edgelist <- rep(c(1,2), time=sum(ones_vector))
        row_ones_edgelist[row_ones_edgelist==1] <- single_gene
        row_ones_edgelist[row_ones_edgelist==2] <- row_tfs[ones_vector==1]

        row_zeros_edgelist <- rep(c(1,2), time=sum(zeros_vector))
        row_zeros_edgelist[row_zeros_edgelist==1] <- single_gene
        row_zeros_edgelist[row_zeros_edgelist==2] <- row_tfs[zeros_vector==1]

    }else{

        sub_adj <- binary_net[row_tfs, row_genes]

        sub_splits <- splits[row_tfs, row_genes]    
        filterSplit <- (sub_splits==cluster_code)
        numericFilterSplit <- apply(filterSplit, c(1,2), as.numeric)

        ones_sub_adj <- sub_adj*numericFilterSplit
        zeros_sub_adj <- (1-sub_adj)*numericFilterSplit

        row_ones_edgelist <- edgelist_from_adj(ones_sub_adj, row_genes)
        row_zeros_edgelist <- edgelist_from_adj(zeros_sub_adj, row_genes)

    }

    return(list(
        "ones" = row_ones_edgelist,
        "zeros" = row_zeros_edgelist
    ))

}

build_graph_from_edgelist <- function(ones_edgelist, zeros_edgelist){

    allGenes <- c((ones_edgelist[c(TRUE, FALSE)]), (zeros_edgelist[c(TRUE, FALSE)]))
    allGenes <- allGenes[!is.na(allGenes)]

    allTfs <- c((ones_edgelist[c(FALSE, TRUE)]), (zeros_edgelist[c(FALSE, TRUE)]))
    allTfs <- allTfs[!is.na(allTfs)]

    uniqueGenes <- unique(allGenes)
    uniqueTfs <- unique(allTfs)

    g <- make_empty_graph(directed=F)

    #--- ADD VERTICES ---#
    # Red: tfs, Green: genes
    g <- add_vertices(g, length(uniqueTfs), color="red", type=TRUE)
    g <- add_vertices(g, length(uniqueGenes), color="green", type=FALSE)
    # Add name attribute to the vertices
    V(g)[ V(g)$color=="red" ]$name <- uniqueTfs
    V(g)[ V(g)$color=="green" ]$name <- uniqueGenes

    #--- ADD EDGES ---#
    # 0-edges and 1-edges with different colors
    g <- add_edges(g, ones_edgelist, color="black")
    g <- add_edges(g, zeros_edgelist, color="red")

    g <- delete.vertices(g, degree(g)==0)

    return(g)

}

get_edgelist_from_row_idx <- function(idx){

    row <- filtered_triclusters[filtered_triclusters$index==idx,]

    cluster <- row$diseases
    sample_dis <- strsplit(cluster, ',')[[1]][1]

    binary_net <- binary_networks[[sample_dis]]
    rownames(binary_net) <- tfs

    # For identifying the split-codes
    cluster_code <- diseases2code(cluster)
    if(digitsum(cluster_code)<7){
        cluster_code <- complement_code(cluster_code)
    }
    
    row_tfs <- unlist(strsplit(row$tfs, ','))
    row_genes <- unlist(strsplit(row$genes, ','))

    row_edgelist <- compute_row_edgelist(row_tfs, row_genes, binary_net, splits, cluster_code)

    return(row_edgelist)

}

get_graph_of_row <- function(idx, include_zeros=FALSE){

    row_edgelist <- get_edgelist_from_row_idx(idx)

    if(include_zeros){
        g <- build_graph_from_edgelist(row_edgelist$ones, row_edgelist$zeros)
    }else{
        g <- build_graph_from_edgelist(row_edgelist$ones, c())
    }
    
    return(g)

}

get_graph_of_multiple_rows <- function(idx_list, include_zeros=FALSE){

    ones_edgelist <- c()
    zeros_edgelist <- c()

    for(idx in idx_list){
        row_edgelist <- get_edgelist_from_row_idx(idx)
        ones_edgelist <- c(ones_edgelist, row_edgelist$ones)
        zeros_edgelist <- c(zeros_edgelist, row_edgelist$zeros)
    }

    if(include_zeros){
        g <- build_graph_from_edgelist(ones_edgelist, zeros_edgelist)
    }else{
        g <- build_graph_from_edgelist(ones_edgelist, c())
    }

    return(g)

}
