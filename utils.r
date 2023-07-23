cat("Importing libraries: data.table, igraph\n")
library(data.table)
library(igraph)

#
###--- Loading ---###
#

load_diseases <- function(){

    return( c('ALL','AML','BL','CLL','CML','DLBCL','FL','HL','MCL','MDS','MM','MZLs','PTCL') )

}

load_triclusters <- function(){

    path <- "/home/PERSONALE/daniele.dallolio3/GEPOncoHematology/EXTRA_WORK/Tricluster_20230721_allRun_noOverlap_0.0.csv"

    triclusters_df <- as.data.frame(fread( path ))
    triclusters_df$V1 <- NULL
    triclusters_df$vol <- triclusters_df$N_tfs*triclusters_df$N_genes

    return(triclusters_df)

}

load_tfs <- function(){

    path <- "/home/PERSONALE/federico.magnani9/gepoh/diseaseNetworks/"

    tfs <- as.data.frame(fread(paste(path, "tfs.csv", sep='')))$TF
    return(tfs)

}

load_genes <- function(){

    path <- "/home/PERSONALE/federico.magnani9/gepoh/diseaseNetworks/"

    genes <- as.data.frame(fread(paste(path, "genes.csv", sep='')))$genes
    return(genes)

}

load_signatures <- function(){

    path <- "/home/PERSONALE/federico.magnani9/gepoh/main_results/"

    signatures <- as.data.frame(fread(paste(path, "signatures.csv", sep=''), colClasses="character"))
    rownames(signatures) <- signatures$V1
    signatures$V1 <- NULL

    # Check 
    genes <- load_genes()
    tfs <- load_tfs()
    stopifnot(names(signatures) == genes)
    stopifnot(rownames(signatures) == tfs)

    return(signatures)

}

load_splits <- function(){

    path <- "/home/PERSONALE/federico.magnani9/gepoh/main_results/"

    splits <- as.data.frame(fread(paste(path, "splits.csv", sep=''), colClasses="character"))
    rownames(splits) <- splits$V1
    splits$V1 <- NULL

    # Check 
    genes <- load_genes()
    tfs <- load_tfs()
    stopifnot(names(splits) == genes)
    stopifnot(rownames(splits) == tfs)

    return(splits)

}

load_binary_networks <- function(){

    path <- "/home/PERSONALE/federico.magnani9/gepoh/diseaseNetworks/"
    extension <- "_PANDAregnet_BINNED__motif.csv"

    diseases <- load_diseases()
    genes <- load_genes()
    tfs <- load_tfs()

    binary_networks <- list()

    for(dis in diseases){

        tmp_df <- as.data.frame(fread( paste(path, dis, extension, sep='') ))

        dim(tmp_df)

        names(tmp_df) <- genes
        rownames(tmp_df) <- tfs

        binary_networks[[dis]] <- tmp_df

    }

    return(binary_networks)

}

load_fisher_wrt_splits <- function(){

    path <- "/home/PERSONALE/federico.magnani9/gepoh/main_results/"

    fisher_wrt_splits <- fread(
        paste(path, "fisher_wrt_splits.csv", sep=''), 
        col.names=c("idx","diseases", "clusterIn", "notClusterIn", "clusterOut", "notClusterOut", "pval"), 
        colClasses=(c("character", "character", "numeric", "numeric", "numeric", "numeric", "numeric"))
    )

    return(fisher_wrt_splits)

}

load_network_stats <- function(){

    path <- "/home/PERSONALE/federico.magnani9/gepoh/main_results/"

    network_stats <- fread(
        paste(path, "network_stats.csv", sep=''), 
        col.names=c('idx','diseases','N_diseases','vol','clusterIn','oneEdges','oneEdges_clusterIn','n_connected_nodes','n_connected_genes','n_connected_tfs'), 
        colClasses=(c("numeric", "character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric"))
    )

    return(network_stats)

}


load_data <- list(
    
    "diseases" = load_diseases,
    "tfs" = load_tfs,
    "genes" = load_genes,
    "binary_networks" = load_binary_networks,
    "triclusters" = load_triclusters,
    "signatures" = load_signatures,
    "splits" = load_splits,
    "fisher_wrt_splits" = load_fisher_wrt_splits,
    "network_stats" = load_network_stats

)

#
###--- Utils ---###
#

get_signature <- function(num_code){
    
    signature <- paste(unlist(as.character(num_code)), collapse='')
    return(signature)
    
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

signature2split <- function(code){
    
    nOnes <- digitsum(code)
    if(nOnes<7){    
        code <- complement_code(code)
    }

    return(code)
}

diseases2code <- function(str_cluster){

    diseases <- load_diseases()

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

#
###--- Signatures and Splits ---###
#

load_disregnet_tensor <- function(){

    diseases <- load_diseases()
    path <- "/home/PERSONALE/federico.magnani9/gepoh/diseaseNetworks/"

    # The first disreg network is created, as df
    idis <- 1
    dis <- diseases[idis]
    tmpdf <- as.data.frame(fread(file.path(paste0(path, dis, "_PANDAregnet_BINNED__motif.csv"))))
    rownames(tmpdf) <- tmpdf$V1
    tmpdf$V1 <- NULL

    # Placeholder array (3 dimensional) with all zeros
    disregnets <- array(0, dim = c(length(diseases),dim(tmpdf)[2],dim(tmpdf)[1]))

    # First disease entry, as matrix
    disregnets[idis,,] <- t(as.array(as.matrix(tmpdf)))
    tfs <- rownames(tmpdf)
    genes <- colnames(tmpdf)

    # Fill tensor, as matrix
    for (idis in 2:length(diseases)) {
    dis <- diseases[idis]
    tmpdf <- as.data.frame(fread(paste0(path, dis, "_PANDAregnet_BINNED__motif.csv")))
    rownames(tmpdf) <- tmpdf$V1
    tmpdf$V1 <- NULL
    disregnets[idis,,] <- t(as.array(as.matrix(tmpdf)))
    }

    return(disregnets)

}

signatures_and_splits <- list(

    "load_disregnet_tensor" = load_disregnet_tensor

)

#
###--- Fisher test ---###
#

perform_fisher_wrt_splits <- function(splits, triclusters){

    # Stats to compute for each cluster:
    #   4 entries of contingency matrix
    clusterIn_col <- c() # Links in cluster with same code of the cluster
    notClusterIn_col <- c() # Links in the cluster with different code than the cluster
    clusterOut_col <- c() # Links out the cluster with same code of the cluster
    notClusterOut_col <- c() # Links out the cluster with different code than the cluster
    #   Cluster name & code
    diseases_col <- c()
    #   Fisher's test
    pval_col <- c()

    totLinks <- length(splits)

    n_iters <- nrow(triclusters_df)
    for(i in 1:n_iters){

        if( i%%1000==0 ){
            message("iteration ", i, "/", n_iters, "...")
        }

        row_i <- triclusters_df[i, ]
        cluster <- row_i$diseases
    
        # Make code for the cluster
        code <- diseases2code(cluster)
        # IMPORTANTE
        # I codici dei link sono stati complementati alla loro versione con più uni
        # I codici dei cluster sono invece rimasti univoci
        # Quindi in poche parole per verificare la coerenza del cluster bisogna prima complementare anch'esso
        #    alla sua versione con più uni
        if(digitsum(code)<7){
            code <- complement_code(code)
        }
        
        # Collect all splits of the cluster   
        row_tfs <- unlist(strsplit(row_i$tfs, ','))
        row_genes <- unlist(strsplit(row_i$genes, ','))
        cluster_splits <- unlist(unname(splits[row_tfs, row_genes]))
       
        allIn <- length(cluster_splits)
        allOut <- totLinks - allIn
        # Cluster In 
        clusterIn <- length(cluster_splits[cluster_splits==code])
        clusterIn_col[i] <- clusterIn
        # Not Cluster In
        notClusterIn <- allIn-clusterIn
        notClusterIn_col[i] <- notClusterIn
    
        totCluster <- length(splits[splits==code])
        # Cluster Out
        clusterOut <- totCluster - clusterIn
        clusterOut_col[i] <- clusterOut
        # Not Cluster Out
        notClusterOut <- totLinks - allIn - clusterOut
        notClusterOut_col[i] <- notClusterOut

        # Index
        diseases_col[i] <- cluster

        # Fisher's test p-val
        m <- matrix(
            c(  #cluster    not cluster
                clusterIn, notClusterIn,    # In cluster
                clusterOut, notClusterOut   # Out cluster
            ),
            nrow = 2
        )
        pval <- fisher.test(
            m,
            or = 1, 
            alternative = "two.sided"
        )$p.value
        pval_col[i] <- pval
    
    }

    fisher_test <- data.frame(

        "diseases" = as.character(diseases_col),
        "clusterIn" = as.numeric(clusterIn_col),
        "notClusterIn" = as.numeric(notClusterIn_col),
        "clusterOut" = as.numeric(clusterOut_col),
        "notClusterOut" = as.numeric(notClusterOut_col),
        "pval" = as.numeric(pval_col)
    
    )

    return(fisher_test)

}

aggregate_triclusters <- function(significative_triclusters, significative_fisher_test){

    #--- Aggregate significative clusters with respect to diseases ---#

    # The volume of the aggregating rows is summed
    aggregate_df <- aggregate(significative_triclusters$vol, by=list(diseases=significative_triclusters$diseases), FUN=sum)
    aggregate_df$vol <- aggregate_df$x
    aggregate_df$x <- NULL

    # Recover the number of diseases from the cluster name
    N_diseases <- list()
    for (i in 1:nrow(aggregate_df)){    
        N_diseases <- c(N_diseases, length(strsplit(aggregate_df[i,]$diseases, ',')[[1]]))
    }     
    aggregate_df["N_diseases"] <- unlist(N_diseases)

    clusterIn <- c()
    for(cluster in aggregate_df$diseases){

        ftest_filtered <- significative_fisher_test[significative_fisher_test$diseases==cluster, ]$clusterIn
        clusterIn <- c(clusterIn, sum(ftest_filtered))

    }
    aggregate_df$clusterIn <- as.numeric(clusterIn)
    aggregate_df <- aggregate_df[, c("diseases", "N_diseases", "vol", "clusterIn")]

    return(aggregate_df)

}

compute_network_stats <- function(aggregate_df, significative_triclusters, binary_networks, splits){

    n_oneEdges_col <- list()
    n_oneEdges_clusterIn_col <- list()

    tot_iterations <- nrow(aggregate_df)
    iter <- 1

    for(cluster in aggregate_df$diseases){

        message(iter, "/", tot_iterations)
        iter <- iter+1

        #--- EDGES ---#
        sample_dis <- strsplit(cluster, ',')[[1]][1]
        sample_df <- binary_networks[[sample_dis]]
        rownames(sample_df) <- tfs

        # For identifying the clusterIn
        cluster_code <- diseases2code(cluster)
        if(digitsum(cluster_code)<7){
            cluster_code <- complement_code(cluster_code)
        }

        n_clusterIn <- 0
        n_oneEdges <- 0
        n_oneEdges_clusterIn <- 0
        for (i in 1:nrow(significative_triclusters[significative_triclusters$diseases==cluster,])){

            row = significative_triclusters[significative_triclusters$diseases==cluster,][i,]

            unique_tfs <- unlist(strsplit(row$tfs, ','))
            unique_genes <- unlist(strsplit(row$genes, ','))
            bipartite_adj <- sample_df[unique_tfs, unique_genes] 

            filterSplit <- (splits[unique_tfs,unique_genes]==cluster_code)
            n_clusterIn <- n_clusterIn + length( filterSplit[filterSplit] )

            filterOnes <- (bipartite_adj[unique_tfs, unique_genes]==1)
            n_oneEdges <- n_oneEdges + length( filterOnes[filterOnes] )

            # Filter one-valued edges AND perfect splits
            filterBoth <- filterSplit&filterOnes
            n_oneEdges_clusterIn <- n_oneEdges_clusterIn + length(filterBoth[filterBoth])

        }

        n_oneEdges_col[cluster] <- n_oneEdges
        n_oneEdges_clusterIn_col[cluster] <- n_oneEdges_clusterIn

        # Safety check
        stopifnot( n_clusterIn==aggregate_df[aggregate_df$diseases==cluster,]$clusterIn )

    }

    network_stats <- aggregate_df
    network_stats$oneEdges <- n_oneEdges_col
    network_stats$oneEdges_clusterIn <- n_oneEdges_clusterIn_col

    return(network_stats)

}

compute_network_stats <- function(aggregate_df, significative_triclusters, binary_networks, splits){

    n_oneEdges_col <- list()
    n_oneEdges_perfectSplit_col <- list()
    n_connected_nodes <- list()
    n_connected_tfs <- list()
    n_connected_genes <- list()

    tot_iterations <- nrow(aggregate_df)
    iter <- 1

    for(cluster in aggregate_df$diseases){

        message(iter, "/", tot_iterations)
        iter <- iter+1

        # The graph is undirected
        g <- make_empty_graph(directed=F)

        sample_dis <- strsplit(cluster, ',')[[1]][1]
        sample_df <- binary_networks[[sample_dis]]
        rownames(sample_df) <- tfs

        # For identifying the perfect splits
        cluster_code <- diseases2code(cluster)
        if(digitsum(cluster_code)<7){
            cluster_code <- complement_code(cluster_code)
        }

        n_perfectSplit <- 0
        n_oneEdges <- 0
        n_oneEdges_perfectSplit <- 0

        #--- VERTICES ---#
        unique_genes <- c()
        unique_tfs <- c()

        #--- EDGES ---#
        edges <- c()

        for (i in 1:nrow(significative_triclusters[significative_triclusters$diseases==cluster,])){

            row = significative_triclusters[significative_triclusters$diseases==cluster,][i,]

            row_tfs <- unlist(strsplit(row$tfs, ','))
            row_genes <- unlist(strsplit(row$genes, ','))

            unique_genes <- unique( c(unique_genes, row_genes) )
            unique_tfs <- unique( c(unique_tfs, row_tfs) )

            sub_splits <- splits[row_tfs, row_genes]
            sub_adj <- sample_df[row_tfs, row_genes]

            filterSplit <- (sub_splits==cluster_code)
            n_perfectSplit <- n_perfectSplit + length( filterSplit[filterSplit] )

            filterOnes <- (sub_adj==1)
            n_oneEdges <- n_oneEdges + length( filterOnes[filterOnes] )

            # Filter one-valued edges AND perfect splits
            filterBoth <- filterSplit&filterOnes
            n_oneEdges_perfectSplit <- n_oneEdges_perfectSplit + length(filterBoth[filterBoth])

            filterBoth <- data.frame(filterBoth)
            names(filterBoth) <- row_genes

            sub_edges <- c()
            # Edges of this sub-matrix
            for(colname in row_genes){
        
                tmp_tfs <- row_tfs[ filterBoth[,colname] ] 

                # Trick for making an array in which there are the tfs names alternating with always the same gene name
                # This is because we want an edgelist, each couple an edge: gene,tfs, gene,tfs, gene,tfs, ...
                tmp_edges <- rep(c(colname,1), time=length(tmp_tfs))
                tmp_edges[tmp_edges==1] <- tmp_tfs
    
                sub_edges <- c(sub_edges, tmp_edges)

            }

            edges <- c(edges, sub_edges)

        }

        n_oneEdges_col[cluster] <- n_oneEdges
        n_oneEdges_perfectSplit_col[cluster] <- n_oneEdges_perfectSplit
        # Safety check
        stopifnot( n_perfectSplit==aggregate_df[aggregate_df$diseases==cluster,]$clusterIn )

        #--- ADD VERTICES ---#
        # Red: tfs, Green: genes
        g <- add_vertices(g, length(unique_tfs), color="red", type=TRUE)
        g <- add_vertices(g, length(unique_genes), color="green", type=FALSE)
        # Add name attribute to the vertices
        V(g)[ V(g)$color=="red" ]$name <- unique_tfs
        V(g)[ V(g)$color=="green" ]$name <- unique_genes

        #--- ADD EDGES ---#
        # Eventually add 0-edges and 1-edges with different colors
        g <- add_edges(g, edges, color="black")

        # Let's delete the un-connected nodes. They happen since we are considering only the 1-edges.
        g <- delete.vertices(g, degree(g)==0)

        n_connected_nodes[cluster] <- length(V(g))

        connected_genes <- V(g)[V(g)$color=="red"]$name
        n_connected_genes[cluster] <- length(connected_genes)

        connected_tfs <- V(g)[V(g)$color=="green"]$name
        n_connected_tfs[cluster] <- length(connected_tfs)

    }

    network_stats <- aggregate_df
    network_stats$oneEdges <- as.numeric(n_oneEdges_col)
    network_stats$oneEdges_clusterIn <- as.numeric(n_oneEdges_perfectSplit_col)
    network_stats$n_connected_nodes <- as.numeric(n_connected_nodes)
    network_stats$n_connected_genes <- as.numeric(n_connected_genes)
    network_stats$n_connected_tfs <- as.numeric(n_connected_tfs)

    return(network_stats)

}

fisher_wrt_splits <- list(

    "diseases2code" = diseases2code,
    "perform_fisher" = perform_fisher_wrt_splits,
    "aggregate_triclusters" = aggregate_triclusters,
    "network_stats"  = compute_network_stats

)


cat(
    "
    Utils are loaded! ^_^ \n
    Useful objects:\n
    \tload_data\n
    \tsignatures_and_splits\n
    \tfisher_wrt_splits
    "
)