cat(
    "
    Importing libraries: data.table\n
    "
)
library(data.table)

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

    for(dis in diseases[[1]]){

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

load_data <- list(
    
    "diseases" = load_diseases,
    "tfs" = load_tfs,
    "genes" = load_genes,
    "binary_networks" = load_binary_networks,
    "triclusters" = load_triclusters,
    "signatures" = load_signatures,
    "splits" = load_splits,
    "fisher_wrt_splits" = load_fisher_wrt_splits

)

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

signatures_and_splits <- list(

    "load_disregnet_tensor" = load_disregnet_tensor,
    "get_signature" = get_signature,
    "digitsum" = digitsum, 
    "complement_code" = complement_code,
    "signature2split" = signature2split

)

#
###--- Fisher test ---###
#

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

fisher_wrt_splits <- list(

    "diseases2code" = diseases2code,
    "perform_fisher" = perform_fisher_wrt_splits

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
