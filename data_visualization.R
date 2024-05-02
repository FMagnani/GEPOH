library(ggplot2)
library(data.table)
library(gridExtra)
library(ggpubr)
library(grid)

load("C:\\Users\\fede\\Desktop\\gepoh\\DrugRepurposing\\HierarchicalClustering\\DATA\\data_pack_1.RData")
load("C:\\Users\\fede\\Desktop\\gepoh\\DrugRepurposing\\HierarchicalClustering\\DATA\\data_pack_2.RData")
load("C:\\Users\\fede\\Desktop\\gepoh\\DrugRepurposing\\HierarchicalClustering\\DATA\\data_pack_3.RData")

diseases <- data_pack_2$diseases
genes <- data_pack_2$genes
tfs <- data_pack_2$tfs
tfs_ens <- data_pack_2$tfs_ens

tfs_ens2use <- tfs_ens[tfs_ens %in% genes]
genes_not_tfs <- genes[!(genes %in% tfs_ens2use)]

sj <- data_pack_3$sj
node_outliers <- data_pack_3$node_outliers
mj <- data_pack_3$mj

sij <- data_pack_1$sij
edge_outliers <- data_pack_1$edge_outliers
mij <- data_pack_1$mij

dimnames(sj)[[1]] <- diseases
dimnames(sj)[[2]] <- genes
dimnames(node_outliers)[[1]] <- diseases
dimnames(node_outliers)[[2]] <- genes
names(mj) <- genes

dimnames(sij)[[1]] <- diseases
dimnames(sij)[[2]] <- genes
dimnames(sij)[[3]] <- tfs
dimnames(edge_outliers)[[1]] <- diseases
dimnames(edge_outliers)[[2]] <- genes
dimnames(edge_outliers)[[3]] <- tfs
dimnames(mij)[[1]] <- genes
dimnames(mij)[[2]] <- tfs

#-------------------------------------------------------------------------------

scorelist <- c()
dislist <- c()
genelist <- c()
for(dis in diseases[1:13]){
  scorelist <- c(scorelist, sj[dis, ])
  dislist <- c(dislist, rep(dis, dim(sj)[2]))
  genelist <- c(genelist, genes)
}
sj_long <- data.frame("score"=scorelist, "dis"=dislist, "gene"=genelist)

sj_long$dis <- factor(sj_long$dis, levels=sort(unique(sj_long$dis), decreasing=TRUE))

p_genes <- ggplot(
  sj_long[sj_long$gene %in% genes_not_tfs, ]
) + theme_minimal() + geom_vline(
  xintercept=c(-2,2), color="red", lty="dashed", linewidth=1
) + geom_boxplot(
  aes(score, dis, fill=dis), notch=TRUE, outlier.size=2, show.legend=FALSE
) + labs(
  title="Specificity score of gene expressions",
  y='', x=""
) + theme(
  axis.ticks.y=element_blank(), axis.text.y=element_blank(),
  axis.text.x = element_text(size=13), plot.title = element_text(size=15, hjust=0.5)
)
p_genes

p_tfs <- ggplot(
  sj_long[sj_long$gene %in% tfs_ens2use, ]
) + theme_minimal() + geom_vline(
  xintercept=c(-2,2), color="red", lty="dashed", linewidth=1
) + geom_boxplot(
  aes(score, dis, fill=dis), notch=TRUE, outlier.size=2, show.legend=FALSE
) + labs(
  title="Specificity score of TF expressions",
  y='', x=""
) + theme(
  axis.ticks.y=element_blank(), axis.text.y=element_blank(),
  axis.text.x = element_text(size=13), plot.title = element_text(size=15, hjust=0.5)
)
p_tfs

p_common_label <- ggplot(
  data.frame("dis" = factor(diseases, levels=sort(unique(sj_long$dis), decreasing=FALSE)), "pos" = 1:length(diseases))
) + scale_y_reverse() + geom_text(
  aes(x=0, y=pos, label=dis, color=dis), size=5, fontface="bold", show.legend=FALSE
) + theme_void()
p_common_label

figure <- ggarrange(
  p_genes, p_common_label, p_tfs,
  ncol=3, nrow=1, align = "h", widths=c(1,0.15,1)
)

figure <- annotate_figure(
  figure, 
  top = text_grob(
    "Distribution of specificity scores for genes and TFs expression", 
    color = "black", face = "bold", size = 14
  )
)

figure

#-------------------------------------------------------------------------------

node_multiplicities <- function(select){

  dis_n1 <- list()
  dis_n2 <- list()
  dis_n3 <- list()
  for(idis in diseases){
    n1 <- sum(sj[idis, select]>2 & (mj[select]==1))
    n2 <- sum(sj[idis, select]>2 & (mj[select]==2))
    n3 <- sum(sj[idis, select]>2 & (mj[select]==3))
    dis_n1[[idis]] <- n1
    dis_n2[[idis]] <- n2
    dis_n3[[idis]] <- n3
    message(idis)
  }
  
  mult_by_dis <- rbind(
    data.frame("dis" = diseases, "n" = unlist(dis_n1), "multiplicity" = as.factor(1)),
    data.frame("dis" = diseases, "n" = unlist(dis_n2), "multiplicity" = as.factor(2)),
    data.frame("dis" = diseases, "n" = unlist(dis_n3), "multiplicity" = as.factor(3))
  )
  
  return(mult_by_dis)
  
}

multiplicities_by_disease <- function(outliers_data){
  
  dis_n1 <- list()
  dis_n2 <- list()
  dis_n3 <- list()
  for(idis in diseases){
    n1 <- sum(outliers_data[idis, , ] & (mij==1))
    n2 <- sum(outliers_data[idis, , ] & (mij==2))
    n3 <- sum(outliers_data[idis, , ] & (mij==3))
    dis_n1[[idis]] <- n1
    dis_n2[[idis]] <- n2
    dis_n3[[idis]] <- n3
    message(idis)
  }
  
  mult_by_dis <- rbind(
    data.frame("dis" = diseases, "n" = unlist(dis_n1), "multiplicity" = as.factor(1)),
    data.frame("dis" = diseases, "n" = unlist(dis_n2), "multiplicity" = as.factor(2)),
    data.frame("dis" = diseases, "n" = unlist(dis_n3), "multiplicity" = as.factor(3))
  )
  
  return(mult_by_dis)
  
}

gene_mults <- node_multiplicities(genes_not_tfs)
tfs_mults <- node_multiplicities(tfs_ens2use)
edge_mults <- multiplicities_by_disease(edge_outliers)

edge_mults$dis <- factor(edge_mults$dis, levels=sort(unique(edge_mults$dis), decreasing=TRUE))
tfs_mults$dis <- factor(tfs_mults$dis, levels=sort(unique(tfs_mults$dis), decreasing=TRUE))
gene_mults$dis <- factor(gene_mults$dis, levels=sort(unique(gene_mults$dis), decreasing=TRUE))

p_edges <- ggplot(
  edge_mults
) + geom_bar(
  aes(x=n, y=dis, fill=multiplicity), 
  color="black", stat="identity", 
  show.legend=TRUE, position=position_stack(reverse=TRUE)
) + scale_fill_manual(
  values = c("#0066FF", "#00FFFF", "#00FF66")
   #values = c("#060199", "#0066FF", "#00FFFF")
) + labs(
  title="Disease specificity in regulation",
  y = ''
) + theme_minimal() + theme(
  text=element_text(size=13),
  axis.text.x = element_text(size=13),
  axis.text.y = element_text(size=15),  
  plot.title = element_text(size=15, hjust=0.5, face="bold"),
  legend.position=c(.7,.2)
) + scale_x_continuous(
  expression(paste("Number of specific regulations x ", 10^{4})), 
  breaks=c(0, 250000, 500000, 750000, 1000000), 
  labels=c(0, 25, 50, 75, 100)
)
p_edges

#-------------------------------------------------------------------------------

p_genes <- ggplot(
  gene_mults
) + geom_bar(
  aes(x=n, y=dis, fill=multiplicity), 
  color="black", stat="identity", 
  show.legend=TRUE, position=position_stack(reverse=TRUE)
) + scale_y_discrete(position="right") + scale_x_reverse(
) + scale_fill_manual(
  values = c("#0066FF", "#00FFFF", "#00FF66")
  #values = c("#710000", "#AD1E1E", "#F35656")
) + labs(
  title="Disease specificity in expression of genes",
  x="Number of specific genes",
  y = ''
) + theme_minimal() + theme(
  text=element_text(size=13),
  axis.text.x = element_text(size=13),
  plot.title = element_text(size=15, hjust=0.5, face="bold"),
  axis.ticks.y=element_blank(), axis.text.y=element_blank(),
  legend.position=c(.3,.7)
)
p_genes

p_tfs <- ggplot(
  tfs_mults[tfs_mults$multiplicity==1,]
) + geom_bar(
  aes(x=n, y=dis, fill=multiplicity), 
  color="black", stat="identity", 
  show.legend=TRUE, position=position_stack(reverse=TRUE)
) + scale_fill_manual(
  values = c("#0066FF", "#00FFFF", "#00FF66")
  #values = c("#0B4C01", "#31AD1E", "#57EC40")
) + labs(
  title="Disease specificity in expression of TFs",
  x="Number of specific TFs",
  y = ''
) + theme_minimal() + theme(
  text=element_text(size=13),
  axis.text.x = element_text(size=13),
  plot.title = element_text(size=15, hjust=0.5, face="bold"),
  axis.ticks.y=element_blank(), axis.text.y=element_blank(),
  legend.position=c(.7,.7)
)
p_tfs

#-------------------------------------------------------------------------------

p_common_label <- ggplot(
  data.frame("dis"=factor(diseases, levels=sort(unique(sj_long$dis), decreasing=FALSE)), "pos"=1:length(diseases))
) + scale_y_reverse() + geom_text(
  aes(x=0, y=pos, label=dis), size=6#, fontface="bold"
) + theme_void()
p_common_label

figure <- ggarrange(
  p_genes, p_common_label, p_tfs,
  ncol=3, nrow=1, align = "h", widths=c(1,0.15,1)
  #  common.legend=TRUE, legend="right"
  #  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top")
)
figure

#-------------------------------------------------------------------------------

# Common specific regulations
common_specificities_matrix <- function(){

  m <- matrix(0, 13,13)
  for(i in 1:length(diseases)){
    for(j in i:length(diseases)){
      n_common_edges <-  sum(edge_outliers[i,,]&edge_outliers[j,,])
      m[i,j] <- n_common_edges
      m[j,i] <- n_common_edges
      message(diseases[i], '-', diseases[j])
    }
  }
    
  return(m)
  
}

commonMat <- common_specificities_matrix()
colnames(commonMat) <- diseases
rownames(commonMat) <- diseases

no_diagMat <- commonMat
diag(no_diagMat) <- 0
#image(no_diagMat)

sorting <- sort(no_diagMat[upper.tri(no_diagMat)], decreasing=FALSE)
couples <- c()
for(s in sorting){
  for(i  in 1:13){
    for(j in i:13){
      if(commonMat[i,j]==s){
        couples <- c(couples, paste(diseases[i], '-', diseases[j], sep=''))
        next
      }
    }
  }
}
sorting <- 100*sorting/sum(mij>1)

p <- ggplot(
  data.frame("value"=sorting[sorting>1], "names"=couples[sorting>1])
) + geom_bar(
  aes(x=value, y=factor(names, levels=names)), 
  color="white", fill='black', stat="identity"
) + labs(
  title="Common specific regulations",
  x="Percentage of specific regulations with multiplicity greater than one",
  y = ''
) + theme_minimal() + theme(
  text=element_text(size=13),
  axis.text.x = element_text(size=13, color='black'),
  axis.text.y = element_text(size=12),  
  plot.title = element_text(size=15, hjust=0.5, face="bold")
)
p
