setwd("../../media/rachel/Windows/Users/All Users/Documents/project_jnu/K-DESeq2/Final_Shot_25_04_2018/Stress/Network/")
library(WGCNA)
library(DESeq2)
library(edgeR)
library(preprocessCore)
library(ComplexHeatmap)
library(flashClust)
all_counts = read.csv("all_htseq-count.csv", row.names = 1, sep = ",")
all_counts = all_counts[1:33977,]
colnames(all_counts) = c("1A","1B","2A","2B","3A","3B","4A","4B","5A","5B","6A","6B")

#tfs_clust = read.csv("DEGenes Network/tfs-in-clust.csv", row.names = 1)
norm_vals = data.frame(normalize.quantiles(log2(cpm(all_counts + 1))))
rownames(norm_vals) = rownames(all_counts)
colnames(norm_vals) = c("1A","1B","2A","2B","3A","3B","4A","4B","5A","5B","6A","6B")

# Genes as nodes
net_nodes = read.table("../../Network/stress/network_nodes_sig_clusters", header = T)
#tpm_vals = read.csv("Arabidopsis_TPM.csv", row.names = 1)
net_data = norm_vals[match(net_nodes$GeneID, rownames(norm_vals)),]
gene_number = 1659

Cor_mat = data.frame(matrix(NA, nrow = gene_number, ncol = gene_number))
for(i in 1:gene_number)
{
  for (j in 1:gene_number)
  {
    Cor_mat[i,j] = cor(as.numeric(net_data[i,]),as.numeric(net_data[j,]), method = "pearson")
  }
}
#Cor_mat = (Cor_mat +1)/2
#Cor_mat = Cor_mat^3

powers = c(c(1:10), seq(from = 12, to=40, by=2))
sft = pickSoftThreshold(t(net_data), powerVector = powers, verbose = 5)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit,signed R^2',
     type='n', main = paste('Scale independence'));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=1,col='red'); abline(h=0.90,col='red')

adj_matrix <- adjacency.fromSimilarity(as.matrix(Cor_mat), power=12, type='signed')
rownames(adj_matrix) = rownames(net_data)
colnames(adj_matrix) = rownames(net_data)

Heatmap(as.data.frame(adj_matrix), cluster_rows = T, cluster_columns = T, show_column_names = F, show_row_names = F, name = "correlation")

gene_tree <- flashClust(as.dist(1 - adj_matrix), method="average")
plot.new()
par(mfrow=c(1,1))
plot(gene_tree)
w = as.dist(1 - adj_matrix)
#module_labels <- cutreeDynamicTree(dendro=gene_tree, minModuleSize=15, deepSplit=F)
module_labels <- cutreeDynamic(dendro = gene_tree, distM = w, deepSplit = F, pamRespectsDendro = FALSE,
                               minClusterSize = 25, method = "tree")
#assign module colours
module_colors = labels2colors(module_labels)

#plot the dendrogram and corresponding colour bars underneath
plotDendroAndColors(gene_tree, module_colors, 'Module colours', dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main='')
#module_colors <- labels2colors(module_labels)
##################################################################

dissTOM =TOMdist(adj_matrix) 
geneTree = flashClust(as.dist(dissTOM),method="average")
# here we define the modules by cutting branches
moduleLabelsManual1=cutreeDynamic(dendro=geneTree,distM=dissTOM,
                                  method="hybrid",deepSplit=2,pamRespectsDendro=F,minClusterSize=30)


# Relabel the manual modules so that their labels
# match those from our previous analysis
moduleLabelsManual2= matchLabels(moduleLabelsManual1,moduleLabelsAutomatic)
# Convert labels to colors for plotting
moduleColorsManual2=labels2colors(moduleLabelsManual2) 
################### Exporting the network  ########################

export_network_to_graphml <- function (adj_mat, filename=NULL, weighted=TRUE,
                                       threshold=0.5, max_edge_ratio=3,
                                       nodeAttr=NULL, nodeAttrDataFrame=NULL,
                                       edgeAttributes=NULL, verbose=FALSE) {
  library('igraph')
  
  # Determine filename to use
  if (is.null(filename)) {
    filename='network.graphml'
  }
  # TODO 2015/04/09
  # Add option to rescale correlations for each module before applying
  # threshold (this is simpler than the previous approach of trying to
  # determine a different threshold for each module)
  #
  # Still, modules with very low correlations should be given somewhat
  # less priority than those with very high correlations.
  
  #module_colors <- unique(nodeAttrDataFrame$color)
  #module_genes <- which(nodeAttrDataFrame$color == color)
  #module_adjmat <- adj_mat[module_genes,]
  #num_genes <- length(module_genes)
  
  # Adjust threshold if needed to limit remaining edges
  max_edges <- max_edge_ratio * nrow(adj_mat)
  
  edge_to_total_ratio <- max_edges / length(adj_mat)
  edge_limit_cutoff <- as.numeric(quantile(abs(adj_mat), 1 - edge_to_total_ratio))
  
  # Also choose a minimum threshold to make sure that at least some edges
  # are left
  min_threshold <- as.numeric(quantile(abs(adj_mat), 0.9999))
  
  threshold <- min(min_threshold, max(threshold, edge_limit_cutoff))
  
  # Remove edges with weights lower than the cutoff
  adj_mat[abs(adj_mat) < threshold] <- 0
  
  # Drop any genes with no edges (TODO: Make optional)
  orphaned <- (colSums(adj_mat) == 0)
  adj_mat <- adj_mat[!orphaned, !orphaned]
  
  # Also remove annotation entries
  if (!is.null(nodeAttr)) {
    nodeAttr <- nodeAttr[!orphaned]
  }
  if (!is.null(nodeAttrDataFrame)) {
    nodeAttrDataFrame <- nodeAttrDataFrame[!orphaned,]
  }
  
  # Keep track of non-positive edges and rescale to range 0,1
  is_zero     <- adj_mat == 0
  is_negative <- adj_mat < 0
  
  adj_mat <- (abs(adj_mat) - threshold) / (max(adj_mat) - threshold)
  adj_mat[is_zero] <- 0
  adj_mat[is_negative] <- -adj_mat[is_negative]
  
  if (verbose) {
    message(sprintf("Outputting matrix with %d nodes and %d edges", 
                    nrow(adj_mat), sum(adj_mat > 0)))
  }
  
  # Create a new graph and add vertices
  # Weighted graph
  if (weighted) {
    g <- graph.adjacency(adj_mat, mode='undirected', weighted=TRUE, diag=FALSE)
  } else {
    adj_mat[adj_mat != 0] <- 1
    g <- graph.adjacency(adj_mat, mode='undirected', diag=FALSE)
  }
  
  # Add single node annotation from vector
  if (!is.null(nodeAttr)) {
    g <- set.vertex.attribute(g, "attr", value=nodeAttr)
  }
  
  # Add node one or more node annotations from a data frame
  if (!is.null(nodeAttrDataFrame)) {
    for (colname in colnames(nodeAttrDataFrame)) {
      g <- set.vertex.attribute(g, colname, value=nodeAttrDataFrame[,colname])
    }
  }
  
  edge_correlation_negative <- c()
  
  # neg_correlations[edge_list]
  edge_list <- get.edgelist(g)
  
  for (i in 1:nrow(edge_list)) {
    from <- edge_list[i, 1]    
    to   <- edge_list[i, 2]    
  }
  
  # Save graph to a file
  write.graph(g, filename, format='graphml')
  
  # return igraph
  return(g)
}
source("https://bioconductor.org/biocLite.R")
biocLite("org.At.tair.db")
source("https://bioconductor.org/biocLite.R")
biocLite("arabidopsis.db0")
library(org.At.tair.db)
keytypes(org.At.tair.db)
gene_ids <- rownames(adj_matrix)
gene_info <- select(org.At.tair.db, keytype='TAIR', keys=gene_ids, 
                    columns=c('GO', 'EVIDENCE', 'ONTOLOGY', 'PATH','GOALL'))
gene_info <- gene_info[!duplicated(gene_info$TAIR),]
# use OrganismDb to retrieve gene annotations
gene_info <- cbind(gene_info, module=module_colors)

#colnames(gene_info) <- c('gene_id', 'description', 'chr', 'strand')

# for now, just grab the description for the first transcript
#gene_info <- gene_info[!duplicated(gene_info$gene_id),]

library(gplots)
# Include RGB versions of module colors for better assignment in Cytoscape
gene_info$color_rgb <- col2hex(gene_info$module)

# first, it's a good idea to check the distribution of edges weights in our
# correlation matrix. This will help us choose a reasonable cutoff for
# exporting the network.
g <- export_network_to_graphml(adj_matrix, filename='../../Network/stress/network-sig_clusters-5-05-2018.graphml',
                               threshold=0.7, nodeAttrDataFrame=gene_info)


MEs = moduleEigengenes(t(net_data), colors = module.colours, excludeGrey = FALSE)$eigengenes

######################################################################################
# https://www.bioconductor.org/packages/devel/bioc/vignettes/CVE/inst/doc/WGCNA_from_TCGA_RNAseq.html
# contruction of correlation matrix
s = abs(bicor(t(tfs_data)))

# picking threshold

powers = c(c(1:10), seq(from = 12, to=40, by=2))
sft = pickSoftThreshold(t(tfs_data), powerVector = powers, verbose = 5)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab='Soft Threshold (power)',ylab='Scale Free Topology Model Fit,signed R^2',
     type='n', main = paste('Scale independence'));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=1,col='red'); abline(h=0.90,col='red')

#calculation of adjacency matrix
beta = 18
a = s^beta

#dissimilarity measure
w = 1-a

# identification of modules

#create gene tree by average linkage hierarchical clustering 
geneTree = flashClust(as.dist(w), method = 'average')

#module identification using dynamic tree cut algorithm
modules = cutreeDynamic(dendro = geneTree, distM = w, deepSplit = 4, pamRespectsDendro = FALSE,
                        minClusterSize = 30, method = "tree")
#assign module colours
module.colours = labels2colors(modules)

#plot the dendrogram and corresponding colour bars underneath
plotDendroAndColors(geneTree, module.colours, 'Module colours', dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main='')
#calculate eigengenes
MEs = moduleEigengenes(t(dg_cipk_mk), colors = module.colours, excludeGrey = FALSE)$eigengenes

#calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);

#cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = 'average');

#plot the result with phytools package
par(mar=c(2,2,2,2))
plot.phylo(as.phylo(METree),type = 'cladogram',show.tip.label = FALSE, main='')
tiplabels(frame = 'circle',col='black', text=rep('',length(unique(modules))), bg = levels(as.factor(module.colours)))

ME_cor =  data.frame(matrix(NA, nrow = 14, ncol = 14))
colnames(ME_cor) = colnames(MEs)
rownames(ME_cor) = colnames(MEs)
for(i in 1:14)
{
  for (j in 1:14)
  {
    ME_cor[i,j] = cor(as.numeric(MEs[,i]),as.numeric(MEs[,j]), method = "pearson")
  }
}
Heatmap(as.data.frame(ME_cor), cluster_rows = T, cluster_columns = T, show_column_names = T, show_row_names = T, name = "correlation")


#https://github.com/iscb-dc-rsg/2016-summer-workshop/blob/master/3B-Hughitt-RNASeq-Coex-Network-Analysis/tutorial/README.md
