#setwd("/media/rachel/Part2/OneDrive - south.du.ac.in/project_jnu/K-DESeq2/Final_Shot_25_04_2018/Stress/2020/")
setwd('E:/Rachel/Rachel/')
#data = read.csv('counts_cipk.csv', row.names = 1)
#data$sd = rowSds(as.matrix(data))
#data = data[data$sd!=0,]
#net_dat_degs = read.csv("E:\\Rachel\\Rachel\\counts.csv", row.names = 1, stringsAsFactors = F)
net_dat_degs = read.csv('cv0_4.csv', row.names = 1)
net_dat_degs = net_dat_degs[,1:8]
net_dat_degs = t(scale(t(net_dat_degs), center = T, scale = T))
#net_dat_degs = cluster_matrix
datExpr = t(net_dat_degs)
# is_na <- is.na(datExpr)
# datExpr[is_na] <- 0

#net_dat_degs = read.csv("../all_htseq-count.csv", row.names = 1)
#net_dat_degs = net_dat_degs[1:33977,]
library(WGCNA)
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=25, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, corOptions = list(use = 'pairwise.complete.obs'))
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,1));
cex1 = 1.5;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
#write.csv(t(datExpr), "dattest.csv", na = "NA")
POW = 12
#write.csv(net_dat_degs, file = "data_network_pids_10_12_2019.csv")
# here we define the adjacency matrix using soft thresholding with beta=6
sim_mat = cor(datExpr, method = 'pearson')
ADJ1 = adjacency.fromSimilarity(sim_mat, power=POW, type='unsigned')
#ADJ1=abs(cor(datExpr,use="p"))^POW
# When you have relatively few genes (<5000) use the following code
# k=as.vector(apply(ADJ1,2,sum, na.rm=T))
# When you have a lot of genes use the following code
k=softConnectivity(datExpr = datExpr,power=POW)
# Plot a histogram of k and a scale free topology plot
sizeGrWindow(10,5)
par(mfrow=c(1,2))
hist(k)
scaleFreePlot(k, main="Check scale free topology\n")

# Turn adjacency into a measure of dissimilarity
#dissADJ=1-ADJ1
dissTOM=TOMdist(ADJ1)
#dissTOM=TOMdist(abs(cor(datExpr,use="p"))^POW)
collectGarbage()
hierTOM = hclust(as.dist(dissTOM),method="average");
#hierTOM = hclust(as.dist(TOMdist(abs(cor(datExpr,use="p"))^POW)),method="average");

# The reader should vary the height cut-off parameter h1
# (related to the y-axis of dendrogram) in the following
colorStaticTOM = as.character(cutreeStaticColor(hierTOM, cutHeight=.99, minSize=200))
colorDynamicTOM = labels2colors(cutreeDynamic(hierTOM,method="tree", minClusterSize = 200, deepSplit = 1))
colorDynamicHybridTOM = labels2colors(cutreeDynamic(hierTOM, distM= dissTOM , cutHeight = 0.998,
                                                   deepSplit=1, pamRespectsDendro = FALSE))
# Now we plot the results
 sizeGrWindow(10,5)
 plotDendroAndColors(hierTOM,
                     colors=data.frame(colorStaticTOM,
                                       colorDynamicTOM, colorDynamicHybridTOM),
                     dendroLabels = FALSE, marAll = c(1, 8, 3, 1), cex1 = 1.5, 
                     main = "Gene dendrogram and module colors, TOM dissimilarity")

# sizeGrWindow(10,5)
# plotDendroAndColors(hierTOM,
#                     colors=data.frame(colorDynamicTOM),
#                     dendroLabels = FALSE, marAll = c(1, 8, 3, 1),
#                     main = "Gene dendrogram and module colors, TOM dissimilarity")


#dynamicMods = cutreeDynamic(hierTOM,method="tree", minClusterSize = 20, deepSplit = 4)
dynamicMods = colorDynamicTOM
table(dynamicMods)
mods<- table(dynamicMods)
write.csv(mods, "modules.csv")
dynamicColors = labels2colors(dynamicMods)
moduleColors = dynamicColors
table(dynamicColors)

MEList = moduleEigengenes(datExpr, colors = dynamicColors)
ME_unmerged = MEList
ME_unmerged = MEList$eigengenes

#library(ComplexHeatmap)
#Heatmap(abs(cor(datExpr,use="p"))^POW, show_row_names = F, show_column_names = F)
######################################

MEDiss = 1-cor(ME_unmerged);
METree = hclust(as.dist(MEDiss), method = "average");
MEDissThres = 0.25
par(mfrow=c(1,1))
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors;
mergedMEs = merge$newMEs;
#moduleColors = mergedColors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(mergedColors, colorOrder)-1;
MEs = mergedMEs;


sizeGrWindow(10,5)
plotDendroAndColors(hierTOM,
                    colors=(cbind(moduleColors,mergedColors)),
                    dendroLabels = FALSE, marAll = c(1, 8, 3, 1),
                    main = "Gene dendrogram and merged module colors, TOM dissimilarity")

moduleColors_copy = moduleColors
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
moduleColors = mergedColors
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

trait_table = data.frame(rownames(datExpr))
# trait_table$normal = c(rep(1,12), rep(0,57))
# trait_table$cancer = c(rep(0,12), rep(1,57))
trait_table$wtpk = c(rep(1,2), rep(0,6))
trait_table$wtmk = c(rep(0,2), rep(1,2), rep(0,4))
trait_table$ck9pk = c(rep(0,4), rep(1,2), rep(0,2))
trait_table$ck9mk = c(rep(0,6), rep(1,2))
rownames(trait_table) = trait_table$rownames.datExpr.
trait_table = trait_table[,-1]
datTraits = trait_table
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
sizeGrWindow(10,6)
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 12, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               #xLabels = names(datTraits),
               xLabels = c("WT +K", "WT -K", "CIPK +K", "CIPK -K"),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1.2,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
#cancerstatus = as.data.frame(datTraits$cancer)
cipkminus = as.data.frame(datTraits$ck9mk)
stress = as.data.frame(datTraits)
#names(cancerstatus) = "Cancerstatus"
names(cipkminus) = "cipkmk"
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, cipkminus, use = "p"));
geneTraitSignificance_2 = as.data.frame(cor(datExpr, stress, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
GSPvalue_2 = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance_2), nSamples));
names(geneTraitSignificance) = paste("GS.", names(cipkminus), sep="");
names(GSPvalue) = paste("p.GS.", names(cipkminus), sep="");

###
names(geneTraitSignificance_2) = paste("GS.", names(stress), sep="");
names(GSPvalue_2) = paste("p.GS.", names(stress), sep="");
genes =data.frame(rownames(net_dat_degs), stringsAsFactors = F)
module = "green"
column = match(module, modNames);
moduleGenes = moduleColors==module
green_genes = genes[moduleGenes,]
column = match("royalblue", modNames);
moduleGenes = moduleColors=="royalblue"
royalblue_genes = genes[moduleGenes,]
column = match("grey60", modNames);
moduleGenes = moduleColors=="grey60"
grey60_genes = genes[moduleGenes,]
# column = match("lightcyan", modNames);
# moduleGenes = moduleColors=="lightcyan"
# lightcyan_genes = genes[moduleGenes,]
# column = match("cyan", modNames);
# moduleGenes = moduleColors=="cyan"
# cyan_genes = genes[moduleGenes,]
# column = match("grey", modNames);
# moduleGenes = moduleColors=="grey"
# grey_genes = genes[moduleGenes,]


write.csv(cbind(green_genes,royalblue_genes,grey60_genes), "./network_genes.csv")
write.csv(green_genes, "./green_genes.csv")
write.csv(royalblue_genes, "./rb_genes.csv")
write.csv(grey60_genes, "./grey_genes.csv")
write.csv(GSPvalue_2, "./gene_Sig.csv")
write.csv(geneModuleMembership,"./modmem.csv")
write.csv(geneTraitSignificance_2,"./genetraitsign.csv")

###################################3

# module membership vs gene significance p value
library(ggplot2)
scatter_data = data.frame(cbind(abs(geneModuleMembership$MMgrey60),abs(geneTraitSignificance$`GS.datTraits$ck9mk`)))
#scatter_data = data.frame(cbind((geneModuleMembership$MMcyan),(geneTraitSignificance$GS.Cancerstatus)))
#cyan_g = read.csv("cyan genes", header = F)
x = geneModuleMembership[match(grey60_genes, rownames(geneModuleMembership)),]
y = geneTraitSignificance[match(grey60_genes, rownames(geneTraitSignificance)),]
scatter_genes = data.frame(cbind(abs(x$MMgrey60),abs(y)))


ggplot(scatter_data, aes(x=X1, y=X2)) +
  geom_point(size=2, shape = 16, alpha = 1/3) +
  xlab("Module membership for Grey60 module") +
  ylab("Gene significance for cipk -k") +
  geom_point(mapping = aes(X1, X2) ,data = scatter_genes, colour = "orange", shape = 16)

#####################################################################################
################### Exporting the network  ########################
is_na <- is.na(ADJ1)
ADJ1[is_na] <- 0
export_network_to_graphml <- function (adj_mat, filename=NULL, weighted=TRUE,
                                       threshold=0.9, max_edge_ratio=3,
                                       nodeAttr=NULL, nodeAttrDataFrame=NULL,
                                       edgeAttributes=NULL, verbose=TRUE) {
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
  adj_mat[abs(adj_mat) < 0.9] <- 0
  
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
  
  #  adj_mat <- (abs(adj_mat) - threshold) / (max(adj_mat) - threshold)
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
gene_ids <- rownames(ADJ1)
gene_info <- data.frame(cbind(gene_ids, module=mergedColors))
library(gplots)
# Include RGB versions of module colors for better assignment in Cytoscape
gene_info$color_rgb <- col2hex(gene_info$module)

# first, it's a good idea to check the distribution of edges weights in our
# correlation matrix. This will help us choose a reasonable cutoff for
# exporting the network.
g <- export_network_to_graphml(ADJ1, filename='./network-March-2020_unsigned.graphml',
                               threshold=0.9, nodeAttrDataFrame=gene_info)

edge_list_cancer = data.frame(get.edgelist(g))
write.csv(edge_list_cancer, "./edge_list_Jan_2020.csv")
edge_list_cancer$Correlation = "NULL"
edge_list_cancer$p_n = "NULL"
i=1
Cor_mat_orig = cor(datExpr)
for (i in 1:512) {
  edge_list_cancer[i,3] = Cor_mat_orig[match(edge_list_cancer[i,1],rownames(Cor_mat_orig)), match(edge_list_cancer[i,2],colnames(Cor_mat_orig))]  
  if (edge_list_cancer[i,3] < 0) {
    edge_list_cancer[i,4] = -1  
  }   
  else{edge_list_cancer[i,4] = 1}
}
write.csv(edge_list_cancer, "./network_edges_nov.csv")
write.csv(row.names(net_dat_degs), "./network_vertices_nov.csv")

write.csv(cbind(rownames(net_dat_degs), colorDynamicTOM), "unmerged modules.csv")

write.csv(res_table_final, "table_final.csv")
