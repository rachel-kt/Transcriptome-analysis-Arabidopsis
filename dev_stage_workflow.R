setwd("../Final_Shot_25_04_2018/Dev_Stage/")
library(DESeq2)
library(edgeR)
library(Biobase)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
ht_seq_counts = read.csv("all_htseq-count.csv", sep = ",", row.names = 1)
col_info = read.csv("condition.csv", row.names = 1)
ht_seq_counts = ht_seq_counts[1:33977,]
colnames(ht_seq_counts) = c("1A","1B","2A","2B","3A","3B","5A","5B")
head(ht_seq_counts)
tail(ht_seq_counts)
ht_seq_counts$cov = apply(ht_seq_counts,1, function(x) sd(x)/mean(x))
ht_seq_counts$threshold <- as.logical((ht_seq_counts$cov > 0.4))
count_mat = ht_seq_counts[which(ht_seq_counts$threshold),1:8]
head(count_mat)
cpm.mat <- log2(cpm(count_mat + 1))
mean.vec <- apply(cpm.mat, 1, mean) # row-wise mean
sdvec <- apply(cpm.mat, 1, sd)
plot(mean.vec, sdvec, pch=".", main="sd vs mean", ylab="sd", xlab="Average logCPM")
dds <- DESeqDataSetFromMatrix(countData = count_mat,
                              colData = col_info, design = ~ condition)
dds$condition <- factor(dds$condition, levels = c("Control", "Control_ck9",
                                                  "WTplusK", "MutantplusK"))
design(dds) <- formula(~ condition)
dds <- DESeq(dds, test = "LRT", reduced = ~1)
res_final <- results(dds)
res_table_final = mcols(dds,use.names=TRUE)
head(res_table_final)
rld <- rlogTransformation(dds, blind=TRUE)
print(plotPCA(rld, intgroup=c("condition")))
p.threshold <- 0.01

res_final$threshold <- as.logical((res_final$padj < p.threshold))
genes.deseq <- row.names(res_final)[which(res_final$threshold)]
genes_for_clustering = data.frame(genes.deseq)
index = match(genes.deseq, rownames(cpm.mat))
cluster_matrix = cpm.mat[index,]
head(cluster_matrix)
head(is.infinite(cluster_matrix))
cluster_matrix[cluster_matrix == "-Inf" ]
cluster_matrix_merged = data.frame(cbind(rowMeans(cluster_matrix[,1:2]), +
                                           rowMeans(cluster_matrix[,3:4]), +
                                           rowMeans(cluster_matrix[,5:6]), +
                                           rowMeans(cluster_matrix[,7:8])))
colnames(cluster_matrix_merged) = c("Control", "Control_ck9", "WTplusK", "MutantplusK")
head(cluster_matrix_merged)
write.csv(cluster_matrix_merged, "Dev_Stage_25-04-2018.csv")
########### Heatmap of Clustered Genes
clust_values = read.csv("clust_values.csv", row.names = 1, sep = "\t")
cluster_genes = read.csv("clust_results.csv", row.names = 1, sep = "\t")
heatmap_data = clust_values[match(rownames(cluster_genes), rownames(clust_values)),]
heatmap_data$Cluster = cluster_genes[match(rownames(heatmap_data),rownames(cluster_genes)),]
f3 = colorRamp2(seq(min(heatmap_data[1:4]), max(heatmap_data[1:4]),length = 3),
                c("blue", "white", "red"), space = "LAB")
colnames(heatmap_data) = c("Col-0","cipk9","Col-0 +K","cipk9 +K", "Cluster")
Heatmap(heatmap_data[1:4], col = f3, cluster_rows = T, show_row_dend = F,
        cluster_columns = F, show_row_names = F, name = "z-score", split = heatmap_data$Cluster, column_order = c(1,2,3,4))
# Heatmap saved as Heatmap_of_clustered_genes.pdf

cpm_clustered_genes = cluster_matrix_merged[match(rownames(heatmap_data),rownames(cluster_matrix_merged)),]
cpm_clustered_genes$Cluster = cluster_genes[match(rownames(cpm_clustered_genes),rownames(cluster_genes)),]
write.csv(cpm_clustered_genes, file = "cpm_clustered_genes.csv")

tfs = read.csv("transcription_factors.csv", row.names = 1)
tfs_hm = heatmap_data[match(paste("gene", rownames(tfs), sep = ":"), rownames(heatmap_data)),]
f4 = colorRamp2(seq(min(tfs_hm[1:4]), max(tfs_hm[1:4]),length = 3),
                c("blue", "white", "red"), space = "LAB")

Heatmap(tfs_hm[1:4], col = f4, cluster_rows = T, show_row_dend = F, row_names_gp = gpar(fontsize=4),
        cluster_columns = F, show_row_names = F, name = "z-score", split = tfs_hm$Cluster, column_order = c(1,2,3,4))
# Saved as Heatmap_of_TFs
############################ DE Genes ###############################
deg_conditions = resultsNames(dds)
comb_condition = t(combn(deg_conditions[-1],2))
# comparisons <- list()                                                                                                                                          
# for (i in 1:nrow(comb_condition)) {                                                                                                                                     
#   comparisons[[i]] <- as.character(comb_condition[i,])                                                                                                      
# } 
#vector to store significant genes
up_genes <- data.frame(matrix(NA, nrow = 10000, ncol = 6))
down_genes <- data.frame(matrix(NA, nrow = 10000, ncol = 6))

for (i in 1:nrow(comb_condition)) {
  # generate string contrast formula, "infLM24 - infLM4"
  contrast <- as.list(comb_condition[i,])
  contrast_result <- results(dds, contrast=contrast)
  contrast_2 <- as.list(deg_conditions[i+1])
  contrast2_result <- results(dds, contrast=contrast_2)
  
  contrast_result$threshold_up <- as.logical((contrast_result$padj < p.threshold) & ((contrast_result$log2FoldChange)>2))
  gene_list <- data.frame(row.names(contrast_result)[which(contrast_result$threshold_up)])
  up_genes[,i] = gene_list[1:10000,]
  
  contrast_result$threshold_dn <- as.logical((contrast_result$padj < p.threshold) & ((contrast_result$log2FoldChange)<(-2)))
  gene_list <- data.frame(row.names(contrast_result)[which(contrast_result$threshold_dn)])
  down_genes[,i] = gene_list[1:10000,]
  
  contrast2_result$threshold_up <- as.logical((contrast2_result$padj < p.threshold) & ((contrast2_result$log2FoldChange)>2))
  gene_list <- data.frame(row.names(contrast2_result)[which(contrast2_result$threshold_up)])
  up_genes[,i+3] = gene_list[1:10000,]
  
  contrast2_result$threshold_dn <- as.logical((contrast2_result$padj < p.threshold) & ((contrast2_result$log2FoldChange)<(-2)))
  gene_list <- data.frame(row.names(contrast2_result)[which(contrast2_result$threshold_dn)])
  down_genes[,i+3] = gene_list[1:10000,]
  
}
colnames(up_genes) = c("cipk_WTpk","cipk_cipkpk", "WTpk_cipkpk", "cipk", "WTpk","cipkpk")
colnames(down_genes) = c("cipk_WTpk","cipk_cipkpk", "WTpk_cipkpk", "cipk", "WTpk","cipkpk")
######## cipk9

write.csv(up_genes, file = "up_regulated_genes_4_fc.csv")
write.csv(down_genes, file = "down_regulated_genes_4_fc.csv")

################# Volcano plot ###################

vol_contrast = list("condition_WTplusK_vs_Control")
vol_con_results <- data.frame(results(dds, contrast=vol_contrast))
vol_plot_data =data.frame(cbind(vol_con_results$log2FoldChange, vol_con_results$padj))
rownames(vol_plot_data) = rownames(vol_con_results)
colnames(vol_plot_data) = c("log2FoldChange", "pvalue")

vol_con_results$threshold <- as.logical((vol_con_results$padj < p.threshold & abs(vol_con_results$log2FoldChange)>1))
vol_sig_gene <- data.frame(row.names(vol_con_results)[which(vol_con_results$threshold)], row.names = 1)
head(vol_con_results)
head(vol_plot_data)

#vol_x = vol_plot_data[match(rownames(vol_sig_gene), rownames(vol_plot_data)),1] 
dat_vol = vol_plot_data[match(rownames(vol_sig_gene), rownames(vol_plot_data)),]
plot(x=vol_plot_data$log2FoldChange, y = -log10(vol_plot_data$pvalue), pch=20)
points(x=vol_plot_data[match(rownames(vol_sig_gene), rownames(vol_plot_data)),1],
     y = -log10(vol_plot_data[match(rownames(vol_sig_gene), rownames(vol_plot_data)),2]), 
     col="pink")

########################### PCA PLOT OF FILTERED GENES

pca_data = count_mat[match(rownames(cluster_genes),rownames(count_mat)),]
dds_pca <- DESeqDataSetFromMatrix(countData = pca_data,
                                  colData = col_info, design = ~ condition)
dds_pca$condition <- factor(dds$condition, levels = c("Control", "Control_ck9",
                                                      "WTplusK", "MutantplusK"))
rld_pca <- rlogTransformation(dds_pca, blind=TRUE)
print(plotPCA(rld_pca, intgroup=c("condition")))

########################## COR PLOT
M_lfc_plusK = res_table_final$condition_MutantplusK_vs_Control
M_lfc_minusK = res_table_final$condition_MutantminusK_vs_Control
WT_lfc_plusK = res_table_final$condition_WTplusK_vs_Control
WT_lfc_minusK = res_table_final$condition_WTminusK_vs_Control
X_cor_plot = M_lfc_plusK - WT_lfc_plusK
Y_cor_plot = M_lfc_minusK - WT_lfc_minusK
cor_plot = data.frame(X_cor_plot, Y_cor_plot)
colnames(cor_plot) = c("Log2(CIPK9 +K/WT +K)", "Log2(CIPK -K/WT -K)")
cor_plot$pc <- predict(prcomp(~X_cor_plot+Y_cor_plot, cor_plot))[,1]
ggplot(cor_plot, aes(X_cor_plot, Y_cor_plot, color = pc)) +
  geom_point(shape = 21, size = 2, show.legend = F) +
  theme_minimal() +
  scale_color_gradient(low = "#0091ff", high = "#f0650e") +
  labs(x = "Log2(CIPK9 +K/WT +K)", y = "Log2(CIPK9 -K/WT -K)") +
  geom_smooth(method="lm", show.legend = F, se = F, color = "blue") +
  geom_abline(slope = 1, show.legend = F, color = "black")
cor.test(cor_plot$`Log2(CIPK9 +K/WT +K)`, cor_plot$`Log2(CIPK -K/WT -K)`,
         method = "pearson", conf.level = 0.95)

###################### Volcano Plot 2

vol_contrast = list("condition_Control_ck9_vs_Control")
vol_con_results <- data.frame(results(dds, contrast=vol_contrast))
vol_plot_data =data.frame(cbind(vol_con_results$log2FoldChange, vol_con_results$padj))
rownames(vol_plot_data) = rownames(vol_con_results)
colnames(vol_plot_data) = c("log2FoldChange", "pvalue")

vol_con_results$threshold <- as.logical((vol_con_results$padj < p.threshold & abs(vol_con_results$log2FoldChange)>1))
vol_sig_gene <- data.frame(row.names(vol_con_results)[which(vol_con_results$threshold)], row.names = 1)
head(vol_con_results)
head(vol_plot_data)

#vol_x = vol_plot_data[match(rownames(vol_sig_gene), rownames(vol_plot_data)),1] 
dat_vol = vol_plot_data[match(rownames(vol_sig_gene), rownames(vol_plot_data)),]
plot(x=vol_plot_data$log2FoldChange, y = -log10(vol_plot_data$pvalue), pch=20)
points(x=vol_plot_data[match(rownames(vol_sig_gene), rownames(vol_plot_data)),1],
       y = -log10(vol_plot_data[match(rownames(vol_sig_gene), rownames(vol_plot_data)),2]), 
       col="pink")
