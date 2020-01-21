#mergingThresh = 0.15
net = blockwiseModules(t(net_data),corType="pearson",maxBlockSize=300,deepSplit = 4,networkType="signed",power=10,minModuleSize=25, mergeCutHeight=mergingThresh, numericLabels=TRUE,saveTOMs=TRUE,pamRespectsDendro=FALSE,saveTOMFileBase="StressTOM") 
moduleLabelsAutomatic=net$colors
# Convert labels to colors for plotting
moduleColorsAutomatic = data.frame(labels2colors(moduleLabelsAutomatic)) 
MEsAutomatic=net$MEs
rownames(MEsAutomatic) = c("1A","1B","2A","2B","3A","3B","4A","4B","5A","5B","6A","6B")
Heatmap(MEsAutomatic, show_row_names = T, cluster_rows = F, cluster_columns = F)
