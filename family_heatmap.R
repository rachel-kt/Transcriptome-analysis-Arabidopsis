cluster_matrix_hm = read.csv("tf_heatmap_9.csv", stringsAsFactors = F, row.names = 2)

cluster_matrix_merged = data.frame(cbind(  rowMeans(cluster_matrix_hm[,2:3]), +
                                           rowMeans(cluster_matrix_hm[,4:5]), +
                                           rowMeans(cluster_matrix_hm[,6:7]), +
                                           rowMeans(cluster_matrix_hm[,8:9])), stringsAsFactors = F)
#cluster_matrix_merged$geneid = cluster_matrix_hm$geneid
cluster_matrix_merged$family = cluster_matrix_hm$Family

colnames(cluster_matrix_merged) = c("WTplusK", "WTminusK","MutantplusK", "MutantminusK", "family")
#library(ComplexHeatmap)
matrix_norm = t(scale(t(cluster_matrix_merged[,1:4]), center = T, scale = T))
matrix_norm[is.na(matrix_norm)] <- 0
cn = c("WT", "WT-K", "cipk9", "cipk9-K")
hm = Heatmap(matrix_norm, split = cluster_matrix_merged$family, 
        cluster_rows = T, name = " ", 
        row_names_gp = gpar(fontsize = 6), show_column_dend = F,
        show_row_names = T,
        column_names_gp = gpar(fontsize=12, fontface="italic"),
        bottom_annotation = HeatmapAnnotation(text = anno_text(cn, rot = 45, just = "right", gp=gpar(fontface="italic")),annotation_height = max_text_width(cn)),
        show_row_dend = F,
        show_column_names = F,
        left_annotation = rowAnnotation(pt = anno_empty(border = F, width = unit(1, "mm")))
        )
draw(hm, padding = unit(c(2, 0, 10, 2), "mm"))

        
