###pull out files for COMET analysis
matrix_cometsc <- Seurat::GetAssayData(WT_combined)
write.table(as.matrix(matrix_cometsc), file= "markers.txt", row.names=TRUE, col.names=TRUE, sep = "\t", quote = FALSE)
#in terminal: sed '1s/.*/\t&/' markers.txt > markers2.txt
write.table(WT_combined@active.ident, file='Cell_Cluster.txt', quote=FALSE, sep='\t', col.names = FALSE)
write.table(WT_combined@reductions$umap@cell.embeddings, file='UMAP.txt', quote=FALSE, sep='\t', col.names = FALSE)
###
