library(monocle) #this is monocle version 2, but there is a version 3
library(Seurat) #Seurat V3
library("DDRTree")
library("pheatmap")

#make sure your data has been clustered with Seurat and each cluster has a name
DimPlot(object = GKO_timecourse_ecto, label = TRUE)
#Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(GKO_timecourse_ecto@assays$RNA@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = GKO_timecourse_ecto@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
#Construct monocle cds
HSMM_GKO_timecourse_ecto <- newCellDataSet(data,
                                           phenoData = pd,
                                           featureData = fd,
                                           #lowerDetectionLimit = 0.5,
                                           expressionFamily = uninormal())# since I have already normalized, thresholded and scaled in Seurat
#Run ordering algorithm
#takes the 2000 most variable features between clusters and uses these to order the pseudotime process
#some of these steps will take a while to compute
var_genes <- GKO_timecourse_ecto[["RNA"]]@var.features
ordering_genes <- var_genes
HSMM_GKO_timecourse_ecto<- setOrderingFilter(HSMM_GKO_timecourse_ecto, ordering_genes)
## reduce dimension - do not normalize or include pseudo count. Use monocle scaling
HSMM_GKO_timecourse_ecto <- reduceDimension(HSMM_GKO_timecourse_ecto,norm_method="none", 
                                            reduction_method="DDRTree",
                                            max_components=4,
                                            scaling=TRUE,
                                            verbose=TRUE,
                                            pseudo_expr=0)
## order cells along the trajectory change colors and theta to match your plot if desired
HSMM_GKO_timecourse_ecto <- orderCells(HSMM_GKO_timecourse_ecto)
#plot the results
#can color by sample, Seurat identity, etc.
#plot by identity by changing colo_by="GKO_meso@active.ident"
plot_cell_trajectory(HSMM_GKO_timecourse_ecto,
                     color_by = "time",
                     show_branch_points = TRUE,
                     show_tree = TRUE,
                     cell_size = 1,  label=TRUE)

