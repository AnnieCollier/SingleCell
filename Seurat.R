library(Seurat)
library(stringr)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
library(ggrepel)
library(ggplot2)
library(dplyr)
library(sctransform)

#######
#Seurat analysis for 1 sample
#######

#load your file from CellRanger
#first specify the folder which contains the barcodes, features, and matrix files
input_files="~/Box/Ann_Collier/Single Cell/scRNA/D50_WT1_filtered_feature_bc_matrix/"
D50_WT1.data <- Read10X(data.dir = input_files)
D50_WT1.data <- CreateSeuratObject(counts = D50_WT1.data, min.features = 100, min.cells = 5, project = "D50_WT1")
D50_WT1.data
#add mitochondrial read percentages to the metadata
D50_GKO1.data[["percent.mt"]] <- PercentageFeatureSet(D50_GKO1.data, pattern = "^MT-")
VlnPlot(D50_GKO1.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#remove low quality cells
D50_GKO1.data <- subset(D50_GKO1.data, subset = nFeature_RNA > 1500 & nFeature_RNA < 8000 & percent.mt < 10)
D50_GKO1.data <- NormalizeData(D50_GKO1.data)
D50_GKO1.data <- FindVariableFeatures(D50_GKO1.data, selection.method = "vst")
#cell cycle scoring is optional, but can help remove clusters that are called based on cycling cells 
s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes
D50_GKO1.data <- CellCycleScoring(D50_GKO1.data, g2m.features=g2m_genes, s.features=s_genes)
#scale data and regress out cell cycle effects. If you are not scoring cell cycle just remove the vars.to.regress section
D50_GKO1.data <- ScaleData(D50_GKO1.data, vars.to.regress = c("S.Score", "G2M.Score"))
#run and plot PCA and elbow plot
#wherever the elbow plot bends is a good place to start for dimenstions in the UMAP section
D50_GKO1.data <- RunPCA(D50_GKO1.data, features = VariableFeatures(D50_WT1.data), ndims.print = 1:10, nfeatures.print = 10)
DimPlot(D50_GKO1.data, reduction = "pca")
ElbowPlot(D50_GKO1.data, ndims = 30)
#run clustering and make UMAPs
#here you can change the number of dimensions and the resolution. May need to play around with these until clusters make sense. I usually use 10 dimensions.
D50_GKO1.data <- RunUMAP(D50_GKO1.data, reduction = "pca", dims = 1:10)
D50_GKO1.data <- FindNeighbors(D50_GKO1.data, reduction = "pca", dims = 1:10)
D50_GKO1.data <- FindClusters(D50_GKO1.data, resolution = 0.15)
#plot the UMAP or PCA
DimPlot(D50_GKO1.data, reduction = "umap")
DimPlot(D50_GKO1.data, reduction = "pca")

#remove or subset clusters
D50_GKO1.data_subset <- subset(D50_GKO1.data, idents = c(0, 1, 2))
#rename clusters
D50.new.cluster.ids <- c("Epidermis", "Dermis", "Cartilage")
names(D50.new.cluster.ids) <- levels(D50_GKO1.data)
D50_GKO1.data <- RenameIdents(D50_GKO1.data, D50.new.cluster.ids)
#quantify the proportion of cells in each cluster
table(D50_integ$sample)
table(Idents(D50_integ), D50_integ$sample)
prop.table(table(Idents(D50_integ), D50_integ$celltype), 2)


#######
#Seurat analysis for multiple samples merged together
#######

#merge samples after input
D50_D7_WT <- merge(x = WT1.data,y = c(WT2.data, D50_WT1.data, D50_WT2.data))
D50_D7_WT[["percent.mt"]] <- PercentageFeatureSet(D50_D7_WT, pattern = "^MT-")
#plot your samples to make sure everything looks correct
VlnPlot(D50_D7_WT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
D50_D7_WT <- subset(D50_D7_WT, subset = nFeature_RNA > 1000 & nFeature_RNA < 8500 & percent.mt < 10)
#add identifiers to the metadata if you want...helpful down the road
D50_D7_WT$sample <- NA
D50_D7_WT$time <- NA
D50_D7_WT$sample[which(str_detect(D50_D7_WT$orig.ident, "D7_WT1"))] <- "D7_WT1"
D50_D7_WT$sample[which(str_detect(D50_D7_WT$orig.ident, "D7_WT2"))] <- "D7_WT2"
D50_D7_WT$sample[which(str_detect(D50_D7_WT$orig.ident, "D50_WT1"))] <- "D50_WT1"
D50_D7_WT$sample[which(str_detect(D50_D7_WT$orig.ident, "D50_WT2"))] <- "D50_WT2"
D50_D7_WT$time[which(str_detect(D50_D7_WT$orig.ident, "WT"))] <- "Day7"
D50_D7_WT$time[which(str_detect(D50_D7_WT$orig.ident, "D50"))] <- "Day50"
#everything else is the same as above from here on out
#to separate samples in UMAP or PCA plots try split.by = "celltype" or group.by = "celltype" in the DimPlot lines

#######
#Seurat integration analysis
#######

#if you want to compare different samples (i.e. WT and KO, treated and untreated) you need to integrate the samples on a set of common anchors
#this is useful if you want to identify cell types present in both samples, or see how certain cell types change between samples

#first input, merge, and name your samples as described above
#now split the samples into a list of objects depending on how you want to compare
D50_split_seurat <- SplitObject(D50_WT_GKO_combined, split.by = "sample")
#normalize each sample separately
D50_split_seurat <- D50_split_seurat[c("WT1","WT2","GKO1",'GKO2')]
for (i in 1:length(D50_split_seurat)) {
  D50_split_seurat[[i]] <- NormalizeData(D50_split_seurat[[i]], verbose = TRUE)
  D50_split_seurat[[i]] <- FindVariableFeatures(D50_split_seurat[[i]], selection.method = "vst", nfeatures = 2000)}
#identify integration anchors
#i like to integrate on all genes in order to regress out cell cycle effects. if you don't want this just remove the features.to.integrate option.
all_genes <- as.vector(unique(rownames(D50_WT_GKO_combined@assays$RNA)))
D50_anchors <- FindIntegrationAnchors(object.list = D50_split_seurat, dims = 1:10) ###before 40 changed to 30
D50_integ <- IntegrateData(anchorset = D50_anchors, dims = 1:10, features.to.integrate = all_genes)
#now perform all of the usual steps on the integrated data
D50_integ <- CellCycleScoring(D50_integ, g2m.features=g2m_genes, s.features=s_genes)
D50_integ <- ScaleData(D50_integ, vars.to.regress = c("S.Score", "G2M.Score"))
D50_integ <- RunPCA(D50_integ, npcs = 30, verbose = FALSE)
PCAPlot(D50_integ,split.by = "sample") 
ElbowPlot(D50_integ)
D50_integ <- RunUMAP(D50_integ, reduction = "pca", dims = 1:10)  ##LJ used 10 dims and 0.15 res
D50_integ <- FindNeighbors(D50_integ, reduction = "pca", dims = 1:10)
D50_integ <- FindClusters(D50_integ, resolution = c(0.05, 0.075, 0.1, 0.125, 0.15, 0.2, 0.3))
#specify that we will perform downstream analyses on the corrected data
DefaultAssay(D50_integ) <- "integrated"
VlnPlot(D50_integ, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#specificy which dimension you want to plot
Idents(D50_integ) <- "integrated_snn_res.0.05"
DimPlot(D50_integ, reduction = "umap")
#can change colors of each cell cype if you want
DimPlot(D50_integ, reduction = "umap", group.by = "celltype", cols=c('#984EA3','#377EB8'))

#######
#RNA expression analyses
#######

#can use these lines on any Seurat object after UMAP analysis
DefaultAssay(D50_integ) <- "RNA"
#plot the expression of a transcript in your UMAP
FeaturePlot(D50_integ, features = c("PAX1"), max.cutoff = 3, cols = c("grey", "red"))
#plot the same thing as a violin plot
VlnPlot(D50_integ, features = c("TRIL"), split.by="celltype")
#find markers that are enriched in each cluster
D50_markers <- FindAllMarkers(object = D50_integ, only.pos = TRUE,logfc.threshold = 0.2)
write.csv(D50_markers,"D50_markers")
#write the top n markers for each cluster
D50_integ_top10 <- D50_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#make a heatmap of the top markers in each cluster
DefaultAssay(D50_integ) <- "integrated"
DoHeatmap(D50_integ, features = D50_integ_top10$gene,group.colors=cols) + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu"))) + guides(color=FALSE)

#######
#Module Scoring
#######

#if you want to score the enrichment of a set of genes in each cluster, you need to add module scores to the meta data
#first, import a list of genes as a text file
Gibbin_targets<-read.table("Gibbin_targets.txt",header=TRUE)
Gibbin_targets<-list("Gibbin_targets"=as.character(Gibbin_targets[,1]))
#add the scores to the metadata and plot the results
D50_integ<-AddModuleScore(D50_integ,feature= Gibbin_targets, name = "Gibbin_targets")
FeaturePlot(D50_integ, features = c("Gibbin_targets1"), cols = c("grey", "red"), min.cutoff = 0)
VlnPlot(D50_integ, features = c("Gibbin_targets1"), pt.size = 0)

