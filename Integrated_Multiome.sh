####just analyze RNA portion
input_files="~/path/to/filtered_feature_bc_matrix/"
WT1.data <- Read10X(data.dir = input_files)
WT1.data <- CreateSeuratObject(counts = WT1.data$`Gene Expression`, min.features = 100, min.cells = 5, project = "D7_WT")
WT1.data

input_files="~/path/to/filtered_feature_bc_matrix/"
KO.data <- Read10X(data.dir = input_files)
KO.data <- CreateSeuratObject(counts = KO.data$`Gene Expression`, min.features = 100, min.cells = 5, project = "D7_KO")
KO.data

#integrate datasets
WT_KO_combined <- merge(x = WT1.data,y = GRHL2_KO.data, add.cell.ids = c("WT", "KO"))
WT_KO_combined[["percent.mt"]] <- PercentageFeatureSet(WT_KO_combined, pattern = "^MT-")
VlnPlot(WT_KO_combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
WT_KO_combined <- subset(WT_KO_combined, subset = nFeature_RNA > 1500 & nFeature_RNA < 8500 & percent.mt < 20)
WT_KO_combined$sample <- NA
WT_KO_combined$celltype <- NA
WT_KO_combined$sample[which(str_detect(WT_KO_combined$orig.ident, "WT"))] <- "WT"
WT_KO_combined$sample[which(str_detect(WT_KO_combined$orig.ident, "KO"))] <- "KO"
WT_KO_combined$celltype[which(str_detect(WT_KO_combined$orig.ident, "WT"))] <- "WT"
WT_KO_combined$celltype[which(str_detect(WT_KO_combined$orig.ident, "KO"))] <- "KO"
split_seurat <- SplitObject(WT_KO_combined, split.by = "sample")
split_seurat <- split_seurat[c("WT","KO")]
for (i in 1:length(split_seurat)) {
  split_seurat[[i]] <- NormalizeData(split_seurat[[i]], verbose = TRUE)
  split_seurat[[i]] <- FindVariableFeatures(split_seurat[[i]], selection.method = "vst", nfeatures = 2000) }
all_genes <- as.vector(unique(rownames(WT_KO_combined@assays$RNA)))
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, dims = 1:10) ###before 40 changed to 30
integ <- IntegrateData(anchorset = integ_anchors, dims = 1:10, features.to.integrate = all_genes)
DefaultAssay(integ) <- "integrated"
s_genes <- cc.genes$s.genes
g2m_genes <- cc.genes$g2m.genes
integ <- CellCycleScoring(integ, g2m.features=g2m_genes, s.features=s_genes)
integ <- ScaleData(integ, vars.to.regress = c("S.Score", "G2M.Score"))
integ <- RunPCA(integ, npcs = 30, verbose = FALSE)
PCAPlot(integ)
integ <- RunUMAP(integ, reduction = "pca", dims = 1:10)
integ <- FindNeighbors(integ, reduction = "pca", dims = 1:10)
integ <- FindClusters(integ, resolution = 0.05)
DimPlot(integ, reduction = "umap")
DimPlot(integ, reduction = "umap", group.by="celltype")
DimPlot(integ, reduction = "umap", split.by="celltype")
integ <- subset(integ, idents = c(0, 1, 2, 4))
new.cluster.ids <- c("cluster1","cluster2","cluster3)
names(new.cluster.ids) <- levels(integ)
integ <- RenameIdents(integ, new.cluster.ids)
PCAPlot(integ, group.by="celltype")
PCAPlot(integ, split.by="celltype")

integ_WT <- subset(integ, celltype == "WT")
integ_KO <- subset(integ, celltype == "KO")


DefaultAssay(integ) <- "RNA"
FeaturePlot(integ, features = c("gene"), max.cutoff = 3, cols = c("grey", "red"), split.by = "celltype")
VlnPlot(integ, features = c("gene"), split.by="celltype")
markers_integ <- FindAllMarkers(object = integ, only.pos = TRUE,logfc.threshold = 0.2)
top10 <- markers_integ %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(integ, features = top10$gene,group.colors=cols) + scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdYlBu"))) + guides(color=FALSE)

#number and proportion of cells in each cluster
table(integ$"active.ident")
table(Idents(integ))
prop.table(table(Idents(integ), integ$sample), 2)


#make a csv of UMAP projections or clusters for the Loupe browser
write.csv(integ@reductions$umap@cell.embeddings, file = paste("integ_loupe.csv",sep=""))
write.csv(integ@active.ident, file = paste("integ_loupe_clusters.csv",sep=""))
