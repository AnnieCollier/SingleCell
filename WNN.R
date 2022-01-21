library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(ggplot2)
library(stringr)

#WNN tutorial from https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html#wnn-analysis-of-10x-multiome-rna-atac-1

# the 10x hdf5 file contains both data types. 
inputdata.10x <- Read10X_h5("WT_AC_filtered_feature_bc_matrix.h5")
# extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks
WT_multiome <- CreateSeuratObject(counts = rna_counts)
WT_multiome[["percent.mt"]] <- PercentageFeatureSet(WT_multiome, pattern = "^MT-")
# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg38"
frag.file <- "AC_atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations)
WT_multiome[["ATAC"]] <- chrom_assay
VlnPlot(WT_multiome, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE)
WT_multiome <- subset(
  x = WT_multiome,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20)
VlnPlot(WT_multiome, features = c("nCount_ATAC", "nCount_RNA","percent.mt"), ncol = 3,
        log = TRUE)


# RNA analysis
DefaultAssay(WT_multiome) <- "RNA"
WT_multiome <- SCTransform(WT_multiome, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')


# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(WT_multiome) <- "ATAC"
WT_multiome <- RunTFIDF(WT_multiome)
WT_multiome <- FindTopFeatures(WT_multiome, min.cutoff = 'q0')
WT_multiome <- RunSVD(WT_multiome)
WT_multiome <- RunUMAP(WT_multiome, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

WT_multiome <- FindMultiModalNeighbors(WT_multiome, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
WT_multiome <- RunUMAP(WT_multiome, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
WT_multiome <- FindClusters(WT_multiome, graph.name = "wsnn", algorithm = 3, verbose = FALSE, resolution = 0.125)

#WT_multiome$celltype <- Idents(WT_multiome)

p1 <- DimPlot(WT_multiome, reduction = "umap.rna", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(WT_multiome, reduction = "umap.atac", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(WT_multiome, reduction = "wnn.umap", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1
p2
p3
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

DefaultAssay(WT_multiome) <- "RNA"
VlnPlot(WT_multiome, features = c("ZIC2"))
FeaturePlot(WT_multiome, features = c("EPCAM"), cols = c("grey", "red"), min.cutoff = 0, reduction = "wnn.umap")
WT_multiome_markers <- FindAllMarkers(object = WT_multiome, only.pos = TRUE,logfc.threshold = 0.2)
top10_WT_multiome_markers <- WT_multiome_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

CoveragePlot(object=WT_multiome, region = 'TP63', features = 'TP63', assay = 'ATAC', expression.assay = 'SCT', peaks = FALSE)

library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

# Scan the DNA sequence of each peak for the presence of each motif, and create a Motif object
DefaultAssay(WT_multiome) <- "ATAC"
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
motif.matrix <- CreateMotifMatrix(features = granges(WT_multiome), pwm = pwm_set, genome = 'hg38', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
WT_multiome <- SetAssayData(WT_multiome, assay = 'ATAC', slot = 'motifs', new.data = motif.object)
# Note that this step can take 30-60 minutes 
WT_multiome <- RunChromVAR(
  object = WT_multiome,
  genome = BSgenome.Hsapiens.UCSC.hg38)
#plot the RNA expression and the motif on the WNN UMAP
motif.name <- ConvertMotifID(WT_multiome, name = 'EPCAM')
gene_plot <- FeaturePlot(WT_multiome, features = "sct_EPCAM", reduction = 'wnn.umap')
motif_plot <- FeaturePlot(WT_multiome, features = motif.name, min.cutoff = 0, cols = c("lightgrey", "darkred"), reduction = 'wnn.umap')
gene_plot | motif_plot
#find motifs that are enriched in each WNN cluster
markers_rna <- presto:::wilcoxauc.Seurat(X = WT_multiome, assay = 'data', seurat_assay = 'SCT')
markers_motifs <- presto:::wilcoxauc.Seurat(X = WT_multiome, assay = 'data', seurat_assay = 'chromvar')
motif.names <- markers_motifs$feature
colnames(markers_rna) <- paste0("RNA.", colnames(markers_rna))
colnames(markers_motifs) <- paste0("motif.", colnames(markers_motifs))
markers_rna$gene <- markers_rna$RNA.feature
markers_motifs$gene <- ConvertMotifID(WT_multiome, id = motif.names)
# a simple function to implement the procedure above
topTFs <- function(celltype, padj.cutoff = 1e-2) {
  ctmarkers_rna <- dplyr::filter(
    markers_rna, RNA.group == celltype, RNA.padj < padj.cutoff, RNA.logFC > 0) %>% 
    arrange(-RNA.auc)
  ctmarkers_motif <- dplyr::filter(
    markers_motifs, motif.group == celltype, motif.padj < padj.cutoff, motif.logFC > 0) %>% 
    arrange(-motif.auc)
  top_tfs <- inner_join(
    x = ctmarkers_rna[, c(2, 11, 6, 7)], 
    y = ctmarkers_motif[, c(2, 1, 11, 6, 7)], by = "gene"
  )
  top_tfs$avg_auc <- (top_tfs$RNA.auc + top_tfs$motif.auc) / 2
  top_tfs <- arrange(top_tfs, -avg_auc)
  return(top_tfs)}
# identify top markers in a given cluster and visualize
#plots the expression of the gene on the left and the accesbility of its motif on the right
head(topTFs("0"), 20)
write.csv(topTFs("1"),"topTFs.csv")
motif.name <- ConvertMotifID(WT_multiome, name = 'TWIST1')
gene_plot <- FeaturePlot(WT_multiome, features = "sct_TWIST1", reduction = 'wnn.umap')
motif_plot <- FeaturePlot(WT_multiome, features = motif.name, min.cutoff = 0, cols = c("lightgrey", "darkred"), reduction = 'wnn.umap')
gene_plot | motif_plot

