devtools::install_github("sqjin/CellChat")
devtools::install_github("jokergoo/ComplexHeatmap")
library(ComplexHeatmap)
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
library(systemfonts)
install.packages("systemfonts")
BiocManager::install("ComplexHeatmap")
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(ComplexHeatmap)
options(stringsAsFactors = FALSE)
library(patchwork)
library(umap)

#######
#Global Cell-Cell Communication Analysis
#######

#start the analysis with a Seurat object that has been clustered.
#make sure each cluster has a name
DimPlot(GKO_combined, reduction = "umap", label = TRUE)

#make a normalized data matrix
data.input <- GetAssayData(GKO_combined, assay = "RNA", slot = "data") 
labels <- Idents(GKO_combined)
#create a dataframe of the cell labels
identity <- data.frame(group = labels, row.names = names(labels))
cellchat_D7_GKO<- createCellChat(object = data.input, met=identity, group.by = "group")
cellchat_D7_GKO <- addMeta(cellchat_D7_GKO, meta = identity, meta.name = "labels")
#set "labels" as default cell identity
cellchat_D7_GKO <- setIdent(cellchat_D7_GKO, ident.use = "labels")
#show factor levels of the cell labels
levels(cellchat_D7_GKO@idents) 
#number of cells in each cell group
groupSize <- as.numeric(table(cellchat_D7_GKO@idents))

#use CellChatDB.human if running on human data
CellChatDB <- CellChatDB.human 
showDatabaseCategory(CellChatDB)
#Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

#there are three options for types of signaling to analyze. pick one and set that to your object
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact")
CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor")
cellchat_D7_GKO@DB <- CellChatDB.use # set the used database in the object

#these few steps are computer intense -- ~15 or so minutes 
#subset the expression data of signaling genes for saving computation cost
cellchat_D7_GKO <- subsetData(cellchat_D7_GKO)
future::plan("multiprocess", workers = 4) # do parallel
cellchat_D7_GKO <- identifyOverExpressedGenes(cellchat_D7_GKO)
cellchat_D7_GKO <- identifyOverExpressedInteractions(cellchat_D7_GKO)
cellchat_D7_GKO <- projectData(cellchat_D7_GKO, PPI.human)
cellchat_D7_GKO <- computeCommunProb(cellchat_D7_GKO)
cellchat_D7_GKO <- computeCommunProbPathway(cellchat_D7_GKO)
cellchat_D7_GKO <- aggregateNet(cellchat_D7_GKO)

#plot the interactions between your clusters
groupSize <- as.numeric(table(cellchat_D50_WT@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_D50_WT@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_D50_WT@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= T, title.name = "Interaction weights/strength")
mat <- cellchat_D50_WT@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = , edge.weight.max = max(mat), title.name = rownames(mat)[i])}

library(NMF)
#dot plot
#useful for global perspective on signaling
selectK(cellchat_WT_d50, pattern = "outgoing") #use this to figure out nPattern number to input below
nPatterns = 3
cellchat_WT_d7 <- identifyCommunicationPatterns(cellchat_WT_d7, pattern = "outgoing", k = nPatterns, heatmap.show = FALSE)
netAnalysis_dot(cellchat_WT_d7, pattern = "outgoing")
netAnalysis_river(cellchat_WT_d50, pattern = "outgoing")
#you can run this same analysis for incoming signaling patterns
selectK(cellchat_WT_d7, pattern = "incoming")
nPatterns = 3  
cellchat_WT_d7 <- identifyCommunicationPatterns(cellchat_WT_d7, pattern = "incoming", k = nPatterns, heatmap.show = FALSE)
netAnalysis_dot(cellchat_WT_d7, pattern = "incoming")
netAnalysis_river(cellchat_WT_d7, pattern = "incoming")

#######
#Pathway Analysis
#######

#now we can zero in on specific signaling pathways
#there are a ton of different pathways you can look at, you need to refrence the website to get a list of them 
pathways.show <- c("ncWNT") 
#a numeric vector, you have to kind of tweak this a bit depending on how many clusters you have 
vertex.receiver = seq(1,2) 
# Hierarchy plot
netVisual_aggregate(cellchat_WT_d50, signaling = pathways.show,  vertex.receiver = vertex.receiver, vertex.size = groupSize)
# Circle plot
netVisual_aggregate(cellchat_WT_d7, signaling = pathways.show, layout = "circle")
#What is the contribution of that patway you're looking at 
netAnalysis_contribution(cellchat_WT_d7, signaling = pathways.show)
#chord plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
#heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat,color.heatmap = "Reds")
#gene expression plot
plotGeneExpression(cellchat_WT_d50, signaling = "ACTIVIN", enriched.only=TRUE)

#Similar to above, this section is dependent on what pathway you put in but gives you some more anlaysis of things 
#the slot 'netP' means the inferred intercellular communication network of signaling pathways
cellchat_D7_GKO <- netAnalysis_computeCentrality(cellchat_D7_GKO, slot.name = "netP")
netAnalysis_signalingRole_network(cellchat_WT_d7, signaling = "WNT")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")

#######
#Comparative Analysis
#######

#if you want to see how signaling changes between to cell types or treamtent groups, we need to marge the objects
#first, perform cellchat on each sample separately
#then, merge together
object.list <- list(WT = cellchat_D7_WT, GKO = cellchat_D7_GKO)
cellchat_merge_d7 <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix=TRUE)
#plot the number and strength of interactions between cell types
gg1 <- compareInteractions(cellchat_merge_d7, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat_merge_d7, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
#circle plot
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat_merge_d7, weight.scale = T, label.edge= T)
netVisual_diffInteraction(cellchat_merge_d7, weight.scale = T, label.edge= T, measure = "weight")
#heatmap
gg1 <- netVisual_heatmap(cellchat_merge_d7)
gg2 <- netVisual_heatmap(cellchat_merge_d7, measure = "weight")
gg1 + gg2
#See if any signaling pathways show major functional changes between cell types, i.e. changes in senders or receivers (for example in the GKO the sender swtiches from mesoderm to ectoderm)
#structural = changes in signaling network wihtout considering senders and receivers (the ligands and receptors are changing their pattern)
cellchat_merge_d7 <- computeNetSimilarityPairwise(cellchat_merge_d7, type = "functional")
cellchat_merge_d7 <- netEmbedding(cellchat_merge_d7, type = "functional")
cellchat_merge_d7 <- netClustering(cellchat_merge_d7, type = "functional")
netVisual_embeddingPairwise(cellchat_merge_d7, type = "functional")
netVisual_embeddingPairwiseZoomIn(cellchat_merge_d7, type = "functional", nCol = 2)
rankSimilarity(cellchat_merge_d7, type = "functional")
cellchat_merge <- computeNetSimilarityPairwise(cellchat_merge, type = "structural")
cellchat_merge <- netEmbedding(cellchat_merge, type = "structural")
cellchat_merge <- netClustering(cellchat_merge, type = "structural")
netVisual_embeddingPairwise(cellchat_merge, type = "structural")
netVisual_embeddingPairwiseZoomIn(cellchat_merge, type = "structural", nCol = 2)
rankSimilarity(cellchat_merge, type = "structural")
#compare outgoing and incoming interaction strengths
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)}
patchwork::wrap_plots(plots = gg)
#information flow plots
gg1 <- rankNet(cellchat_merge_d7, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat_merge_d7, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2
#plot gene expression
cellchat_merge_d50@meta$datasets = factor(cellchat_merge_d50@meta$datasets, levels = c("WT", "GKO")) # set factor level
plotGeneExpression(cellchat_merge_d7, signaling = "CXCL", split.by = "datasets")



