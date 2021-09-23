# SingleCell
All scripts related to 10X Genomics single cell analysis

Use the Seurat script to process CellRanger output files. Use the "filtered_feature_bc_matrix" folder found in the "outs" directory of the CellRanger output.
The Seurat script performs the following tasks on single, merged, or integrated datasets:
- Clustering and UMAPs
- Transcript expression quantification
- Module scoring

![D7_UMAP](https://user-images.githubusercontent.com/90862478/134575326-e0671d0c-5c8d-47d1-8e83-712fe51a053d.png)



Use the CellChat script to perform cell-cell interaction analysis. Using this script you can:
- Quantify interactions between clusters
- Identify signaling pathways between clusters
- Quantify changes in signaling between two cell types or treatment groups

Use the Monocle2 script to perform pseudotime analysis![Uploading D7_UMAP.png…]()
