# SingleCell
All scripts related to 10X Genomics single cell analysis

Use the Seurat script to process CellRanger output files. Use the "filtered_feature_bc_matrix" folder found in the "outs" directory of the CellRanger output.
The Seurat script performs the following tasks on single, merged, or integrated datasets:
- Clustering and UMAPs
- Transcript expression quantification
- Module scoring

![D7_UMAP](https://user-images.githubusercontent.com/90862478/134575326-e0671d0c-5c8d-47d1-8e83-712fe51a053d.png)

![GRHL2_D7_expression](https://user-images.githubusercontent.com/90862478/134575838-55f188ea-24b8-4a81-babe-a901c0332ee8.png)


Use the CellChat script to perform cell-cell interaction analysis. Using this script you can:
- Quantify interactions between clusters
- Identify signaling pathways between clusters
- Quantify changes in signaling between two cell types or treatment groups

![cell_cell_contact_heatmap_D50_diff](https://user-images.githubusercontent.com/90862478/134576096-d021a790-f875-4486-8aa5-5de116da265b.png)
![circle_changes_cellchat](https://user-images.githubusercontent.com/90862478/134576187-49a04ea4-1661-4c58-bc3c-2dde7b042c0c.png)


Use the Monocle2 script to perform pseudotime analysis
![WT_pseudo_ecto](https://user-images.githubusercontent.com/90862478/134575981-9df7256d-e36b-47de-80be-bef37a7f30c3.png)
