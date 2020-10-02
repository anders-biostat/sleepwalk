library(Seurat)
library(tidyverse)

cbmc.rna <- as.sparse(read.csv(file = "~/Downloads/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz", sep = ",", 
                               header = TRUE, row.names = 1))
cbmc.rna <- CollapseSpeciesExpressionMatrix(cbmc.rna)
cbmc.adt <- as.sparse(read.csv(file = "~/Downloads/GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz", sep = ",", 
                               header = TRUE, row.names = 1))
cbmc.adt <- cbmc.adt[setdiff(rownames(x = cbmc.adt), c("CCR5", "CCR7", "CD10")), ]
cbmc <- CreateSeuratObject(counts = cbmc.rna)
cbmc <- NormalizeData(cbmc)
cbmc <- FindVariableFeatures(cbmc)
cbmc <- ScaleData(cbmc)
cbmc <- RunPCA(cbmc, verbose = FALSE)

cbmc <- FindNeighbors(cbmc, dims = 1:25)
cbmc <- FindClusters(cbmc, resolution = 0.8)
cbmc <- RunTSNE(cbmc, dims = 1:25, method = "FIt-SNE")

cbmc.rna.markers <- FindAllMarkers(cbmc, max.cells.per.ident = 100, min.diff.pct = 0.3, only.pos = TRUE)
new.cluster.ids <- c("Memory CD4 T", "CD14+ Mono", "Naive CD4 T", "NK", "CD14+ Mono", "Mouse", "B", 
                     "CD8 T", "CD16+ Mono", "T/Mono doublets", "NK", "CD34+", "Multiplets", "Mouse", "Eryth", "Mk", 
                     "Mouse", "DC", "pDCs")
names(new.cluster.ids) <- levels(cbmc)
cbmc <- RenameIdents(cbmc, new.cluster.ids)

cbmc <- RunUMAP(cbmc, dims = 1:25)

cite_pca <- Embeddings(Reductions(cbmc, "pca"))
cite_tsne <- Embeddings(Reductions(cbmc, "tsne"))
cite_umap <- Embeddings(Reductions(cbmc, "umap"))
cite_data <- Matrix::t(GetAssayData(Assays(cbmc, "RNA"))[VariableFeatures(cbmc), ])

save(cite_pca, cite_umap, cite_tsne, cite_data, file = "cite_seurat.rda")
