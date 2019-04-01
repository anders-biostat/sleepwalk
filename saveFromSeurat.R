library(Seurat)

load(file = "data/cite_seurat.rda")

citeSeq <- list(ndims = 13, pca = cbmc@dr$pca@cell.embeddings, tsne = cbmc@dr$tsne@cell.embeddings,
                umap = cbmc@dr$umap@cell.embeddings, clusters = cbmc@ident, adt = t(cbmc@assay$CITE@scale.data))

save(citeSeq, file = "data/cite_data.rda")
