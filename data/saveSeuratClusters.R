#copy-pasted from https://satijalab.org/seurat/multimodal_vignette.html

library(Seurat)

cbmc.rna <- read.csv("~/Downloads/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz", 
                     sep = ",", header = TRUE, row.names = 1)

# To make life a bit easier going forward, we're going to discard all but
# the top 100 most highly expressed mouse genes, and remove the 'HUMAN_'
# from the CITE-seq prefix
cbmc.rna.collapsed <- CollapseSpeciesExpressionMatrix(cbmc.rna)

# Load in the ADT UMI matrix
cbmc.adt <- read.csv("~/Downloads/GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz", 
                     sep = ",", header = TRUE, row.names = 1)

# To avoid any confusion where genes and proteins might have the same name,
# we'll append 'CITE_' to each of the ADT rownames. This is not strictly
# necessary, but it helps for clarity
cbmc.citeseq <- cbmc.adt
rownames(cbmc.citeseq) <- paste0("CITE_", rownames(cbmc.adt))

# Lastly, we observed poor enrichments for CCR5, CCR7, and CD10 - and
# therefore remove them from the matrix (optional)
cbmc.citeseq <- cbmc.citeseq[setdiff(rownames(cbmc.citeseq), c("CITE_CCR5", 
                                                               "CITE_CCR7", "CITE_CD10")), ]

cbmc <- CreateSeuratObject(raw.data = cbmc.rna.collapsed)

# standard log-normalization
cbmc <- NormalizeData(cbmc)

# choose ~1k variable genes
cbmc <- FindVariableGenes(cbmc, do.plot = FALSE, y.cutoff = 0.5)

# standard scaling (no regression)
cbmc <- ScaleData(cbmc, display.progress = FALSE)

# Run PCA, select 13 PCs for tSNE visualization and graph-based clustering
cbmc <- RunPCA(cbmc, pcs.print = 0)

cbmc <- FindClusters(cbmc, dims.use = 1:13, print.output = FALSE)
cbmc <- RunTSNE(cbmc, dims.use = 1:13)

cbmc <- RunUMAP(cbmc, dims.use = 1:13, spread = 2)

current.cluster.ids <- 0:15
# Note, for simplicity we are merging two CD14+ Mono clusters (that differ
# in the expression of HLA-DR genes), and two NK clusters (that differ in
# cell cycle stage)
new.cluster.ids <- c("CD4 T", "CD14+ Mono", "CD14+ Mono", "NK", "Mouse", "B", 
                     "CD8 T", "CD16+ Mono", "Unknown", "CD34+", "Mk", "Eryth", "DC", "Mouse", 
                     "pDC", "NK")
cbmc@ident <- plyr::mapvalues(x = cbmc@ident, from = current.cluster.ids, to = new.cluster.ids)

# We will define a CITE assay, and store raw data for it. Note that it's
# convenient, but not required, to use the same name as the rowname prefix
# we defined earlier.

# If you are interested in how these data are internally stored, you can
# check out the @assay slot, and the assay class, which is defined in
# multimodal.R Note that RNA data is still stored in its normal slots, but
# can also be accessed using GetAssayData and SetAssayData, using the 'RNA'
# assay

cbmc <- SetAssayData(cbmc, assay.type = "CITE", slot = "raw.data", new.data = cbmc.citeseq)

# Now we can repeat the preprocessing (normalization and scaling) steps that
# we typically run with RNA, but modifying the 'assay.type' argument.  For
# CITE-seq data, we do not recommend typical LogNormalization. Instead, we
# use a centered log-ratio (CLR) normalization, computed independently for
# each gene.  This is a slightly improved procedure from the original
# publication, and we will release more advanced versions of CITE-seq
# normalizations soon.
cbmc <- NormalizeData(cbmc, assay.type = "CITE", normalization.method = "genesCLR")
cbmc <- ScaleData(cbmc, assay.type = "CITE", display.progress = FALSE)

save(cbmc, file = "cite_seurat.rda")
