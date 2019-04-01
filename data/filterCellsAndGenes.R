library(Seurat)
library(biomaRt)
library(tidyverse)

saveAsRDS <-function(dir) {
  exprs <- t(Read10X(data.dir = dir))  
  
  #remove empty droplets (0 UMIs) and non-exepressed genes
  exprs <- exprs[rowSums(exprs) > 0, colSums(exprs) > 0]
  
  ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  useDataset("mmusculus_gene_ensembl", mart = ensembl)
  mt_genes <- intersect((select(ensembl, keys = "MT", keytype = "chromosome_name", 
                     columns = "mgi_symbol"))$mgi_symbol, colnames(exprs))
  mt_frac <- rowSums(exprs[, mt_genes])/rowSums(exprs)
  #remove cells with more than 10% mt genes
  exprs <- exprs[mt_frac < 0.1, ]
  
  #remove mt genes
  exprs <- exprs[, !(colnames(exprs) %in% mt_genes)]
  
  #remove ribosomal genes
  exprs <- exprs[, !grepl("^M?RP[SL]", colnames(exprs), ignore.case = TRUE)]

  #remove cells with less than 3500 or more than 15000 UMIs
  allCounts <- rowSums(exprs)
  exprs <- exprs[allCounts > 3500 & allCounts < 15000, ]
  
  exprs <- exprs[, colSums(exprs) > 0]
  
  so <- CreateSeuratObject(raw.data = t(exprs))
  
  write_rds(so, paste0(dir, ".rds"), compress = "gz")
}

#saveAsRDS("e11_A")
saveAsRDS("e13_A")
saveAsRDS("e13_B")
#saveAsRDS("e15_A")
saveAsRDS("e14_A")
