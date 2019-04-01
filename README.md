# sleepwalk
Exploring dimension-reduced embeddings

## Code for paper

This branch of the repository contains the R code for all the figures and for the supplementary material of the Sleepwalk paper.

## List of scripts

* `saveSeuratClusters.R` - reads in and normalizes the cord-blood single-cell data, performs clustering, calculates t-SNE and 
UMAP embeddings. Mostly copied from the Seurat [multi-modal vignette](https://satijalab.org/seurat/multimodal_vignette.html).
Before running this script one needs to download the 
[ADT](ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz) 
and [RNA](ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE100nnn/GSE100866/suppl/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz) data.
Saves a Seurat object for further use. The resulting object can be downloaded from [figshare.com](10.6084/m9.figshare.7908059).

* `saveFromSeurat.R` - reads in the [Seurat object](10.6084/m9.figshare.7908059) and saves only the necessary data for
further use. The resulting list can be downloaded from [figshare.com](10.6084/m9.figshare.7910504).

* `filterCellsAndGenes.R` - reads the CellRanger output files of the developing murine cerebellum data, performs genes and
cells filtering, and saves the result as several Seurat objects. Input files as well as the resulting objects can be
downloaded from [figshare.com](10.6084/m9.figshare.7910483).

* `Fig_A.R` - generates Figure 1 and also saves it to be used in the supplement. Requires `cite_data.rda` (can be 
downloaded [here](10.6084/m9.figshare.7910504)).

* `Fig_B.R` - generates Figure 2 and corresponding life version for the supplement. Requires `cite_data.rda` (can be 
downloaded [here](10.6084/m9.figshare.7910504)).

* `Fig_C.R` - generates Figure 3 and corresponding life version for the supplement. Requires `e13_A.rds` (can be 
downloaded [here](10.6084/m9.figshare.7910483)). Saves a UMAP embedding for further use.

* `Fig_Ca.R` - generates a UMAP-t-SNE comparison for cord-blood data. Requires `cite_data.rda` (can be 
downloaded [here](10.6084/m9.figshare.7910504)).

* `Fig_D.R` - generates Figure 4, Supplement Figure 1 and corresponding life version for the supplement. Requires input files from 
[here](10.6084/m9.figshare.7910483) as well as the UMAP embedding from `fig_C.R` (can be downloaded 
[here](10.6084/m9.figshare.7910504)).

* `Fig_E.R` - generates a distance comparison for cord-blood data. Requires `cite_data.rda` (can be 
downloaded [here](10.6084/m9.figshare.7910504)).
