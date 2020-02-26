library( Rtsne )
library( irlba )
library( umap )
library( sleepwalk )
# devtools::load_all("w/repos/sleepwalk/")

# These two file can be downloaded from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE100866
countsRNA_filename <- "~/Downloads/GSE100866_CBMC_8K_13AB_10X-RNA_umi.csv.gz"
countsADT_filename <- "~/Downloads/GSE100866_CBMC_8K_13AB_10X-ADT_umi.csv.gz"

# Load the file. (This takes a while as it is a large file.)
countsRNA <- as.matrix( read.csv( gzfile( countsRNA_filename ), row.names = 1) )
countsADT <- as.matrix( read.csv( gzfile( countsADT_filename ), row.names = 1) )


# Calculate for each cell ratio of molecules mapped to human genes
# versus mapped to mouse genes
human_mouse_ratio <- 
  colSums( countsRNA[ grepl( "HUMAN" , rownames(countsRNA) ), ] ) / 
  colSums( countsRNA[ grepl( "MOUSE" , rownames(countsRNA) ), ] )

# Keep only the cells with at least 10 times more human than mouse genes
# and keep only the counts for the human genes.
countsRNA <- countsRNA[ 
  grepl( "HUMAN" , rownames(countsRNA) ), 
  human_mouse_ratio > 10 ]

# Remove the "HUMAN_" prefix from the gene names
rownames(countsRNA) <- sub( "HUMAN_", "", rownames(countsRNA) )

# Subset the ADT matrix to the same cells as in the RNA matrix
countsADT <- countsADT[ , colnames(countsRNA) ]


# Calculate size factors
exprsRNA <- matrix( nrow = nrow(countsRNA), ncol = ncol(countsRNA), dimnames = dimnames(countsRNA) )
for( j in seq.int( ncol(countsRNA) ) )
  exprsRNA[,j] <- log2( 1 + countsRNA[,j] / sum(countsRNA[,j]) )

# Run a PCA on the expression data
pca <- prcomp_irlba( t(exprsRNA), n=50 )

# Set the rownames manually (as IRLBA doesn't do that for us)
rownames(pca$x) <- colnames(exprsRNA)
rownames(pca$rotation) <- rownames(exprsRNA)

# Run t-SNE on the data
tsneRNA <- Rtsne( pca$x, pca=FALSE, verbose=TRUE )
rownames(tsneRNA$Y) <- rownames(pca$x)

# Explore the result with sleepwalk
sleepwalk( tsneRNA$Y, pca$x, 0.07 )

# How much difference did it make that we allowed so many PCs? Let's try with only 10
pca_B <- prcomp_irlba( t(exprsRNA), n=10 )
tsneRNA_B <- Rtsne( pca_B$x, pca=FALSE, verbose=TRUE )

# And compare side by side
sleepwalk( 
  list( tsneRNA$Y, tsneRNA_B$Y ), 
  list( pca$x, pca_B$x ), 
  c( 0.07, 0.07 ) )

# How about UMAP?
umapres <- umap( pca$x )

# Now compare these
sleepwalk( 
  list( tsneRNA$Y, umapres$layout), 
  pca$x,
  0.07)

#norms <- pca$x/sqrt(rowSums(pca$x^2))

sleepwalk(tsneRNA$Y,
          list(norms, norms), 
          metric = c("euclid", "cosine"),
          compare = "distances")

# Does one actually need to normalize?
pca_U <- prcomp_irlba( t(log2(countsRNA+1)), n=10 )
tsneRNA_U <- Rtsne( pca_U$x, pca=FALSE, verbose=TRUE )

sleepwalk( 
  list( tsneRNA_B$Y, tsneRNA_U$Y ), 
  list( pca_B$x, pca_U$x ) )


# We still have the ADT data
tsneADT <- Rtsne( t( log10( 1 + countsADT ) ), pca=FALSE, verbose=TRUE )

sleepwalk( 
  list( tsneRNA$Y, tsneADT$Y ), 
  list( pca$x, t( log10( 1 + countsADT ) ) ), 
  c( 0.07, 2 ) )

# It is possible to use distances instead of features 
# Let's look at genes instead of cells
means <- apply( exprsRNA, 1, mean )
vars <- apply( exprsRNA, 1, var )
plot( means, vars/means, pch = ".", log = "xy" )
abline( h = 1e-3 )
goodGenes <- names( which( vars/means > 1e-3 ) )

distGenes <- acos( cor( t(exprsRNA[ goodGenes, ]), method = "spearman") ) / pi
tsneDual <- Rtsne( distGenes, is_distance = TRUE )

sleepwalk( tsneDual$Y, distances = distGenes, maxdists = 0.7, pointSize = 3)

# To demonstrate the other multi view mode, let's divide up the RNA data
cellGroup1 <- sample( colnames(countsRNA), ncol(countsRNA)/3 )
cellGroup2 <- setdiff( colnames(countsRNA), cellGroup1 )
sleepwalk( 
  list( tsneRNA$Y[cellGroup1,], tsneRNA$Y[cellGroup2,] ), 
  list( pca$x[cellGroup1,], pca$x[cellGroup2,] ), 
  .07,
  same = "features" )

tsne_G1 <- Rtsne( pca$x[cellGroup1,], pca = FALSE, verbose=TRUE )
tsne_G2 <- Rtsne( pca$x[cellGroup2,], pca = FALSE, verbose=TRUE )

sleepwalk( 
  list( tsneRNA$Y, tsne_G1$Y, tsne_G2$Y ), 
  list( pca$x, pca$x[cellGroup1,], pca$x[cellGroup2,] ), 
  c(.05, .05, 0.01),
  same = "features")

# now the sleepwalk with dual 
means <- apply( exprsRNA, 1, mean )
vars <- apply( exprsRNA, 1, var )
plot( means, vars/means, pch=".", log="xy" )
abline( h=1e-3 )
goodGenes <- names( which( vars/means > 1e-3 ) )

distGenes <- acos( cor( t(exprsRNA[ goodGenes, ]), method="spearman") ) / pi
tsneDual <- Rtsne( acos( genecors )/pi, is_distance = TRUE )

sleepwalk( tsneDual$Y, , .003 )
