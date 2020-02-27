# calculate diffusion distance
# this function uses the 'destiny' package
# data - nxm matrix where n is number of samples, m is number of features
# n - number of timesteps to spread diffusion
# n_neigh - number of nearest neighbours taken into account
# n_local - which nearest neighbour will determine the local sigma
getDifusionDist <- function(data, n = 16, n_neigh = 10, n_local = seq(to = min(7, n_neigh))) {
  knn <- destiny:::get_knn(data, NULL, n_neigh)
  sigmas <- destiny:::get_sigmas(data, knn$dist, "local", n_local)
  sigmas <- sigmas@optimal_sigma
  
  L <- destiny:::no_censoring(knn$dist_mat, sigmas)
  L[cbind(1:nrow(L), 1:nrow(L))] <- 1
  
  d <- rowSums(L, na.rm = TRUE)
  
  M <- L / d
  M <- drop0(M)
  
  max_d <- floor(log2(n))
  ds <- list()
  ds$d0 <- M
  
  res <- Diagonal(nrow(M))
  for(i in 1:max_d)
    ds[[paste0("d", i)]] <- ds[[paste0("d", (i - 1))]] %*% ds[[paste0("d", (i - 1))]]
  
  while(n > 0) {
    deg <- floor(log2(n))
    res <- res %*% ds[[paste0("d", deg)]]
    n <- n - 2^deg
  }
  
  Dist(as.matrix(t(t(res) / sqrt(d/sum(d)))))
}

library(readr)
library(Seurat)
library(irlba)
library(Rfast)
library(cowplot)
library(sleepwalk)

e13_A <- read_rds("data/e13_A.rds")
normAndVargenes <- function(obj) {
  obj <- NormalizeData(obj)
  FindVariableGenes(obj, do.plot = F, y.cutoff = 0.5)
}

e13_A <- normAndVargenes(e13_A)
um13_A <- read_rds("data/e13A_umap.rds") #use embedding from the figure C
data13_A <- t(e13_A@data[e13_A@var.genes, ]) #use only variable genes

pca <- prcomp_irlba(data13_A, n = 50)

#distance is rounded to make the resulting html file smaller
eu_dist <- round(Dist(pca$x), 0)
eu_dist_all <- round(Dist(as.matrix(data13_A)), 0)

difD_all <- round(getDifusionDist(data13_A), 1)
difD <- round(getDifusionDist(pca$x), 1)


# cosine distance 
# cos_dist <- acos(apply(pca$x, 1, function(row) {
#   colSums(row * t(pca$x))/sqrt(sum(row^2) * rowSums(pca$x^2))
# }))/pi

sleepwalk(list(um13_A, um13_A, um13_A, um13_A), 
          distances = list(eu_dist_all, eu_dist, difD_all, difD), 
          compare = "distances", 
          titles = c("Euclidean distance", "Euclidean distance after PCA", "Diffusion distance", "Diffusion distance after PCA"),
          saveToFile = "supplement/src/Fig_F.html")

file <- readLines("supplement/src/Fig_F.html")
line <- which(grepl("^set_up_chart\\(", file))
newFile <- c(file[1:(line - 1)], 
             "width = window.innerWidth/(n_charts) * 0.9 * 2;",
             file[line:length(file)])
writeLines(newFile, "supplement/src/Fig_F.html")

sleepwalk(list(um13_A, um13_A, um13_A, um13_A), 
          distances = list(eu_dist_all, eu_dist, difD_all, difD), 
          compare = "distances", 
          titles = c("Euclidean distance", "Euclidean distance after PCA", "Diffusion distance", "Diffusion distance after PCA"),
          maxdists = c(40, 27, 2.7, 4.5))

cell <- 1354
plots <- slw_snapshot(cell, 3, returnList = TRUE)

vert_um <- (max(um13_A[, 2]) - min(um13_A[, 2]))/25
hor_um <- (max(um13_A[, 1]) - min(um13_A[, 1]))/25

plots <- lapply(plots, function(plot) {
  plot + geom_segment(aes(x = um13_A[cell, 1] + 1.2 * hor_um, xend = um13_A[cell, 1] + 0.2 * hor_um, 
                   y = um13_A[cell, 2] + 1.2 * vert_um, yend = um13_A[cell, 2] + 0.2 * vert_um), size = 1,
               colour = "red", arrow = arrow(length = unit(0.03, "npc")))
})
figF <- plot_grid(plotlist = plots, labels = c("A", "B", "C", "D"))

png(filename = "figures/Fig_F.png", 
    width = 900, height = 1000)
print(figF)
dev.off()
  