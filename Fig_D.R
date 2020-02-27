library(Seurat) # v2.3.4
library(tidyverse)

e13_A <- read_rds("data/e13_A.rds")
e13_B <- read_rds("data/e13_B.rds")
e14 <- read_rds("data/e14_A.rds")

normAndVargenes <- function(obj) {
  obj <- NormalizeData(obj)
  FindVariableGenes(obj, do.plot = F, y.cutoff = 0.5)
}

e13_A <- normAndVargenes(e13_A)
e13_B <- normAndVargenes(e13_B)
e14 <- normAndVargenes(e14)

data13_A <- t(e13_A@data[e13_A@var.genes, ])
data13_B <- t(e13_B@data[e13_B@var.genes, ])
data14 <- t(e14@data[e14@var.genes, ])

library(irlba)
library(uwot)

getEmbedding <- function(data) {
  set.seed(55555)
  pca <- prcomp_irlba(data, n = 50)
  umap(pca$x, spread = 7, n_sgd_threads = 1)
}

um13_A <- read_rds("data/e13A_umap.rds") #use embedding from the figure C
um13_B <- getEmbedding(data13_B)
um14 <- getEmbedding(data14)

commonGenes <- intersect(intersect(e13_A@var.genes, e13_B@var.genes), e14@var.genes)

set.seed(55555)
pca <- prcomp_irlba(rbind(data13_A[, commonGenes], data13_B[, commonGenes], data14[, commonGenes]), n = 50)
comFeatures13_A <- data13_A[, commonGenes] %*% pca$rotation
comFeatures13_B <- data13_B[, commonGenes] %*% pca$rotation
comFeatures14 <- data14[, commonGenes] %*% pca$rotation

library(sleepwalk)

sleepwalk(list(um13_A, um13_B, um14), list(comFeatures13_A, comFeatures13_B, comFeatures14), same = "features", pointSize = 2.5,
          titles = c("e13.5_A", "e13.5_B", "e14.5"), nrow = 1)

cell <- 3105
plots <- slw_snapshot(cell, 2, TRUE)

library(ggforce)
library(r2d2)

#use "Msx3" for progenitors, "Meis2" for Glutamatergic cells, "Lhx5" for GABAergic progenitors
#classification is very unprecise, required only to ruoghly show the cell types
#we'll leave only cells that are close to center of mass of the corresponding cluster
getTypes <- function(data, emb) {
  typeMarkers <- cbind(data[, "Msx3"] > 2, data[, "Meis2"] > 1.5, data[, "Lhx5"] > 1)
  types <- apply(typeMarkers, 1, function(row) {
    if(sum(row) == 0) return("unknown")
    if(sum(row) > 1) return("unknown")
    c("Early progenitors", "Glutamatergic", "GABAergic")[row]
  })
  
  means <- apply(emb, 2, function(col) {
    tapply(col, as.factor(types), mean)
  })
  for(type in rownames(means)) {
    confs <- conf2d(emb[types == type, 1], emb[types == type, 2], level = 0.9, plot = F)
    types[types == type][!confs$inside] <- "unknown"
  #  dists <- rowSums(t(t(emb[types == type, ]) - means[type, ])^2)
  #  cutoff <- quantile(dists, probs = 0.95)
  #  types[types == type][dists > cutoff] <- "unknown"
  }
  types
}

types13_A <- getTypes(data13_A, um13_A)
types13_B <- getTypes(data13_B, um13_B)
types14 <- getTypes(data14, um14)

addLabels <- function(plot, types, emb, arrow = FALSE) {
  vert <- (max(emb[, 2]) - min(emb[, 2]))/25
  hor <- (max(emb[, 1]) - min(emb[, 1]))/25
  
#  plot$layers <- c(geom_mark_hull(aes(x = emb[, 1], y = emb[, 2], group = types, label = types, filter = types != "unknown"), 
#                                     linetype = 2, concavity = 3, size = 1, colour = "#555555", label.fontface = "plain",
#                                  label.fontsize = 20, label.buffer = unit(4, "mm"), con.cap = unit(1, "mm"), con.type = "straight"), plot$layers)
  plot$layers <- c(geom_mark_ellipse(aes(x = emb[, 1], y = emb[, 2], group = types, label = types, filter = types != "unknown"), 
                                     linetype = 2, size = 1, colour = "#555555", label.fontface = "plain",
                                  label.fontsize = 20, label.buffer = unit(2, "mm"), con.cap = unit(1, "mm"), con.type = "straight",
                                  radius = unit(2, "mm"), expand = -0.0025), plot$layers)
  if(arrow) {
    plot <- plot + geom_segment(aes(x = emb[cell, 1] + hor, xend = emb[cell, 1], 
                            y = emb[cell, 2] + vert, yend = emb[cell, 2]), size = 1,
                        colour = "red", arrow = arrow(length = unit(0.03, "npc")))
    
  }
  
  plot + theme(legend.text = element_text(size = 21), plot.title = element_text(size = 24))
}


plots[[1]] <- addLabels(plots[[1]], types13_A, um13_A)
plots[[2]] <- addLabels(plots[[2]], types13_B, um13_B, TRUE)
plots[[3]] <- addLabels(plots[[3]], types14, um14)

library(cowplot) #v0.9.4

figD <- plot_grid(plotlist = plots, nrow = 1, labels = c("A", "B", "C"), label_size = 24)

getGeneExpr <- function(emb, data, title) {
   lapply(c("Msx3", "Meis2", "Lhx5"), function(gene)
      ggplot() + geom_point(aes(x = emb[, 1], y = emb[, 2], colour = data[, gene])) + 
        labs(colour = gene, x = "", y = "", title = title) + 
        theme(legend.text = element_text(size = 15), plot.title = element_text(size = 24), 
              legend.title = element_text(size = 20), legend.key.height = unit(8, "mm")))
}

figDsup <- plot_grid(
  plot_grid(plotlist = getGeneExpr(um13_A, data13_A, "E13.5_A"), nrow = 3),
  plot_grid(plotlist = getGeneExpr(um13_B, data13_B, "E13.5_B"), nrow = 3),
  plot_grid(plotlist = getGeneExpr(um14, data14, "E14.5"), nrow = 3), ncol = 3)

png(filename = "figures/Fig_D.png", 
    width = 1500, height = 575)
print(figD)
dev.off()

png(filename = "figures/Fig_Dsup.png", 
    width = 1500, height = 1500)
print(figDsup)
dev.off()


sleepwalk(list(um13_A, um13_B, um14), list(comFeatures13_A, comFeatures13_B, comFeatures14), same = "features", pointSize = 2.5,
          titles = c("e13.5_A", "e13.5_B", "e14.5"), nrow = 1, saveToFile = "supplement/src/Fig_D.html")

file <- readLines("supplement/src/Fig_D.html")
line <- which(grepl("^set_up_chart\\(", file))
newFile <- c(file[1:(line - 1)], 
             "width = window.innerWidth/(n_charts) * 0.9;",
             file[line:length(file)])
writeLines(newFile, "supplement/src/Fig_D.html")