set.seed(1234)

library(Seurat)
library(cowplot)
library(tidyverse)

e135 <- read_rds("data/e13_A.rds")

e135 <- NormalizeData(e135)
e135 <- FindVariableGenes(e135, do.plot = FALSE, y.cutoff = 0.5)

#save normalized expression of only variable genes
data <- t(as.matrix(e135@data[e135@var.genes, ]))

library(irlba)
pca <- prcomp_irlba(data, n = 50)

library(Rtsne)
tsne <- Rtsne(pca$x, pca = F)

library(uwot)
um <- umap(pca$x, spread = 7)

write_rds(um, "data/e13A_umap.rds", compress = "gz")

#used to draw arrows
vert_tsne <- (max(tsne$Y[, 2]) - min(tsne$Y[, 2]))/25
hor_tsne <- (max(tsne$Y[, 1]) - min(tsne$Y[, 1]))/25
vert_um <- (max(um[, 2]) - min(um[, 2]))/25
hor_um <- (max(um[, 1]) - min(um[, 1]))/25

library(sleepwalk)
sleepwalk(list(tsne$Y, um), pca$x, pointSize = 2.5, titles = c("t-SNE", "UMAP"))

cell <- 1914
snsh <- slw_snapshot(cell, 1, TRUE)
figC <- plot_grid(snsh[[1]] + geom_segment(aes(x = tsne$Y[cell, 1] + hor_tsne, xend = tsne$Y[cell, 1], 
                                       y = tsne$Y[cell, 2] + vert_tsne, yend = tsne$Y[cell, 2]), size = 1,
                                   colour = "red", arrow = arrow(length = unit(0.03, "npc"))),
          snsh[[2]] + geom_segment(aes(x = um[cell, 1] - hor_um, xend = um[cell, 1], 
                                       y = um[cell, 2] - vert_um, yend = um[cell, 2]), size = 1,
                                   colour = "red", arrow = arrow(length = unit(0.03, "npc"))),
          labels = c("a", "b"))

png(filename = "figures/Fig_C.png", 
    width = 1000, height = 575)
print(figC)
dev.off()

sleepwalk(list(tsne$Y, um), pca$x, pointSize = 2.5, saveToFile = "supplement/src/Fig_C.html", titles = c("t-SNE", "UMAP"))
file <- readLines("supplement/src/Fig_C.html")
line <- which(grepl("^set_up_chart\\(", file))
newFile <- c(file[1:(line - 1)], 
             "width = window.innerWidth/(n_charts) * 0.9;",
             file[line:length(file)])
writeLines(newFile, "supplement/src/Fig_C.html")
