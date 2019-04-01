library(cowplot)
library(tidyverse)
library(ggrepel)
library(sleepwalk)

load(file ="cite_data.rda")

features <- citeSeq$pca[, 1:citeSeq$ndims]
vert <- (max(citeSeq$tsne[, 2]) - min(citeSeq$tsne[, 2]))/25
hor <- (max(citeSeq$tsne[, 1]) - min(citeSeq$tsne[, 1]))/25

getPlot <- function(cell, maxdist) {
  sleepwalk(citeSeq$tsne, features, maxdist)
  slw_snapshot(cell) + geom_segment(aes(x = citeSeq$tsne[cell, 1] + hor, xend = citeSeq$tsne[cell, 1], 
                                        y = citeSeq$tsne[cell, 2] - vert, yend = citeSeq$tsne[cell, 2]), size = 1,
                                    colour = "red", arrow = arrow(length = unit(0.03, "npc")))
}

plotA <- getPlot(7289, 9)
plotB <- getPlot(1659, 9)
plotC <- getPlot(8176, 12)
plotD <- getPlot(5618, 30)

figB <- plot_grid(plotA, plotB, plotC, plotD, labels = c("a", "b", "c", "d"))

png(filename = "../Fig_B.png", 
    width = 900, height = 1000)
print(figB)
dev.off()

sleepwalk(citeSeq$tsne, features, 9, saveToFile = "../supplement/Fig_B.html")
file <- readLines("../supplement/Fig_B.html")
line <- which(grepl("^set_up_chart\\(", file))
newFile <- c(file[1:(line - 1)], 
             "width = window.innerWidth/(n_charts + 1) * 0.9;",
             file[line],
             'd3.select("#embeddings")',
             '.select("tr")',
             '.append("td")',
             '.append("img")',
             '.attr("src", "Fig_A.png")',
             '.attr("width", width)',
             '.attr("height", width + 75);',
             file[(line + 1):length(file)])
writeLines(newFile, "../supplement/Fig_B.html")