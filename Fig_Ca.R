library(sleepwalk)

load("data/cite_data.rda")

sleepwalk(list(citeSeq$tsne, citeSeq$umap), citeSeq$pca[, 1:citeSeq$ndims], saveToFile = "supplementary/src/Fig_Ca.html", 
          titles = c("t-SNE", "UMAP"))

file <- readLines("supplementary/src/Fig_Ca.html")
line <- which(grepl("^set_up_chart\\(", file))
newFile <- c(file[1:(line - 1)], 
             "width = window.innerWidth/(n_charts) * 0.9;",
             file[line:length(file)])
writeLines(newFile, "supplement/src/Fig_Ca.html")