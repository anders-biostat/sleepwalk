load("data/cite_data.rda")

library(sleepwalk)
sleepwalk(list(citeSeq$tsne, citeSeq$tsne), list(citeSeq$adt, citeSeq$pca[, 1:13]),
          compare = "distances", pointSize = 2, titles = c("Epitope distances", "mRNA distances"),
          saveToFile = "supplement/src/Fig_E.html")
file <- readLines("supplement/src/Fig_E.html")
line <- which(grepl("^set_up_chart\\(", file))
newFile <- c(file[1:(line - 1)], 
             "width = window.innerWidth/(n_charts) * 0.9;",
             file[line:length(file)])
writeLines(newFile, "supplement/src/Fig_E.html")