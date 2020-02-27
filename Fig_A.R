load(file = "data/cite_data.rda")

library(tidyverse)

citeSeq$tsne %>%
  as.data.frame() %>%
  rownames_to_column("cell") %>%
  mutate(cluster = citeSeq$clusters[cell]) %>%
  group_by(cluster) %>%
  summarise(tSNE_1 = mean(tSNE_1), tSNE_2 = mean(tSNE_2)) -> tsne_labels

citeSeq$tsne %>%
  as.data.frame() %>%
  rownames_to_column("cell") %>%
  mutate(cluster = citeSeq$clusters[cell]) %>%
  ggplot(aes(x = tSNE_1, y = tSNE_2)) + geom_point(aes(colour = cluster)) +
  geom_text(aes(label = cluster), size = 5, data = tsne_labels) +
  theme(legend.position = "bottom", axis.title = element_blank(),
        axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text = element_blank(), legend.title = element_blank(),
        panel.background = element_blank(),
        legend.text = element_text(size = 15)) + 
  scale_colour_hue(labels = c("CD4+ T Cells", "CD14+ Monocytes", "Natural Killer Cells", "Mouse Fibroblasts",
                              "B Cells", "CD8+ T Cells", "CD16+ Monocytes", "Unknown", "Stem Cells", "Megakaryocytes",
                              "Erythrocytes", "Dendritic cells", "Plasmacytoid Dendritic Cells")) +
  guides(colour = guide_legend(ncol = 3)) -> figA


png(filename = "figures/Fig_A.png", 
    width = 550, height = 550)
print(figA)
dev.off()

png(filename = "supplement/src/Fig_A.png", 
    width = 500, height = 575)
print(figA)
dev.off()