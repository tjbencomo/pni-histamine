## Description: Plot marker gene scores for borad UMAP broad cell types across all 3 conditions

library(Seurat)
library(readr)
library(dplyr)
library(ggplot2)
library(scCustomize)

data_dir <- "data"
cellsfp <- file.path(data_dir, "cells.rds")
metadatafp <- file.path(data_dir, "broad_labels.csv")
fullMetfp <- file.path(data_dir, "cell_metadata.csv.gz")
outfp <- file.path("figures", "manuscript", "Fig1_Broad_UMAP_Marker_Genes.eps")

cells <- readRDS(cellsfp)
met <- read_csv(metadatafp)
cellMet <- read_csv(fullMetfp)

new_cluster_labels <- met$cell_type
names(new_cluster_labels) <- met$cluster_id
cells <- RenameIdents(cells, new_cluster_labels)

cellMet <- cellMet %>%
  arrange(match(full_barcode, colnames(cells)))
stopifnot(all(colnames(cells) == cellMet$full_barcode))
cells$broad <- cellMet$broad_type
cells$fine <- cellMet$fine_type

# keepCells <- cellMet %>%
#   filter(Poor_Quality == 0) %>%
#   group_by(broad_type) %>%
#   slice_head(prop = .) %>%
#   ungroup()

# cells <- cells[, keepCells$full_barcode]
# table(keepCells$broad_type)
# table(cellMet$broad_type)

doubletIds <- cellMet %>% filter(Poor_Quality == 1) %>% pull(full_barcode)
doubletIdx <- which(colnames(cells) %in% doubletIds)
cells <- cells[, -doubletIdx]


# doubletIds <- cellMet %>% filter(Poor_Quality == 1) %>% pull(full_barcode)
# doubletIdx <- which(colnames(cells) %in% doubletIds)
# cells <- cells[, -doubletIdx]
gc()
print("Finished subsetting")
# umap_plot <- DimPlot(cells)
markers <- c("KRT5", "CD3E", "LYZ", "VWF", "LUM", "PMEL", "TPSAB1", "HTN3", "CD79A")
marker_plot <- FeaturePlot_scCustom(cells, features = markers)
# print(umap_plot)
# print(marker_plot)
ggsave(filename = outfp, plot = marker_plot, width = 10, height = 6)

# Make individual panels so Ankit can manipulate them in CorrelDraw
for (i in 1:length(markers)) {
  gene <- markers[i]
  print(paste("Saving", gene, "graph"))
  p <- FeaturePlot_scCustom(cells, features = gene)
  fn <- paste0("Fig1_Broad_UMAP_Marker_Genes_", gene, ".eps")
  fp <- file.path("figures", "manuscript", fn)
  ggsave(filename = fp, plot = p, width = 6, height = 5)
}
