## Description: Plot UMAP showing broad cell types across all 3 conditions

library(Seurat)
library(readr)
library(dplyr)
library(ggplot2)

data_dir <- "data"
cellsfp <- file.path(data_dir, "cells.rds")
metadatafp <- file.path(data_dir, "broad_labels.csv")
fullMetfp <- file.path(data_dir, "cell_metadata.csv.gz")
outfp <- file.path("figures", "manuscript", "Fig1_Broad_UMAP.eps")

cells <- readRDS(cellsfp)
met <- read_csv(metadatafp)
cellMet <- read_csv(fullMetfp)

new_cluster_labels <- met$cell_type
names(new_cluster_labels) <- met$cluster_id
cells <- RenameIdents(cells, new_cluster_labels)

doubletIds <- cellMet %>% filter(Poor_Quality == 1) %>% pull(full_barcode)
doubletIdx <- which(colnames(cells) %in% doubletIds)
cells <- cells[, -doubletIdx]

umap_plot <- DimPlot(cells)
print(umap_plot)
ggsave(filename = outfp, plot = umap_plot, width = 7, height = 6)
