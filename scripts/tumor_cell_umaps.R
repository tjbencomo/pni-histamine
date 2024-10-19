## Description: Plot UMAPs showing lowrisk and PNI tumor cells.

library(Seurat)
library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(harmony)

integration_pipeline <- function(seu) {
  seu <- seu %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData(vars.to.regress = "percent.mt") %>%
    RunPCA(npcs = 30) %>%
    RunHarmony("patient", max.iter.harmony = 30) %>%
    RunUMAP(reduction = "harmony", dims = 1:30) %>%
    FindNeighbors(reduction = "harmony")
}

data_dir <- "data"
fig_dir <- file.path("figures", "manuscript")
metadatafp <- file.path(data_dir, "cell_metadata.csv.gz")

met <- read_csv(metadatafp)
pni <- readRDS("data/pni_epithelial.rds")
lr <- readRDS("data/lowrisk_epithelial.rds")

pniMet <- met %>%
  filter(condition2 == "PNI+ Tumor", level1_celltype == "Epithelial", Poor_Quality == 0)
cell_levels <- paste0("KC_", c("Basal", "Cycling", "Diff", "TSK", "PSTC-1", "PSTC-2"))
pniMet$level2_celltype <- factor(pniMet$level2_celltype, levels = cell_levels)
lrMet <- met %>%
  filter(condition2 == "PNI- Tumor", level1_celltype == "Epithelial", Poor_Quality == 0) %>%
  filter(level2_celltype != "Pilo/Eccrine") %>%
  mutate(level2_celltype = factor(level2_celltype, levels = cell_levels))

pni <- pni[, pniMet$full_barcode]
stopifnot(all(colnames(pni) == pniMet$full_barcode))
pni$level2_celltype <- pniMet$level2_celltype
lr <- lr[, lrMet$full_barcode]
stopifnot(all(colnames(lr) == lrMet$full_barcode))
lr$level2_celltype <- lrMet$level2_celltype

# Here we re-integrate to remove outlier cells that were removed when we dropped clusters/low QC cells
pni <- integration_pipeline(pni)
pniPlot <- DimPlot(pni, group.by = "level2_celltype") + scale_color_brewer(palette = "Set2") + ggtitle("")
pniPlot
ggsave(pniPlot, filename = file.path(fig_dir, "Fig1_PNI_Epithelial_UMAP.eps"), width = 6, height = 5)

lr <- integration_pipeline(lr)
lrPlot <- DimPlot(lr, group.by = "fine_type") + scale_color_brewer(palette = "Set2")  + ggtitle("")
ggsave(lrPlot, filename = file.path(fig_dir, "FigS2_LR_Epithelial_UMAP.eps"), width = 6, height = 5)
print(lrPlot + pniPlot)
