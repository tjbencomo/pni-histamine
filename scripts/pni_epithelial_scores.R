## Description: Plot scores for Basal/Cycling/Diff/TSK/Stress/Immune evasion
## signatures on PNI UMAP. Immune evasion comes from googling the
## Don't eat me signature and represents immune suppression/evasion.
## Likely going in figure 1/S1

library(Seurat)
library(readr)
library(dplyr)
library(harmony)
library(stringr)
library(Nebulosa)
library(scCustomize)
library(patchwork)

integration_pipeline <- function(seu) {
  seu <- seu %>%
    # NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData(vars.to.regress = "percent.mt") %>%
    RunPCA(npcs = 30) %>%
    RunHarmony("patient", max.iter.harmony = 30) %>%
    RunUMAP(reduction = "harmony", dims = 1:30) %>%
    FindNeighbors(reduction = "harmony")
}

fig_dir <- file.path("figures", "manuscript")
cells <- readRDS("data/pni_epithelial.rds")
met <- read_csv("data/cell_metadata.csv.gz")
pniMet <- met %>% filter(condition2 == "PNI+ Tumor", broad_type == "Epithelial", Poor_Quality == 0)
cells <- cells[, pniMet$full_barcode]

DimPlot(cells)
cells <- integration_pipeline(cells)
DimPlot(cells)
pniMet <- pniMet %>%
  mutate(subpop = case_when(
    fine_type == "KC_PSK-2" ~ "KC_Basal",
    fine_type == "KC_Inflam" ~ "KC_Basal",
    TRUE ~ fine_type
  )) %>%
  mutate(subpop = str_replace(subpop, "KC_", ""))
cells$subpop <- pniMet$subpop

DimPlot(cells, group.by = "subpop")

## Score with Basal/TSK/Diff/Stress/
basalGenes <- c("KRT15", "CCL2", "COL17A1", "CXCL14", "DST", "CCN1", "FTH1", "MT2A", "IGFBP5", "THBS2")
tskGenes <- c("MMP10", "MMP1", "PTHLH", "FEZ1", "IL24", "KCNMA1", "INHBA", "MAGEA4", "NT5E", "LAMC2", "SLITRK6")
diffGenes <- c("KRTDAP", "KRT10", "KRT1", "S100A7", "SBSN", "DMKN", "KRT60", "LYPD3", "KRT6A", "CALML5")
stressGenes <- scan("data/puram_stress.txt", what = character())
immuneGenes <- c("CD274", "CD47", "CD24")
cyclingGenes <- c("STMN1", "HIST1H4C", "TUBA1B", "PTTG1", "HMGB2", "H2AFZ", "TOP2A", "UBE2C")
genesets <- list("Basal" = basalGenes, "Cycling" = cyclingGenes, "Diff" = diffGenes, 
                 "TSK" = tskGenes, "Stress" = stressGenes, "Immune" = immuneGenes)
cells <- AddModuleScore(cells, features = genesets, name = "Score")
idx <- which(str_detect(colnames(cells@meta.data), "Score[0-9]"))
colnames(cells@meta.data)[idx] <- names(genesets)

# FeaturePlot(cells, c("Basal", "Cycling", "Diff", "TSK", "Stress", "Immune"))
# plot_density(cells, c("Basal", "Cycling", "Diff", "TSK", "Stress", "Immune"), reduction = "umap")
# 
# FeaturePlot_scCustom(seurat_object = cells, features = "TSK")
# FeaturePlot(cells, features = "TSK")
# 
# FeaturePlot_scCustom(seurat_object = cells, 
#                      features = c("Basal", "Cycling", "Diff", "TSK", "Stress", "Immune"))

## Plot signature scores
basalPlot <- FeaturePlot_scCustom(cells, features = "Basal", na_cutoff = NULL)
cyclingPlot <- FeaturePlot_scCustom(cells, features = "Cycling", na_cutoff = NULL)
diffPlot <- FeaturePlot_scCustom(cells, features = "Diff", na_cutoff = NULL)
tskPlot <- FeaturePlot_scCustom(cells, features = "TSK", na_cutoff = NULL)
stressPlot <- FeaturePlot_scCustom(cells, features = "Stress", na_cutoff = NULL)
# immunePlot <- FeaturePlot_scCustom(cells, features = "Immune", na_cutoff = NULL)
# pdl1Plot <- FeaturePlot_scCustom(cells, features = "CD274")
# ctla4Plot <- FeaturePlot_scCustom(cells, features = "CTLA4")
# cd47Plot <- FeaturePlot_scCustom(cells, features = "CD47")
# cd27Plot <- FeaturePlot_scCustom(cells, features = "CD27")
# ki67Plot <- FeaturePlot_scCustom(cells, features = "MKI67")

fig1E <- tskPlot | stressPlot
figS2D <- basalPlot | cyclingPlot | diffPlot

## Save plots
ggsave(basalPlot, filename = file.path(fig_dir, "Fig1_PNI_Tumor_Basal_Score.eps"), width = 8, height = 6)
ggsave(cyclingPlot, filename = file.path(fig_dir, "Fig1_PNI_Tumor_Cycling_Score.eps"), width = 8, height = 6)
ggsave(diffPlot, filename = file.path(fig_dir, "Fig1_PNI_Tumor_Diff_Score.eps"), width = 6, height = 5)
ggsave(tskPlot, filename = file.path(fig_dir, "Fig1_PNI_Tumor_TSK_Score.eps"), width = 6, height = 5)
ggsave(stressPlot, filename = file.path(fig_dir, "Fig1_PNI_Tumor_Stress_Score.eps"), width = 6, height = 5)
ggsave(immunePlot, filename = file.path(fig_dir, "Fig1_PNI_Tumor_Immune_Evasion_Score.eps"), width = 8, height = 6)
ggsave(pdl1Plot, filename = file.path(fig_dir, "Fig1_PNI_Tumor_PDL1.eps"), width = 8, height = 6)
ggsave(ctla4Plot, filename = file.path(fig_dir, "Fig1_PNI_Tumor_CTLA4.eps"), width = 8, height = 6)
ggsave(ki67Plot, filename = file.path(fig_dir, "Fig1_PNI_Tumor_KI67.eps"), width = 8, height = 6)
ggsave(fig1E, filename = file.path(fig_dir, "Fig1_PNI_Epithelial_Subpopulation_Scores.eps"), width = 8, height = 4)
ggsave(figS2D, filename = file.path(fig_dir, "FigS2_PNI_Epithelial_Subpopulation_Supplemental_Scores.eps"), width = 12, height = 4)

