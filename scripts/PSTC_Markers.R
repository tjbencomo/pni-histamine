## Description: Identify marker genes specific to KC PSTC1/2.
## Make expression plots and marker gene table files


library(Seurat)
library(readr)
library(dplyr)
library(ggplot2)

cells <- readRDS("data/cells.rds")
met <- read_csv("data/cell_metadata.csv.gz")

tumorMet <- met %>%
  filter(condition2 == "PNI+ Tumor", level1_celltype == "Epithelial", Poor_Quality == 0)

cells <- cells[, tumorMet$full_barcode]
stopifnot(all(colnames(cells) == tumorMet$full_barcode))
cells
Idents(cells) <- tumorMet$level2_celltype
markers <- FindAllMarkers(cells, only.pos = T, min.pct = 0.25)
topMarkers <- markers %>%
  filter(p_val_adj < 1e-4) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 10)

topMarkersShort <- markers %>%
  filter(p_val_adj < 1e-4) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5)

topMarkers %>%
  filter(cluster == "KC_PSTC1")

topMarkers %>%
  filter(cluster == "KC_PSTC2")

DoHeatmap(cells, features = topMarkersShort$gene)

tsk_genes <- c("MMP10", "PTHLH", "FEZ1", "IL24", "KCNMA1", "INHBA", "MAGEA4", "NT5E", "LAMC2", "SLITRK6")
basal_genes <- c("KRT15", "CCL2", "COL17A1", "CXCL14", "DST", "CCN1", "FTH1", "MT2A", "IGFBP5", "THBS2")
diff_genes <-c("KRTDAP", "KRT10", "KRT1", "S100A7", "SBSN", "DMKN", "KRT60", "LYPD3", "KRT6A", "CALML5")
cycling_genes <- c("STMN1", "HIST1H4C", "TUBA1B", "PTTG1", "HMGB2", "H2AFZ", "TOP2A", "UBE2C")

VlnPlot(cells, basal_genes, pt.size = 0)
VlnPlot(cells, cycling_genes, pt.size = 0)
VlnPlot(cells, diff_genes, pt.size = 0)
VlnPlot(cells, tsk_genes, pt.size = 0)

VlnPlot(cells, c("IL20", "CCL2", "ZBTB16", "SEMA5A", "NFIL3", "DUSP2", "HES1", "EGR2", "MAFF", "SOCS3"), pt.size = 0)
VlnPlot(cells, c("CLCA4", "MMP7", "SOSTDC1", "EN1", "NES", "PRICKLE1", "BBOX1", "CYP1B1", "AGR2", "RHOV"), pt.size = 0)

plot_genes <- c("CXCL14", "COL17A1", "STMN1", "TOP2A", "KRTDAP", "KRT10", "PTHLH", "INHBA", "MAFF", "SOCS3", "SOSTDC1", "BBOX1")
markerPlot <- VlnPlot(cells, plot_genes, pt.size = 0)
markerPlot

ggsave(
  plot = markerPlot,
  filename = "figures/manuscript/Tumor_Subpop_Marker_Genes.svg",
  width = 10,
  height = 10
)
write_csv(topMarkers, "data/PNI_Tumor_Subpops_Marker_Genes.csv")
write_csv(markers, "data/PNI_Tumor_All_Subpops_Marker_Genes.csv")
# 
# 
# cells <- ScaleData(cells, features = rownames(cells))
# DoHeatmap(cells, features = plot_genes)
# 
# 
# stress <- scan("data/puram_stress.txt", what = character())
# 
# cells <- AddModuleScore(cells, features = list(stress), name = "Stress")
# VlnPlot(cells, "Stress1", pt.size=0)
