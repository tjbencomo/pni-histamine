## Description: Make plots showing KITLG expression in tumor cells/fibroblasts
## and KIT expression in mast cells. Probably going in figure 2/S2

library(Seurat)
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)


data_dir <- file.path("data")
fig_dir <- file.path("figures", "manuscript")
pnifp <- file.path(data_dir, "pni_epithelial.rds")
sccfp <- file.path(data_dir, "lowrisk_epithelial.rds")
myeloidfp <- file.path(data_dir, "pni_myeloid.rds")
fibroblastfp <- file.path(data_dir, "integrated_fibroblasts.rds")
metadatafp <- file.path(data_dir, "cell_metadata.csv.gz")

keepTumorFineType <- c("KC_Basal", "KC_Cycling", "KC_TSK", "KC_PSK-2", "KC_Inflam")

met <- read_csv(metadatafp)
pniMet <- met %>% 
  filter(
    PNI == "PNI+", 
    broad_type == "Epithelial", 
    Poor_Quality == 0, 
    fine_type %in% keepTumorFineType
  )
sccMet <- met %>% 
  filter(
    condition2 == "PNI- Tumor", 
    broad_type == "Epithelial", 
    fine_type %in% keepTumorFineType, 
    Poor_Quality == 0
  )
myeloidMet <- met %>% filter(PNI == "PNI+", broad_type == "Myeloid", Poor_Quality == 0)
fibMet <- met %>% filter(condition == "Tumor", broad_type == "Fibroblast", Poor_Quality == 0)

pni <- readRDS(pnifp)
pni <- pni[, pniMet$full_barcode]
pni$fine_type <- pniMet$fine_type
pni$subpop <- case_when(
  pni$fine_type %in% c("KC_Inflam", "KC_PSK-2") ~ "KC_Basal",
  TRUE ~ pni$fine_type
) %>% str_replace("KC_", "")

scc <- readRDS(sccfp)
scc <- scc[, sccMet$full_barcode]
scc$fine_type <- sccMet$fine_type
scc$subpop <- case_when(
  scc$fine_type %in% c("KC_Inflam", "KC_PSK-2") ~ "KC_Basal",
  TRUE ~ scc$fine_type
) %>% str_replace("KC_", "")

myeloid <- readRDS(myeloidfp)
myeloid <- myeloid[, myeloidMet$full_barcode]
myeloid$subpop <- myeloidMet$fine_type

fibs <- readRDS(fibroblastfp)
fibs <- fibs[, fibMet$full_barcode]

## Draw PNI plots
tumorPlot <- VlnPlot(pni, "KITLG", group.by = "subpop", pt.size = 0)
myeloidPlot <- VlnPlot(myeloid, "KIT", group.by = "subpop", pt.size = 0)

## Draw comparative plots for KITLG
fibCompPlot <- VlnPlot(fibs, "KITLG", group.by = "PNI", pt.size = 0) +
  scale_fill_brewer(palette = "Accent", direction = -1)
tumors <- merge(scc, pni)
subpop_order <- c("Basal", "Cycling", "TSK")
tumors$subpop <- factor(tumors$subpop, levels = subpop_order)
tumorCompPlot <- VlnPlot(tumors, "KITLG", group.by = "subpop", split.by = "PNI", 
        pt.size = 0) + scale_fill_brewer(palette = "Accent", direction = -1)
# Get p-values for KITLG comparisons
FindMarkers(subset(tumors, subpop == "Basal"), features = c("KITLG"), logfc.threshold = 0,
            ident.1 = "PNI+", ident.2 = "PNI-", group.by = "PNI")
FindMarkers(subset(tumors, subpop == "Cycling"), features = c("KITLG"), logfc.threshold = 0,
            ident.1 = "PNI+", ident.2 = "PNI-", group.by = "PNI")
FindMarkers(subset(tumors, subpop == "TSK"), features = c("KITLG"), logfc.threshold = 0,
            ident.1 = "PNI+", ident.2 = "PNI-", group.by = "PNI")
FindMarkers(fibs, features = c("KITLG"), logfc.threshold = 0,
            ident.1 = "PNI+", ident.2 = "PNI-", group.by = "PNI")


## Save plots
ggsave(tumorPlot, filename = file.path(fig_dir, "Fig2_PNI_Tumor_KITLG_Subpop.eps"), width = 8, height = 6)
ggsave(myeloidPlot, filename = file.path(fig_dir, "Fig2_PNI_Myeloid_KIT.eps"), width = 8, height = 6)
ggsave(fibCompPlot, filename = file.path(fig_dir, "Fig2_Fibroblast_KITLG_Comparison.eps"), width = 6, height = 6)
ggsave(tumorCompPlot, filename = file.path(fig_dir, "Fig2_Tumor_KITLG_Comparison.eps"), width = 8, height = 6)

