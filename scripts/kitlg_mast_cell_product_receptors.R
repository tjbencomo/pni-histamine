## Description: Make plots showing expression of KITLG, and receptors for mast cell degranulation
## products (HRH1 etc). Also make plot for KIT expression in myeloid cells.

library(Seurat)
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(muscat)



fix_colors <- function(plots) {
    plots <- lapply(
        X = plots,
        FUN = function(p) p + scale_fill_brewer(palette = "Accent", direction = -1) + labs(x = "")
    )
    CombinePlots(plots = plots, legend = 'right')
}


data_dir <- file.path("data")
fig_dir <- file.path("figures", "manuscript")
pnifp <- file.path(data_dir, "pni_epithelial.rds")
sccfp <- file.path(data_dir, "lowrisk_epithelial.rds")
myeloidfp <- file.path(data_dir, "pni_myeloid.rds")
fibroblastfp <- file.path(data_dir, "integrated_fibroblasts.rds")
metadatafp <- file.path(data_dir, "cell_metadata.csv.gz")

keepTumorFineType <- c("KC_Basal", "KC_Cycling", "KC_TSK")

met <- read_csv(metadatafp)
pniMet <- met %>% 
  filter(
    PNI == "PNI+", 
    level1_celltype == "Epithelial", 
    Poor_Quality == 0, 
    level2_celltype %in% keepTumorFineType
  )
sccMet <- met %>% 
  filter(
    condition2 == "PNI- Tumor", 
    level1_celltype == "Epithelial", 
    level2_celltype %in% keepTumorFineType, 
    Poor_Quality == 0
  )
myeloidMet <- met %>% filter(PNI == "PNI+", level1_celltype == "Myeloid", Poor_Quality == 0)
fibMet <- met %>% filter(condition == "Tumor", level1_celltype == "Fibroblast", Poor_Quality == 0)

pni <- readRDS(pnifp)
pni <- pni[, pniMet$full_barcode]
stopifnot(all(colnames(pni) == pniMet$full_barcode))
pni$subpop <- str_replace(pniMet$level2_celltype, "KC_", "")

scc <- readRDS(sccfp)
scc <- scc[, sccMet$full_barcode]
stopifnot(all(colnames(scc) == sccMet$full_barcode))
scc$subpop <- str_replace(sccMet$level2_celltype, "KC_", "")

tumors <- merge(scc, pni)
subpop_order <- c("Basal", "Cycling", "TSK")
tumors$subpop <- factor(tumors$subpop, levels = subpop_order)


myeloid <- readRDS(myeloidfp)
myeloid <- myeloid[, myeloidMet$full_barcode]
stopifnot(all(colnames(myeloid) == myeloidMet$full_barcode))
myeloid$subpop <- myeloidMet$level2_celltype

fibs <- readRDS(fibroblastfp)
fibs <- fibs[, fibMet$full_barcode]
stopifnot(all(colnames(fibs) == fibMet$full_barcode))

# Mast cell degranulation product receptors
histamineGenes <- c("HRH1", "HRH2", "HRH3", "HRH4")
tnfGenes <- c("TNFRSF1A", "TNFRSF1B")
tryptaseGenes <- c("F2RL1")
interleukGenes <- c("IL4R", "IL5RA", "CSF2RB", "IL6R", "IL13RA1", "IL17RA", 
                    "IL17RB", "IL17RC", "IL17RD", "IL17RE")
vegfGenes <- c("FLT1", "KDR", "FLT4")
thyroidPaperGenes <- c("ACKR1", "CXCR1", "CXCR2", "CXCR3")
# Show these in supplemental figure
suppMastGenes <- c("IL17RA", "FLT1", "F2RL1")

genes <- unique(c("KITLG", histamineGenes, tnfGenes, tryptaseGenes, 
                  interleukGenes, vegfGenes, thyroidPaperGenes, suppMastGenes))

## Draw PNI Plot for Myeloid KIT Expression
myeloidKitPlot <-  VlnPlot(myeloid, "KIT", group.by = "subpop", pt.size = 0)
myeloidKitPlot

## Draw comparative plots for KITLG
## Fibroblast KITLG
fibCompPlot <- VlnPlot(fibs, "KITLG", group.by = "PNI", pt.size = 0) +
  scale_fill_brewer(palette = "Accent", direction = -1)
fibCompPlot

## Tumor Subpops KITLG
tumorCompPlot <- VlnPlot(tumors, "KITLG", group.by = "subpop", split.by = "PNI", 
        pt.size = 0) + scale_fill_brewer(palette = "Accent", direction = -1)
tumorCompPlot

## Tumor Subpop Mast Cell Product Receptors
histaminePlot <- VlnPlot(tumors, histamineGenes, group.by = "subpop", split.by = "PNI", pt.size = 0, combine = F)
suppMediatorPlot <- VlnPlot(tumors, suppMastGenes, group.by = "subpop", split.by = "PNI", pt.size = 0, combine = F)
histaminePlot <- fix_colors(histaminePlot)
suppMediatorPlot <- fix_colors(suppMediatorPlot)
histaminePlot
suppMediatorPlot

# Get p-values for gene comparisons comparisons
FindMarkers(subset(tumors, subpop == "Basal"), features = genes, logfc.threshold = 0,
            ident.1 = "PNI+", ident.2 = "PNI-", group.by = "PNI")
FindMarkers(subset(tumors, subpop == "Cycling"), features = genes, logfc.threshold = 0,
            ident.1 = "PNI+", ident.2 = "PNI-", group.by = "PNI")
FindMarkers(subset(tumors, subpop == "TSK"), features = genes, logfc.threshold = 0,
            ident.1 = "PNI+", ident.2 = "PNI-", group.by = "PNI")
FindMarkers(fibs, features = genes, logfc.threshold = 0,
            ident.1 = "PNI+", ident.2 = "PNI-", group.by = "PNI")


## Save plots
ggsave(tumorPlot, filename = file.path(fig_dir, "Fig2_PNI_Tumor_KITLG_Subpop.eps"), width = 8, height = 6)
ggsave(myeloidPlot, filename = file.path(fig_dir, "Fig2_PNI_Myeloid_KIT.eps"), width = 8, height = 6)
ggsave(fibCompPlot, filename = file.path(fig_dir, "Fig2_Fibroblast_KITLG_Comparison.eps"), width = 6, height = 6)
ggsave(tumorCompPlot, filename = file.path(fig_dir, "Fig2_Tumor_KITLG_Comparison.eps"), width = 8, height = 6)

ggsave(histaminePlot, filename = file.path(fig_dir, "Fig3_Tumor_HRH_Expression.eps"), width = 8, height = 4)
ggsave(suppMediatorPlot, filename = file.path(fig_dir, "FigS4_Tumor_Mast_Mediator_Receptor_Expression.pdf"), width = 8, height = 6)

## Pseudobulk tests for HRH1
tumors$PNI_Label <- ifelse(tumors$PNI == "PNI+", "PNI_POS", "PNI_NEG")
sce <- as.SingleCellExperiment(tumors)
sce <- prepSCE(sce, 
               kid = "subpop", # subpopulation assignments
               gid = "PNI_Label",  # group IDs (ctrl/stim)
               sid = "patient",   # sample IDs (ctrl/stim.1234)
               drop = TRUE)
pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))
res <- pbDS(pb)
res$table$PNI_POS$Basal %>%
    filter(gene %in% genes) %>%
    arrange(p_val)
res$table$PNI_POS$Cycling %>%
    filter(gene %in% genes) %>%
    arrange(p_val)
res$table$PNI_POS$TSK %>%
    filter(gene %in% genes) %>%
    arrange(p_val)
