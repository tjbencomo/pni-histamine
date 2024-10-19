## Draw heatmaps showing phosphoproteomics data from control vs histamine vs anti-hstamine

library(readr)
library(dplyr)
library(readxl)
library(pheatmap)
library(stringr)
library(tidyr)
library(ggplot2)
library(patchwork)

fig_dir <- file.path("figures", "manuscript")
earlyfp <- file.path(fig_dir, "FigS6_Phospho_Early_Timepoint_Heatmap_Top.pdf")
latefp <- file.path(fig_dir, "FigS6_Phospho_Late_Timepoint_Heatmap.pdf")
sharedfp <- file.path(fig_dir, "Fig5_Phospho_Early_Timepoint_Comp_Proteins_Heatmap.pdf")
twoCondFp <- file.path(fig_dir, "Fig5_Phospho_Early_Timepoint_Ctrl_Hist_Proteins_Heatmap.pdf")
resdf <- read_csv("data/proteomics/histamine_anti_histamine_10min_120min_phospho_proteomics_results.csv")
proteindf <- read_excel(("data/proteomics/PTM_msstats_quant_results_wide.xlsx"),
                        na = "NA")

## Set colors for treatment conditions
# Makes ggplot boxplots match color order of pheatmap
treatmentColors <- c("#619CFF", "#F8766D", "#00BA38")
twoCondColors <- c("#619CFF", "#F8766D")


## Show all 3 conditions at 10min together
## Play with whether to include all padj < .05 genes or just N top hits by fold change
## Using only top N genes shows better pattern changes between conditions 
earlyTimePointRes <- resdf %>%
  filter(Label %in% c("anti_10min_vs_his_10min", "his_10min_vs_ctrl_10min")) %>%
  group_by(Label) %>%
  filter(adj.pvalue < .05) %>%
  ungroup()
earlySamples <- c(paste0("ctrl_10min_0", 1:3), paste0("his_10min_0", 1:3), paste0("anti_10min_0", 1:3))
earlySampleInfo <- data.frame(
  sampleID = earlySamples
)
rownames(earlySampleInfo) <- earlySampleInfo$sampleID
earlySampleInfo <- earlySampleInfo %>%
  mutate(Treatment = case_when(
    str_detect(sampleID, "ctrl") ~ "Vehicle",
    str_detect(sampleID, "his") ~ "Histamine",
    str_detect(sampleID, "anti") ~ "Anti-Histamine"
  ))



## Shows genes that show complementary DE changes in both ctrl vs hist and hist vs anti-hist
## Only doing for early timepoint
histUp <- resdf %>%
  filter(Label == "his_10min_vs_ctrl_10min", log2FC > 0, adj.pvalue < .05)
histDn <- resdf %>%
  filter(Label == "his_10min_vs_ctrl_10min", log2FC < 0, adj.pvalue < .05)
antiUp <- resdf %>%
  filter(Label == "anti_10min_vs_his_10min", log2FC > 0, adj.pvalue < .05)
antiDn <- resdf %>%
  filter(Label == "anti_10min_vs_his_10min", log2FC < 0, adj.pvalue < .05)
sharedProteins <- c(
  intersect(histUp %>% pull(Protein), antiDn %>% pull(Protein)),
  intersect(histDn %>% pull(Protein), antiUp %>% pull(Protein))
)


# Save genesets for Enrichr analysis
histUpAntiDnUnion <- union(histUp$gene, antiDn$gene)
histUpAntiDnIntersect <- intersect(histUp$gene, antiDn$gene)
histDnAntiUpUnion <- union(histDn$gene, antiUp$gene)
histDnAntiUpIntersect <- intersect(histDn$gene, antiUp$gene)

write(histUpAntiDnUnion, "data/proteomics/hist_up_anti_hist_down_10min_union.txt")
write(histUpAntiDnIntersect, "data/proteomics/hist_up_anti_hist_down_10min_intersect.txt")
write(histDnAntiUpUnion, "data/proteomics/hist_down_anti_hist_up_10min_union.txt")
write(histDnAntiUpIntersect, "data/proteomics/hist_down_anti_hist_up_10min_intersect.txt")


## Heatmap for control + histamine groups
twoCondProteins <- c(histUp$Protein, histDn$Protein)
twoCondSamples <- earlySamples[!str_detect(earlySamples, "anti")]
twoCondInfo <- earlySampleInfo %>% filter(Treatment != "Anti-Histamine")
twoCondMat <- proteindf %>%
  filter(ProteinName %in% twoCondProteins) %>%
  select(all_of(c("ProteinName", twoCondSamples))) %>%
  data.frame() %>%
  tibble::column_to_rownames("ProteinName") %>%
  as.matrix()
twoCondPlot <- pheatmap(
  twoCondMat, 
  show_rownames = F, 
  scale = "row", 
  cluster_cols = F, 
  annotation_col = twoCondInfo["Treatment"],
  border_color = NA,
  treeheight_row = 0,
  annotation_colors = list(Treatment=c(Vehicle="#bbafd1", Histamine="#90c686"))
)


# 
# pdf(twoCondFp)
# twoCondPlot
# dev.off()





