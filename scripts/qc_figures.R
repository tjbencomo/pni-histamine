## Description: Plot quality control (QC) metrics for single cell samples


library(Seurat)
library(readr)
library(dplyr)
library(ggplot2)

data_dir <- file.path("data")
fig_dir <- file.path("figures", "manuscript")
cellsfp <- file.path(data_dir, "raw_cells.rds")
metadatafp <- file.path(data_dir, "cell_metadata.csv.gz")

met <- read_csv(metadatafp)
cells <- readRDS(cellsfp)

mitodf <- cells@meta.data %>%
  tibble::rownames_to_column(var = "full_barcode") %>%
  select(full_barcode, percent.mt)
rm(cells); gc()
met <- met %>%
  inner_join(mitodf) %>%
  mutate(sample_id = paste0(patient, "-", condition))

## Sample ID ordering
pid_order <- unique(met$sample_id)
pid_order <- rev(sort(pid_order))

## Group samples by condition (Ji-Normal/Ji-Tumor/Lee-PNI)
met <- met %>%
  mutate(group_label = case_when(
    condition == "Normal" ~ "Ji-Normal",
    condition == "Tumor" & PNI == "PNI-" ~ "Ji-Tumor",
    condition == "Tumor" & PNI == "PNI+" ~ "Lee-Tumor",
  ))
label_order <- c("Ji-Normal", "Ji-Tumor", "Lee-Tumor")

## Number of cells
## Already done with cell proportions

## Number of genes
genePlot <- met %>%
  mutate(group_label = factor(group_label, levels = c(label_order))) %>%
  # mutate(sample_id = factor(sample_id, levels = pid_order)) %>%
  ggplot(aes(nFeature_RNA, group_label)) +
  geom_boxplot(aes(fill = group_label)) +
  theme_classic() +
  guides(fill = "none") +
  labs(x = "Number of genes", y = "")
print(genePlot)


## Number of UMIs
umiPlot <- met %>%
  mutate(group_label = factor(group_label, levels = c(label_order))) %>%
  # mutate(sample_id = factor(sample_id, levels = pid_order)) %>%
  ggplot(aes(log10(nCount_RNA), group_label)) +
  geom_boxplot(aes(fill = group_label)) +
  theme_classic() +
  guides(fill = "none") +
  labs(x = "Log10 UMIs", y = "")
print(umiPlot)


## MT Reads
mitoPlot <- met %>%
  mutate(group_label = factor(group_label, levels = c(label_order))) %>%
  # mutate(sample_id = factor(sample_id, levels = pid_order)) %>%
  ggplot(aes(percent.mt, group_label)) +
  geom_boxplot(aes(fill = group_label)) +
  theme_classic() +
  guides(fill = "none") +
  labs(x = "% Mitochondrial UMIs", y = "")
print(mitoPlot)

combinedPlot <- genePlot | (umiPlot + theme(axis.text.y = element_blank())) | (mitoPlot + theme(axis.text.y = element_blank()))

## Save plots
ggsave(genePlot, filename = file.path(fig_dir, "FigS1_QC_Gene_Counts.eps"), width = 6, height = 9)
ggsave(umiPlot, filename = file.path(fig_dir, "FigS1_QC_UMI_Counts.eps"), width = 6, height = 9)
ggsave(mitoPlot, filename = file.path(fig_dir, "FigS1_QC_MT_Percent.eps"), width = 6, height = 9)
ggsave(combinedPlot, filename = file.path(fig_dir, "FigS1_QC_Plots.eps"), width = 15, height = 9)
