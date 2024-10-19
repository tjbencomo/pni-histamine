## Description: Draw composition plots for each patient showing the number
## of cells coming from each broad cell type


library(readr)
library(dplyr)
library(ggplot2)
library(patchwork)

data_dir <- "data"
fig_dir <- file.path("figures", "manuscript")
infp <- file.path(data_dir, "cell_metadata.csv.gz")
met <- read_csv(infp)
met$level1_celltype <- factor(met$level1_celltype)
met$patient <- factor(met$patient)
met2 <- met %>% filter(Poor_Quality == 0) %>%
  mutate(group_label = case_when(
    condition == "Normal" ~ "Ji-Normal",
    condition == "Tumor" & PNI == "PNI-" ~ "Ji-Tumor",
    condition == "Tumor" & PNI == "PNI+" ~ "Lee-Tumor",
  ))
label_order <- c("Ji-Normal", "Ji-Tumor", "Lee-Tumor")

cellCounts <- met2 %>%
  count(group_label, level1_celltype, .drop = F, name = "cells")

totalCounts <- met2 %>%
  count(group_label, name = "total_cells")

props <- cellCounts %>%
  inner_join(totalCounts) %>%
  mutate(prop = cells / total_cells)

widthSize <- .5
# pid_order <- unique(paste0(props$patient, "-", props$condition))
# pid_order <- rev(sort(pid_order))
propPlot <- props %>%
  mutate(group_label = factor(group_label, levels = c(label_order))) %>%
  # mutate(pid = factor(paste0(patient, "-", condition), levels = pid_order)) %>%
  ggplot(aes(prop, group_label)) +
  geom_col(position = 'fill', aes(fill = level1_celltype), width = widthSize) +
  theme_classic() +
  labs(x = "Proportion", y = "", fill = "") +
  scale_fill_brewer(palette = "Set2")
print(propPlot)

totalPlot <- totalCounts %>%
  mutate(group_label = factor(group_label, levels = c(label_order))) %>%
  # mutate(pid = factor(paste0(patient, "-", condition), levels = pid_order)) %>%
  ggplot(aes(total_cells, group_label)) +
  geom_col(fill = "white", color = "black", width = widthSize) +
  theme_classic() +
  labs(x = "Number of cells", y = "")
print(totalPlot)

combinedPlot <- propPlot + (totalPlot + theme(axis.text.y = element_blank())) + 
  plot_layout(guides = 'collect')
print(combinedPlot)


ggsave(propPlot, filename = file.path(fig_dir, "FigS1_Patient_Composition.eps"), width = 6, height = 8)
ggsave(totalPlot, filename = file.path(fig_dir, "FigS1_Patient_Total_Cells.eps"), width = 6, height = 8)
ggsave(combinedPlot, filename = file.path(fig_dir, "FigS1_Patient_Composition_Totals_Combined.eps"), width = 8, height = 10)

