## Description: Score tumor cells using DvP signature from Bailey 2023

library(Seurat)
library(readr)
library(readxl)
library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)



met <- read_csv("data/cell_metadata.csv.gz")
# tumorMet <- met %>% filter(condition == "Tumor", broad_type == "Epithelial", Poor_Quality == 0)
cellMet <- met %>% filter(broad_type == "Epithelial", Poor_Quality == 0)
cells <- readRDS("data/cells.rds")
cells <- cells[, cellMet$full_barcode]
# stopifnot(all(colnames(cells) == cellMet$full_barcode))
# cells$fine_type <- cellMet$fine_type
# cells$condition2 <- cellMet$condition2


dvp_df <- read_excel("data/Bailey_2023_Signatures.xlsx", sheet = "DvP_Signature")
dvp_df <- dvp_df %>% filter(Gene_symbol %in% rownames(cells))

dvpMat <- cells@assays$RNA@data[dvp_df$Gene_symbol, ]
dvpScores <-  colSums(dvpMat * dvp_df$Coefficient)
stopifnot(all(names(dvpScores) == colnames(cells)))
cells$DvP_Score <- scale(dvpScores)[, 1]

cellMet <- cellMet %>%
  inner_join(
    cells@meta.data %>%
      tibble::rownames_to_column("full_barcode") %>%
      select(full_barcode, DvP_Score)
  )


## Check PNI- vs PNI+ single cell level
wilcox.test(DvP_Score ~ PNI, data=cellMet)
t.test(DvP_Score ~ PNI, data=cellMet)
## Check PNI- vs PNI+ mixed effects level
f <- lmer(DvP_Score ~ condition2 + (1 | patient), data=cellMet)
summary(f)

## Set PNI+ tumor as reference level to get p-values for each comparison
cellMet <- cellMet %>%
  mutate(group = factor(condition2, levels = c("PNI+ Tumor", "PNI- Tumor", "Normal")))
f2 <- lmer(DvP_Score ~ group + (1 | patient), data=cellMet)
summary(f2)  

cellMet %>%
  group_by(patient, PNI) %>%
  summarize(avg_dvp = mean(DvP_Score)) %>%
  ungroup() %>%
  ggplot(aes(PNI, avg_dvp)) +
  geom_boxplot(aes(fill = PNI)) +
  theme_classic() +
  labs(x = "", y = "Mean DvP Score", title = "DvP Score") +
  theme(plot.title = element_text(hjust = .5)) +
  guides(fill = "none") +
  scale_fill_brewer(palette = "Accent", direction = -1)

dvpPlot <- cellMet %>%
  group_by(patient, condition2) %>%
  summarize(avg_dvp = mean(DvP_Score)) %>%
  ungroup() %>%
  ggplot(aes(condition2, avg_dvp)) +
  geom_boxplot(aes(fill = condition2)) +
  theme_classic() +
  labs(x = "", y = "Mean DvP Score", title = "DvP Score") +
  theme(plot.title = element_text(hjust = .5)) +
  guides(fill = "none") +
  scale_fill_brewer(palette = "Accent", direction = -1)
print(dvpPlot)

ggsave(
  filename = file.path("figures", "manuscript", "DvP_Normal_Tumor_PNI_Boxplot.svg"),
  plot = dvpPlot,
  width = 6,
  height = 4
)

