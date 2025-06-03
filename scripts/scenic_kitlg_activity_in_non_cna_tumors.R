## Description: Compare KITLG regulator activity from SCENIC in PNI+ KITLG amplified vs diploid tumors 
## and vs PNI- tumors
## Amplified samples found in HMM_CNV_predictions.HMMi6.hmm_mode-samples.Pnorm_0.5.pred_cnv_genes.dat in
## data/infercnv/results/atwood-control/PNI_TSK_Results

library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(lme4)
library(pbkrtest)
library(patchwork)

fig_dir <- file.path("figures", "manuscript")
auc <- read_csv("data/scenic/LR_PNI_Tumor_Cells_SCENIC.csv")
regulons <- read_csv("data/scenic/LR_PNI_Tumor_Cell_Regulons.csv")
met <- read_csv("data/cell_metadata.csv.gz")
ampPatients <- c("LEE01", "LEE03", "LEE05")
dipPatients <- c("LEE02", "LEE04")

df <- auc %>%
  mutate(full_barcode = str_replace(full_barcode, "-[0-9]+-[0-9]+", "")) %>%
  mutate(full_barcode = str_c(full_barcode, "-1", sep = "")) %>%
  inner_join(met) %>%
  mutate(subpop = case_when(
    fine_type %in% c("KC_Inflam", "KC_PSK-2") ~ "KC_Basal",
    TRUE ~ fine_type
  )) %>%
  mutate(subpop = str_replace(subpop, "KC_", "")) %>%
  mutate(CNA_Status = ifelse(patient %in% ampPatients, "Amplified", "Diploid")) %>%
  mutate(group = case_when(
    CNA_Status == "Amplified" ~ "KITLG Amplified\nPNI+",
    CNA_Status == "Diploid" & PNI == "PNI+" ~ "KITLG Diploid\nPNI+",
    CNA_Status == "Diploid" & PNI == "PNI-" ~ "KITLG Diploid\nPNI-"
  )) %>%
  mutate(group = factor(group, levels = c("KITLG Diploid\nPNI-", "KITLG Diploid\nPNI+", "KITLG Amplified\nPNI+")))


## Highlight TFs that control KTILG
colorTfs <- c("E2F7", "E2F1", "IRF2", "STAT1", "IRF1")
regPlot <- regulons %>%
  filter(str_detect(TargetGenes, "KITLG")) %>%
  mutate(colorLabel = ifelse(TF %in% colorTfs, "highlight", "other")) %>%
  ggplot(aes(reorder(TF, -NES), NES)) +
  geom_boxplot(aes(fill = colorLabel)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = .5),
    plot.title = element_text(hjust = .5)
  ) +
  labs(x = "", y = "Normalized Enrichment Score", title = "SCENIC Predicted KITLG Regulators") +
  scale_fill_manual(values = c("highlight" = "#dc6961", "other" = "gray")) +
  guides(fill = "none")
print(regPlot)

regPlotShort <- regulons %>%
  filter(str_detect(TargetGenes, "KITLG")) %>%
  filter(TF %in% c("E2F7", "E2F1", "IRF2", "STAT1", "IRF1")) %>% # to shorten - Ankit ask
  mutate(colorLabel = ifelse(TF %in% colorTfs, "highlight", "other")) %>%
  ggplot(aes(reorder(TF, -NES), NES)) +
  geom_boxplot(aes(fill = colorLabel)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = .5),
    plot.title = element_text(hjust = .5)
  ) +
  labs(x = "", y = "Normalized Enrichment Score", title = "SCENIC Predicted KITLG Regulators") +
  # scale_fill_manual(values = c("highlight" = "#dc6961", "other" = "gray")) +
  guides(fill = "none")
print(regPlotShort)

ggsave(
  filename = file.path(fig_dir, "SCENIC_Top_KITLG_Regulator_Ranks.svg"),
  plot = regPlotShort,
  width = 6,
  height = 5
)

## Only look at PNI+ tumors
df %>%
  filter(PNI == "PNI+") %>%
  ggplot(aes(CNA_Status, `IRF2(+)`)) +
  geom_violin(aes(fill = CNA_Status)) +
  theme_classic()
df %>%
  filter(PNI == "PNI+") %>%
  ggplot(aes(CNA_Status, `IRF2(+)`)) +
  geom_violin(aes(fill = CNA_Status)) +
  theme_classic() +
  facet_wrap(~subpop)
df %>%
  filter(PNI == "PNI+", subpop == "TSK") %>%
  ggplot(aes(CNA_Status, `IRF2(+)`)) +
  geom_violin(aes(fill = CNA_Status)) +
  theme_classic()
t.test(`IRF2(+)` ~ CNA_Status, data = df %>% filter(PNI == "PNI+", subpop == "TSK"))
wilcox.test(`IRF2(+)` ~ CNA_Status, data = df %>% filter(PNI == "PNI+", subpop == "TSK"))
df %>%
  filter(PNI == "PNI+", subpop == "TSK") %>%
  group_by(patient, CNA_Status) %>%
  summarize(IRF2 = mean(`IRF2(+)`)) %>%
  ungroup() %>%
  ggplot(aes(CNA_Status, IRF2)) +
  geom_violin(aes(fill = CNA_Status)) +
  theme_classic()
df %>%
  filter(PNI == "PNI+", subpop == "TSK") %>%
  group_by(patient, CNA_Status) %>%
  summarize(IRF2 = mean(`IRF2(+)`)) %>%
  ungroup() %>%
  ggplot(aes(CNA_Status, IRF2)) +
  geom_point(aes(color = CNA_Status)) +
  theme_classic()
pni.cna.fit <- lmer(`IRF2(+)` ~ CNA_Status + (1 | patient), data = df %>% filter(PNI == "PNI+", subpop == "TSK"))
summary(pni.cna.fit)

## Also compare to PNI- tumors
## This is currently used in the manuscript

## Change condition ordering for easier comparison
df$group2 <- factor(df$group, levels = c("KITLG Diploid\nPNI+", "KITLG Diploid\nPNI-", "KITLG Amplified\nPNI+"))
df %>%
  ggplot(aes(group, `IRF2(+)`)) +
  geom_violin(aes(fill = group)) +
  theme_classic() +
  facet_wrap(~subpop)
p1 <- df %>%
  filter(subpop == "TSK") %>%
  mutate(z = scale(`IRF2(+)`)) %>%
  ggplot(aes(group, z)) +
  geom_violin(aes(fill = group)) +
  geom_boxplot(width=0.1, color="black", alpha=0.2, outlier.shape = NA) +
  theme_classic() +
  labs(x = "", y = "IRF2 Activity Z-Score", title = "IRF2 Activity") +
  guides(fill = "none") +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 14)
  )
print(p1)
irf2.mem.fit <- lmer(`IRF2(+)` ~ group2 + (1 | patient), data=df %>% filter(subpop == "TSK"))
summary(irf2.mem.fit)
anova(irf2.mem.fit)

p2 <- df %>%
  filter(subpop == "TSK") %>%
  mutate(z = scale(`STAT1(+)`)) %>%
  ggplot(aes(group, z)) +
  geom_violin(aes(fill = group)) +
  geom_boxplot(width=0.1, color="black", alpha=0.2, outlier.shape = NA) +
  theme_classic() +
  labs(x = "", y = "STAT1 Activity Z-Score", title = "STAT1 Activity") +
  guides(fill = "none") +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 14)
  )
print(p2)
stat1.mem.fit <- lmer(`STAT1(+)` ~ group2 + (1 | patient), data=df %>% filter(subpop == "TSK"))
summary(stat1.mem.fit)
anova(stat1.mem.fit)



p3 <- df %>%
  filter(subpop == "TSK") %>%
  mutate(z = scale(`E2F7(+)`)) %>%
  ggplot(aes(group, z)) +
  geom_violin(aes(fill = group)) +
  geom_boxplot(width=0.1, color="black", alpha=0.2, outlier.shape = NA) +
  theme_classic() +
  labs(x = "", y = "E2F7 Activity Z-Score", title = "E2F7 Activity") +
  guides(fill = "none") +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 14)
  )
print(p3)
e2f7.mem.fit <- lmer(`E2F7(+)` ~ group2 + (1 | patient), data=df %>% filter(subpop == "TSK"))
summary(e2f7.mem.fit)
anova(e2f7.mem.fit)


p4 <- df %>%
  filter(subpop == "TSK") %>%
  mutate(z = scale(`E2F1(+)`)) %>%
  ggplot(aes(group, z)) +
  geom_violin(aes(fill = group)) +
  geom_boxplot(width=0.1, color="black", alpha=0.2, outlier.shape = NA) +
  theme_classic() +
  labs(x = "", y = "E2F1 Activity Z-Score", title = "E2F1 Activity") +
  guides(fill = "none") +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 14)
  )
print(p4)
e2f1.mem.fit <- lmer(`E2F1(+)` ~ group2 + (1 | patient), data=df %>% filter(subpop == "TSK"))
summary(e2f1.mem.fit)
anova(e2f1.mem.fit)


p5 <- df %>%
  filter(subpop == "TSK") %>%
  mutate(z = scale(`IRF1(+)`)) %>%
  ggplot(aes(group, z)) +
  geom_violin(aes(fill = group)) +
  geom_boxplot(width=0.1, color="black", alpha=0.2, outlier.shape = NA) +
  theme_classic() +
  labs(x = "", y = "IRF1 Activity Z-Score", title = "IRF1 Activity") +
  guides(fill = "none") +
  theme(
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 14)
  )
print(p5)
irf1.mem.fit <- lmer(`IRF1(+)` ~ group2 + (1 | patient), data=df %>% filter(subpop == "TSK"))
summary(irf1.mem.fit)
anova(irf1.mem.fit)



topHitsViolinPlots <- p1 | p2 | p3 | p4 | p5
ggsave(
  filename = file.path(fig_dir, "SCENIC_Top_KITLG_Regulator_PNI_Scores.svg"), 
  plot = topHitsViolinPlots, 
  width = 24,
  height = 4
)

# df %>% ggplot(aes(`IRF2(+)`)) + geom_density(aes(fill = group), alpha = .5) + theme_bw()

irf2.fit <- lmer(`IRF2(+)` ~ group + (1 | patient), data = df %>% filter(subpop == "TSK"))
summary(irf2.fit)
anova(irf2.fit)


stat1.fit <- lmer(`STAT1(+)` ~ group + (1 | patient), data = df %>% filter(subpop == "TSK"))
summary(stat1.fit)
anova(stat1.fit)



ggsave(filename = file.path(fig_dir, "Fig2_IRF2_Activity_By_KITLG_CNA_Status.svg"), plot = p1, width = 5, height = 5)
ggsave(filename = file.path(fig_dir, "Fig2_STAT1_Activity_By_KITLG_CNA_Status.eps"), plot = p2, width = 5, height = 5)
ggsave(filename = file.path(fig_dir, "Fig2_STAT1_Activity_By_KITLG_CNA_Status.pdf"), plot = p2, width = 5, height = 5)

ggsave(filename = file.path(fig_dir, "FigS3_SCENIC_KITLG_Regulators.eps"), plot = regPlot, width = 10, height = 5)
ggsave(filename = file.path(fig_dir, "FigS3_SCENIC_KITLG_Regulators.pdf"), plot = regPlot, width = 10, height = 5)
