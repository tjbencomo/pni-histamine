## Description: Plot figures for KITLG/HRH1 expression, mast cell product receptors,  
## mast cell infiltration, P38 signaling, and Puram stress signature
## in TCGA PNI cohort. 

library(DESeq2)
library(readr)
library(dplyr)
library(ggplot2)
library(GSVA)
library(patchwork)
library(fgsea)



data_dir <- file.path("data")
fig_dir <- file.path("figures", "manuscript")
rnafp <- file.path(data_dir, "xena_tcga_pancancer_vst_normalized_data.rds")
patientsfp <- file.path(data_dir, "TCGA_PNI_Cohort_Info.csv")
mastgenefp <- file.path(data_dir, "mast_cell_markers.txt")
stressfp <- file.path(data_dir, "puram_stress.txt")
biocartafp <- file.path(data_dir, "c2.cp.biocarta.v7.4.symbols.gmt")



rna <- readRDS(rnafp)
patients <- read_csv(patientsfp)
mastGenes <- scan(mastgenefp, what = character())
stressGenes <- scan(stressfp, what = character())
bc <- gmtPathways(biocartafp)

rna <- rna[, patients$sample]

gs <- list('MCI' = mastGenes, 'Stress' = stressGenes, 'P38' = bc$BIOCARTA_P38MAPK_PATHWAY)
es <- gsva(assay(rna), gs, verbose = T)

esdf <- data.frame(
  sample_id = patients$sample,
  hasPNI = patients$perineural_invasion_present,
  cancer_type = patients$`cancer type abbreviation`,
  KITLG = assay(rna)["KITLG", ],
  HRH1 = assay(rna)["HRH1", ],
  HRH2 = assay(rna)["HRH2", ],
  HRH3 = assay(rna)["HRH3", ],
  HRH4 = assay(rna)["HRH4", ],
  IL17RA = assay(rna)["IL17RA", ],
  FLT1 = assay(rna)["FLT1", ],
  F2RL1 = assay(rna)["F2RL1", ],
  MCI = es["MCI", ],
  Stress = es["Stress", ],
  P38 = es['P38', ]
) %>% as_tibble()

## Draw plots
hrh1Plot <- esdf %>%
  ggplot(aes(hasPNI, HRH1)) +
  geom_violin(aes(fill = hasPNI)) +
  geom_boxplot(width=0.1, color="black", alpha=0.2, outlier.shape = NA) +  
  theme_classic() +
  scale_fill_brewer(palette = "Accent", direction = -1) +
  labs(x = "PNI Present", y = "HRH1 mRNA", title = "HRH1") +
  guides(fill = "none") +
  theme(
    plot.title = element_text(hjust = .5)
  )
print(hrh1Plot)
t.test(HRH1 ~ hasPNI, data = esdf)
wilcox.test(HRH1 ~ hasPNI, data = esdf)
lm(HRH1 ~ hasPNI + cancer_type, data = esdf) %>% summary()

## Draw plots
hrh2Plot <- esdf %>%
  ggplot(aes(hasPNI, HRH2)) +
  geom_violin(aes(fill = hasPNI)) +
  geom_boxplot(width=0.1, color="black", alpha=0.2, outlier.shape = NA) +  
  theme_classic() +
  scale_fill_brewer(palette = "Accent", direction = -1) +
  labs(x = "PNI Present", y = "HRH2 mRNA", title = "HRH2") +
  guides(fill = "none") +
  theme(
    plot.title = element_text(hjust = .5)
  )
print(hrh2Plot)
t.test(HRH2 ~ hasPNI, data = esdf)
wilcox.test(HRH2 ~ hasPNI, data = esdf)
lm(HRH2 ~ hasPNI + cancer_type, data = esdf) %>% summary()

hrh3Plot <- esdf %>%
  ggplot(aes(hasPNI, HRH3)) +
  geom_violin(aes(fill = hasPNI)) +
  geom_boxplot(width=0.1, color="black", alpha=0.2, outlier.shape = NA) +  
  theme_classic() +
  scale_fill_brewer(palette = "Accent", direction = -1) +
  labs(x = "PNI Present", y = "HRH3 mRNA", title = "HRH3") +
  guides(fill = "none") +
  theme(
    plot.title = element_text(hjust = .5)
  )
print(hrh3Plot)
t.test(HRH3 ~ hasPNI, data = esdf)
wilcox.test(HRH3 ~ hasPNI, data = esdf)
lm(HRH3 ~ hasPNI + cancer_type, data = esdf) %>% summary()

hrh4Plot <- esdf %>%
  ggplot(aes(hasPNI, HRH4)) +
  geom_violin(aes(fill = hasPNI)) +
  geom_boxplot(width=0.1, color="black", alpha=0.2, outlier.shape = NA) +  
  theme_classic() +
  scale_fill_brewer(palette = "Accent", direction = -1) +
  labs(x = "PNI Present", y = "HRH4 mRNA", title = "HRH4") +
  guides(fill = "none") +
  theme(
    plot.title = element_text(hjust = .5)
  )
print(hrh4Plot)
t.test(HRH4 ~ hasPNI, data = esdf)
wilcox.test(HRH4 ~ hasPNI, data = esdf)
lm(HRH4 ~ hasPNI + cancer_type, data = esdf) %>% summary()

kitlgPlot <- esdf %>%
  ggplot(aes(hasPNI, KITLG)) +
  geom_violin(aes(fill = hasPNI)) +
  geom_boxplot(width=0.1, color="black", alpha=0.2, outlier.shape = NA) +  
  theme_classic() +
  scale_fill_brewer(palette = "Accent", direction = -1) +
  labs(x = "PNI Present", y = "KITLG mRNA") +
  guides(fill = "none")
print(kitlgPlot)
t.test(KITLG ~ hasPNI, data = esdf)
wilcox.test(KITLG ~ hasPNI, data = esdf)
lm(KITLG ~ hasPNI + cancer_type, data = esdf) %>% summary()


mciPlot <- esdf %>%
  ggplot(aes(hasPNI, MCI)) +
  geom_violin(aes(fill = hasPNI)) +
  geom_boxplot(width=0.1, color="black", alpha=0.2, outlier.shape = NA) +  
  theme_classic() +
  scale_fill_brewer(palette = "Accent", direction = -1) +
  labs(x = "PNI Present", y = "Mast Cell Infiltration") +
  guides(fill = "none") +
  scale_y_continuous(limits = c(-.8, .8), breaks = seq(from = -.8, to = .8, by = .4))
print(mciPlot)
t.test(MCI ~ hasPNI, data = esdf)
wilcox.test(MCI ~ hasPNI, data = esdf)
lm(MCI ~ hasPNI + cancer_type, data = esdf) %>% summary()


stressPlot <- esdf %>%
  ggplot(aes(hasPNI, Stress)) +
  geom_violin(aes(fill = hasPNI)) +
  geom_boxplot(width=0.1, color="black", alpha=0.2, outlier.shape = NA) +  
  theme_classic() +
  scale_fill_brewer(palette = "Accent", direction = -1) +
  labs(x = "PNI Present", y = "Stress Score") +
  guides(fill = "none") +
  scale_y_continuous(limits = c(-.7, .8))
print(stressPlot)
t.test(Stress ~ hasPNI, data = esdf)
wilcox.test(Stress ~ hasPNI, data = esdf)
lm(Stress ~ hasPNI + cancer_type, data = esdf) %>% summary()

il17raPlot <- esdf %>%
  ggplot(aes(hasPNI, IL17RA)) +
  geom_violin(aes(fill = hasPNI)) +
  geom_boxplot(width=0.1, color="black", alpha=0.2, outlier.shape = NA) +  
  theme_classic() +
  scale_fill_brewer(palette = "Accent", direction = -1) +
  labs(x = "PNI Present", y = "IL17RA mRNA", title = "IL17RA") +
  guides(fill = "none") +
  theme(
    plot.title = element_text(hjust = .5)
  )
print(il17raPlot)
t.test(IL17RA ~ hasPNI, data = esdf)
wilcox.test(IL17RA ~ hasPNI, data = esdf)
lm(IL17RA ~ hasPNI + cancer_type, data = esdf) %>% summary()

flt1Plot <- esdf %>%
  ggplot(aes(hasPNI, FLT1)) +
  geom_violin(aes(fill = hasPNI)) +
  geom_boxplot(width=0.1, color="black", alpha=0.2, outlier.shape = NA) +  
  theme_classic() +
  scale_fill_brewer(palette = "Accent", direction = -1) +
  labs(x = "PNI Present", y = "FLT1 mRNA", title = "FLT1") +
  guides(fill = "none") +
  theme(
    plot.title = element_text(hjust = .5)
  )
print(flt1Plot)
t.test(FLT1 ~ hasPNI, data = esdf)
wilcox.test(FLT1 ~ hasPNI, data = esdf)
lm(FLT1 ~ hasPNI + cancer_type, data = esdf) %>% summary()

f2rl1Plot <- esdf %>%
  ggplot(aes(hasPNI, F2RL1)) +
  geom_violin(aes(fill = hasPNI)) +
  geom_boxplot(width=0.1, color="black", alpha=0.2, outlier.shape = NA) +  
  theme_classic() +
  scale_fill_brewer(palette = "Accent", direction = -1) +
  labs(x = "PNI Present", y = "F2RL1 mRNA", title = "F2RL1") +
  guides(fill = "none") +
  theme(
    plot.title = element_text(hjust = .5)
  )
print(f2rl1Plot)
t.test(F2RL1 ~ hasPNI, data = esdf)
wilcox.test(F2RL1 ~ hasPNI, data = esdf)
lm(F2RL1 ~ hasPNI + cancer_type, data = esdf) %>% summary()

p38Plot <- esdf %>%
  ggplot(aes(hasPNI, P38)) +
  geom_violin(aes(fill = hasPNI)) +
  geom_boxplot(width=0.1, color="black", alpha=0.2, outlier.shape = NA) +  
  theme_classic() +
  scale_fill_brewer(palette = "Accent", direction = -1) +
  labs(x = "PNI Present", y = "P38 Score", title = "P38-MAPK Pathway") +
  guides(fill = "none") +
  theme(
    plot.title = element_text(hjust = .5)
  )
print(p38Plot)
t.test(P38 ~ hasPNI, data = esdf)
wilcox.test(P38 ~ hasPNI, data = esdf)
lm(P38 ~ hasPNI + cancer_type, data = esdf) %>% summary()


## Combine HRH1 and KITLG into single plot
finalPlot <- kitlgPlot | mciPlot | hrh1Plot
hrhPlots <-  hrh2Plot | hrh3Plot | hrh4Plot
print(hrhPlots)
otherMediatorPlots <- il17raPlot | flt1Plot | f2rl1Plot
print(otherMediatorPlots)
# genePlot <- hrh1Plot | kitlgPlot
# print(genePlot)

## Save plots
ggsave(hrh1Plot, filename = file.path(fig_dir, "Fig3_TCGA_PNI_HRH1_Expression.eps"), width = 4, height = 6)
ggsave(kitlgPlot, filename = file.path(fig_dir, "Fig3_TCGA_PNI_KITLG_Expression.eps"), width = 4, height = 6)
ggsave(mciPlot, filename = file.path(fig_dir, "Fig3_TCGA_PNI_Mast_Infiltration.eps"), width = 4, height = 6)
# ggsave(genePlot, filename = file.path(fig_dir, "Fig3_TCGA_PNI_HRH1_KITLG_Expression.eps"), width = 6, height = 6)
# ggsave(stressPlot, filename = file.path(fig_dir, "FigS2_TCGA_PNI_Stress.eps"), width = 4, height = 6)
ggsave(finalPlot, filename = file.path(fig_dir, "Fig3_TCGA_PNI_HRH1_KITLG_Mast_Infiltration_Plots.eps"), width = 10, height = 4)
ggsave(hrhPlots, filename = file.path(fig_dir, "FigS4_TCGA_PNI_HRH_Plots.eps"), width = 10, height = 4)
ggsave(p38Plot, filename = file.path(fig_dir, "Fig5_TCGA_PNI_P38_Plots.eps"), width = 4, height = 6)
ggsave(stat1Plot, filename = file.path(fig_dir, "Fig5_TCGA_PNI_STAT1_Plot.eps"), width = 4, height = 6)

