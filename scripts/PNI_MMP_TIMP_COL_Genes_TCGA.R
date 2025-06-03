## Look for DE MMP/TIMP/COLLAGEN genes in PNI+ vs PNI- TCGA tumors to help
## Ankit prioritize genes for validation
## See which MMP/COLLAGEN/TIMP genes correlate with histamine response in TCGA tumors

library(DESeq2)
library(readr)
library(dplyr)
library(ggplot2)
library(RegParallel)
library(pheatmap)
library(ggrepel)
library(GSVA)
library(fgsea)



data_dir <- file.path("data")
fig_dir <- file.path("figures", "manuscript")
rnafp <- file.path(data_dir, "tcga-bulk", "xena_pancancer_vst_normalized_data.rds")
patientsfp <- file.path(data_dir, "tcga-bulk", "PNI_Cohort_Info.csv")

rna <- readRDS(rnafp)
patients <- read_csv(patientsfp)

rna <- rna[, patients$sample]

## Get MMP + TIMP + Collagen genes
allGenes <- rownames(rna)
genes <- allGenes[str_detect(allGenes, "^MMP")]
# genes <- allGenes[str_detect(allGenes, "^MMP|^TIMP|^COL")]
genes <- sort(genes)

rna2 <- rna[genes, ]
rnadata <- data.frame(colData(rna), t(assay(rna2)))
stopifnot(all(rnadata$sampleID == patients$sample))
rnadata$PNI <- patients$perineural_invasion_present
rnadata$cancer_type <- patients$`cancer type abbreviation`
rnadata$PNI <- factor(rnadata$PNI, levels = c("NO", "YES"))

res1 <- RegParallel(
  data = rnadata,
  formula = '[*] ~ PNI',
  FUN = function(formula, data)
    lm(formula = formula,
        data = data),
  FUNtype = 'lm',
  variables = genes,
  blocksize = 20,
  p.adjust = 'fdr')
# adjust for cancer type to check for robustness
res2 <- RegParallel(
  data = rnadata,
  formula = '[*] ~ PNI + cancer_type',
  FUN = function(formula, data)
    lm(formula = formula,
       data = data),
  FUNtype = 'lm',
  variables = genes,
  blocksize = 20,
  p.adjust = 'fdr')


rnadata %>%
  ggplot(aes(PNI, MMP10)) +
  geom_violin(aes(fill = PNI), draw_quantiles = c(.5)) +
  theme_classic() +
  scale_fill_brewer(palette = "Accent", direction = -1) +
  guides(fill = "none")
rnadata %>%
  ggplot(aes(PNI, MMP2)) +
  geom_violin(aes(fill = PNI), draw_quantiles = c(.5)) +
  theme_classic() +
  scale_fill_brewer(palette = "Accent", direction = -1) +
  guides(fill = "none")
rnadata %>%
  ggplot(aes(PNI, MMP9)) +
  geom_violin(aes(fill = PNI), draw_quantiles = c(.5)) +
  theme_classic() +
  scale_fill_brewer(palette = "Accent", direction = -1) +
  guides(fill = "none")
# rnadata %>%
#   ggplot(aes(PNI, COL7A1)) +
#   geom_violin(aes(fill = PNI), draw_quantiles = c(.5)) +
#   theme_classic() +
#   scale_fill_brewer(palette = "Accent", direction = -1) +
#   guides(fill = "none")
# rnadata %>%
#   ggplot(aes(PNI, COL17A1)) +
#   geom_violin(aes(fill = PNI), draw_quantiles = c(.5)) +
#   theme_classic() +
#   scale_fill_brewer(palette = "Accent", direction = -1) +
#   guides(fill = "none")

topGenes <- res2 %>% filter(Term == "PNIYES") %>% slice_max(Beta, n = 20) %>% pull(Variable)
rnadata <- rnadata %>%
  arrange(PNI)
pheatmap(t(rnadata[, topGenes]), annotation_col = rnadata["PNI"], show_colnames = F, scale = "row", cluster_cols = F)

p1 <- res2 %>%
  filter(Term == "PNIYES") %>%
  rename(padj = P.adjust) %>%
  mutate(color_label = case_when(
    padj < .05 & Beta > 0 ~ "Up",
    padj < .05 & Beta < 0 ~ "Down",
    TRUE ~ "NS"
  )) %>%
  group_by(Beta > 0) %>%
  arrange(desc(abs(Beta))) %>%
  mutate(brank = row_number()) %>%
  ungroup() %>%
  mutate(text_label = case_when(
    Beta > 0 & padj < .05 & brank < 30 ~ Variable,
    Beta < 0 & padj < .05 & brank < 10 ~ Variable,
    TRUE ~ ""
  )) %>%
  ggplot(aes(Beta, -log10(padj), label = text_label)) +
  geom_point(aes(color = color_label)) +
  geom_text_repel() +
  scale_color_manual(values = c("Blue", "Gray", "Red")) +
  labs(x = "Log2 Fold Change", y = "-Log10 Adjusted P-Value", color = "FDR < 5%") +
  theme_classic()
print(p1)


# Draw MA plot
genedf <- data.frame(
  gene = genes,
  meanExp = colMeans(rnadata[, genes])
)
p2 <- res2 %>%
  inner_join(genedf, by = c("Variable" = "gene")) %>%
  filter(Term == "PNIYES") %>%
  rename(padj = P.adjust) %>%
  mutate(color_label = case_when(
    padj < .05 & Beta > 0 ~ "Up",
    padj < .05 & Beta < 0 ~ "Down",
    TRUE ~ "NS"
  )) %>%
  group_by(Beta > 0) %>%
  arrange(desc(abs(Beta))) %>%
  mutate(brank = row_number()) %>%
  ungroup() %>%
  mutate(text_label = case_when(
    Beta > 0 & padj < .05 & brank < 30 ~ Variable,
    Beta < 0 & padj < .05 & brank < 10 ~ Variable,
    TRUE ~ ""
  )) %>%
  filter(Beta > 0) %>%
  ggplot(aes(meanExp, Beta, label = text_label)) +
  geom_point(aes(color = color_label)) +
  geom_text_repel() +
  scale_color_manual(values = c("Gray", "Red")) +
  labs(x = "Log2 Mean Expression", y = "Log2 Fold Change", color = "FDR < 5%") +
  theme_classic() +
  scale_y_continuous(limits = c(0, 1.2))
print(p2)
ggsave(p2, filename = file.path("figures", "manuscript", "Fig5_TCGA_MMP_DE.eps"), width = 4, height = 6)
ggsave(p2, filename = file.path("figures", "manuscript", "Fig5_TCGA_MMP_DE.pdf"), width = 4, height = 6)


# gobpfp <- file.path(data_dir, "genesets", "c5.go.bp.v7.4.symbols.gmt")
# gobp <- gmtPathways(gobpfp)
# es <- gsva(assay(rna), list('HistResp' = gobp$GOBP_RESPONSE_TO_HISTAMINE), verbose = T)
# 
# esdf <- data.frame(
#   sampleID = patients$sample,
#   HistResp = es['HistResp', ]
# ) %>% inner_join(rnadata)
# 
# rMat <- esdf %>%
#   dplyr::select(HistResp, all_of(genes)) %>%
#   as.matrix() %>%
#   cor(method = 'spearman')

