## Redo TCGA PNI stress geneset enrichment analysis using GSEA instead of GSVA

library(DESeq2)
library(readr)
library(dplyr)
library(ggplot2)
library(RegParallel)
library(fgsea)

plotRes <- function(geneset, ranks) {
  p <- plotEnrichment(geneset, ranks) + 
    labs(
      x = "Rank", y = "Enrichment Score"
    )
  p$layers[[5]]$aes_params$colour <- 'black'
  # p$layers[[5]]$aes_params$size=1.5
  p$layers[[2]]$aes_params$linetype = "solid"
  p$layers[[3]]$aes_params$size = 0
  p$layers[[3]]$aes_params$colour = "white"
  p$layers[[5]]$aes_params$size=1.5
  p$layers[[6]]$aes_params$size = .7
  p <- p + scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) +
    theme(
      axis.title.x = element_text(size = 16),
      axis.title.y = element_text(size = 16),
      axis.text.x = element_text(size = 0, colour = "white"), 
      axis.ticks.x = element_line(size = 0),
      axis.text.y = element_text(size = 14)
    )
  return(p)
}



data_dir <- file.path("data")
fig_dir <- file.path("figures", "manuscript")
rnafp <- file.path(data_dir,  "xena_tcga_pancancer_vst_normalized_data.rds")
patientsfp <- file.path(data_dir,  "TCGA_PNI_Cohort_Info.csv")

rna <- readRDS(rnafp)
patients <- read_csv(patientsfp)

rna <- rna[, patients$sample]

rnadata <- data.frame(colData(rna), t(assay(rna)))
stopifnot(all(rnadata$sampleID == patients$sample))
rnadata$PNI <- patients$perineural_invasion_present
rnadata$cancer_type <- patients$`cancer type abbreviation`
rnadata$PNI <- factor(rnadata$PNI, levels = c("NO", "YES"))
genes <-rownames(rna)[rownames(rna) %in% colnames(rnadata)]


res1 <- RegParallel(
  data = rnadata,
  formula = '[*] ~ PNI',
  FUN = function(formula, data)
    lm(formula = formula,
       data = data),
  FUNtype = 'lm',
  variables = genes,
  blocksize = 500,
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
  blocksize = 500,
  p.adjust = 'fdr')

set.seed(42)
ranks <- res1$t
names(ranks) <- res1$Variable
ranks2 <- res2 %>% filter(Term == "PNIYES") %>% pull(t)
names(ranks2) <- res2 %>% filter(Term == "PNIYES") %>% pull(Variable)
stress <- scan('data/puram_stress.txt', what = character())
gs <- list('Stress' = stress)
x1 <- fgsea(gs, ranks, eps = 0)
x2 <- fgsea(gs, ranks2, eps = 0)

p <- plotRes(stress, ranks)
ggsave(p, filename = file.path(fig_dir, "Fig1_Stress_Enrichment_In_TCGA_GSEA.eps"), width = 6, height = 5)
