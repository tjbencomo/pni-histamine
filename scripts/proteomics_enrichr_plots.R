## Description: Plot enrichr results from phospho-proteomics histamine experiments
## Right now only plots 10 minute experiments

library(readr)
library(dplyr)
library(ggplot2)
library(stringr)

# keep = pathways to include even if not in top N genesets
# Filters pathway by combined score, ranks in plot by adjusted p-value
# Cutoff line at .1
# color can be 'green' or 'purple'
plot_enrichr <- function(resdf, title, topN = 10, keep = NULL, color = "green") {
  if (color == "green") {
    color_val <- "#7fc97f"
  } else if (color == "purple") {
    color_val <- "#beaed4"
  } else {
    print("Color specified not green or purple. Using green as default")
    color_val <- "#7fc97f"
  }
  p <- resdf %>%
    rename(combined_score = `Combined Score`) %>%
    rename(padj = "Adjusted P-value") %>%
    mutate(value = -log10(padj)) %>%
    mutate(term = gsub('(.{1,20})(\\s|$)', '\\1\n', Term)) %>%
    mutate(term = str_replace(term, "WP[0-9]+", "")) %>%
    arrange(desc(combined_score)) %>%
    mutate(rank = row_number()) %>%
    filter(rank < topN | Term %in% keep) %>%
    # slice_max(combined_score, n = topn) %>%
    ggplot(aes(value, reorder(term, value))) +
    geom_col(fill = color_val, width = .75) +
    theme_classic() +
    labs(x = "-log10 Adjusted P-Value", y = "") +
    geom_vline(xintercept = -log10(.1), linetype = "dashed", color = "red") +
    theme(text = element_text(size = 12), plot.title = element_text(hjust = .5)) +
    ggtitle(title)
  return(p)
}


data_dir <- "data"
fig_dir <- file.path("figures", "manuscript")
proteomics_dir <- file.path(data_dir, "proteomics", "enrichr-results")

## Load results tables
histUpWiki <- read_tsv(file.path(proteomics_dir, "Histamine_vs_Control_10min_Phospho_WikiPathway_Results.txt"))
histDnWiki <- read_tsv(file.path(proteomics_dir, "Histamine_vs_Control_10min_DePhospho_WikiPathway_Results.txt"))
histUpEls <- read_tsv(file.path(proteomics_dir, "Histamine_vs_Control_10min_Phospho_Elsevier_Results.txt"))
histDnEls <- read_tsv(file.path(proteomics_dir, "Histamine_vs_Control_10min_DePhospho_Elsevier_Results.txt"))
antiHistUpWiki <- read_tsv(file.path(proteomics_dir, "Antihistamine_vs_Histamine_10min_Phospho_WikiPathway_Results.txt"))
antiHistDnWiki <- read_tsv(file.path(proteomics_dir, "Antihistamine_vs_Histamine_10min_DePhospho_WikiPathway_Results.txt"))
antiHistUpEls <- read_tsv(file.path(proteomics_dir, "Antihistamine_vs_Histamine_10min_Phospho_Elsevier_Results.txt"))
antiHistDnEls <- read_tsv(file.path(proteomics_dir, "Antihistamine_vs_Histamine_10min_DePhospho_Elsevier_Results.txt"))
wikiInt <- read_tsv(file.path(proteomics_dir, "Antihistamine_10min_DePhospho_INTERSECT_Histamine_10min_Phospho_WikiPathway_Results.txt"))
elsInt <- read_tsv(file.path(proteomics_dir, "Antihistamine_10min_DePhospho_INTERSECT_Histamine_10min_Phospho_Elsevier_Results.txt"))


## Make plots
histUpWPlot <- plot_enrichr(histUpWiki, "Histamine vs Control Phospho 10min WikiPathways",
                            keep = c("RAC1/PAK1/p38/MMP2 Pathway WP3303", "EGF/EGFR signaling pathway WP437",
                                     "Kit receptor signaling pathway WP304", 
                                     "Brain-derived neurotrophic factor (BDNF) signaling pathway WP2380",
                                     "MAPK Signaling Pathway WP382"),
                            color = "green")

histDnWPlot <- plot_enrichr(histDnWiki, "Histamine vs Control De-Phospho 10min WikiPathways",
                            keep = c("Focal Adhesion WP306", "Integrin-mediated Cell Adhesion WP185",
                                     "Regulation of Actin Cytoskeleton WP51"),
                            color = "green")

histUpEPlot <- plot_enrichr(histUpEls, "Histamine vs Control Phospho 10min Elsevier",
                            keep = c("Neurotrophic Factor Deprivation in Retinal Ganglion Cell Death",
                                     "VEGF Signaling", "P38 MAPK/MAPK14 Signaling",
                                     "Integrins in Cancer Cell Motility, Invasion and Survival",
                                     "AXL Receptor Tyrosine Kinase Signaling", "Tight Junction Assembly (Claudins)",
                                     "Adherens Junction Assembly (Nectin)"),
                            color = "purple")

histDnEPlot <- plot_enrichr(histDnEls, "Histamine vs Control De-Phospho 10min Elsevier", color = "purple")
antiHistUpWPlot <- plot_enrichr(antiHistUpWiki, "Anti-Histamine vs Histamine Phospho 10min WikiPathways", color = "green")
antiHistDnWPlot <- plot_enrichr(antiHistDnWiki, "Anti-Histamine vs Histamine De-Phospho 10min WikiPathways",
                                keep = c("Focal Adhesion WP306", "RAC1/PAK1/p38/MMP2 Pathway WP3303",
                                         "ErbB signaling pathway WP673", "MAPK Signaling Pathway WP382"),
                                color = "green")

antiHistUpEPlot <- plot_enrichr(antiHistUpEls, "Anti-Histamine vs Histamine Phospho 10min Elsevier", color = "purple")
antiHistDnEPlot <- plot_enrichr(antiHistDnEls, "Anti-Histamine vs Histamine De-Phospho 10min Elsevier",
                                keep = c("Neurotrophic Factor Deprivation in Retinal Ganglion Cell Death",
                                         "Integrins in Cancer Cell Motility, Invasion and Survival",
                                        "Desmosome Assembly", "Tight Junction Assembly (Claudins)"),
                                color = "purple")

## Draw plots showing enrichr results for
## intersection of phosphorylated proteins in histamine vs control 10 min and
## de-phosphorylated proteins form antihistamine vs histamine 10 min
wikiKeep <- c("RAC1/PAK1/p38/MMP2 Pathway WP3303", "Leptin signaling pathway WP2034",
                     "EGFR Tyrosine Kinase Inhibitor Resistance WP4806", 
                     "Regulation of Microtubule Cytoskeleton WP2038",
                     "Alpha 6 Beta 4 signaling pathway WP244", 
                     "IL-6 signaling pathway WP364")
wikiIntPlot <- plot_enrichr(wikiInt %>% filter(Term %in% wikiKeep), 
                            title = "Histamine AND Anti-histamine proteins Wikipath", color = "green")
elsKeep <- c("Neurotrophic Factor Deprivation in Retinal Ganglion Cell Death",
                "Integrins in Cancer Cell Motility, Invasion and Survival",
                "Growth Factor Signaling in Neuroblastoma", "FOXO1 Signaling",
                "Tamoxifen Induced Endometrial Cancer", "NTRK1/2/3 -> Acetylcholine Production")
elsIntPlot <- plot_enrichr(elsInt %>% filter(Term %in% elsKeep), title = "Histamine AND Anti-histamine proteins Elsevier", color = "purple")


## Save plots
# ggsave(filename = file.path(fig_dir, "Histamine_vs_Control_10min_Up_Wiki.eps"), plot = histUpWPlot, width = 7, height = 8)
# ggsave(filename = file.path(fig_dir, "Histamine_vs_Control_10min_Down_Wiki.eps"), plot = histDnWPlot, width = 7, height = 8)
# ggsave(filename = file.path(fig_dir, "Histamine_vs_Control_10min_Up_Elsevier.eps"), plot = histUpEPlot, width = 7, height = 8)
# ggsave(filename = file.path(fig_dir, "Histamine_vs_Control_10min_Down_Elsevier.eps"), plot = histDnEPlot, width = 7, height = 8)
# ggsave(filename = file.path(fig_dir, "Antihistamine_vs_Histamine_10min_Up_Wiki.eps"), plot = antiHistUpWPlot, width = 7, height = 8)
# ggsave(filename = file.path(fig_dir, "Antihistamine_vs_Histamine_10min_Down_Wiki.eps"), plot = antiHistDnWPlot, width = 7, height = 8)
# ggsave(filename = file.path(fig_dir, "Antihistamine_vs_Histamine_10min_Up_Elsevier.eps"), plot = antiHistUpEPlot, width = 7, height = 8)
# ggsave(filename = file.path(fig_dir, "Antihistamine_vs_Histamine_10min_Down_Elsevier.eps"), plot = antiHistDnEPlot, width = 7, height = 8)
ggsave(filename = file.path(fig_dir, "Fig5_Intersection_Histamine_10min_Up_and_Antihistamine_10min_Down_Wikipath.eps"),
       plot = wikiIntPlot, width = 7, height = 8)
ggsave(filename = file.path(fig_dir, "Fig5_Intersection_Histamine_10min_Up_and_Antihistamine_10min_Down_Elsevier.eps"),
       plot = elsIntPlot, width = 7, height = 8)


