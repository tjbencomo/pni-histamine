## Description: Draw information flow plot comparing signaling in LR vs
## PNI tumors. Filter to only include KC involved pathways

library(CellChat)
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(RColorBrewer)


## Plotting and helper functions
# Plot information flow plot comparing LR vs PNI signaling presence
plot_comp <- function(res, pthresh = 1e-3, text_size = 18, title = "") {
  # Filter signaling to only include pathways above threshold
  # Need to use this approach as each signaling pathway has 2 rows (LR + PNI) with different p-values
  # Filtering just by p-value may exclude one row but not the other
  pvaldf <- res %>% 
    select(name, group, pvalues) %>% 
    pivot_wider(names_from = "group", values_from = "pvalues") %>%
    mutate(min_pval = min(LR, PNI)) %>%
    filter(min_pval < pthresh)
  res <- res %>% filter(name %in% pvaldf$name)
  # For some reason my filtering messes up the ordering so I need
  # to determine the ordering myself, ranking from most PNI to most LR 
  name_order <- get_order(res)
  # This code was pulled from the rankNet function used by CellChat
  p <- res %>% 
    mutate(name = factor(name, levels = name_order)) %>%
    filter(pvalues < pthresh) %>% 
    ggplot(aes(name, contribution.scaled, fill=group)) + 
    geom_bar(stat = 'identity', position = 'fill') + 
    theme_classic() + 
    labs(x = "", y = "Relative Information Flow") +
    guides(fill = guide_legend(title = "")) +
    coord_flip() + 
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50", size = 0.5) +
    theme(text = element_text(size = text_size), plot.title = element_text(hjust = .5)) +
    scale_fill_manual(values = c("LR" = "#BEAED4", "PNI" = "#7FC97F")) +
    # scale_fill_manual(values = c("LR" = "#00BFC4", "PNI" = "#F8766D")) +
    ggtitle(title)
  print(p)
  return(p)
}

# Determine ordering of pathways to ensure contributions are properly
# sorted
get_order <- function(res) {
  # determine whether a signaling pathway is all PNI/LR or mixed
  catdf <- res %>%
    select(-contribution) %>%
    pivot_wider(
      names_from = "group", 
      values_from = c("contribution.scaled", "pvalues")
    ) %>%
    rowwise() %>%
    mutate(min_contrib = min(contribution.scaled_LR, contribution.scaled_PNI)) %>%
    ungroup() %>%
    mutate(category = case_when(
      min_contrib == 0 & contribution.scaled_PNI > 0 ~ "PNI",
      min_contrib == 0 & contribution.scaled_LR > 0 ~ "LR",
      TRUE ~ "Middle"
    ))
  # Separate into separate groups for within group sorting
  catdf$name <- as.character(catdf$name)
  pnidf <- catdf %>% filter(category == "PNI")
  lrdf <- catdf %>% filter(category == "LR")
  # Note d needs to be divided and not subtracted because the size/ratio
  # of the contribution effects the bar size, not the difference
  middf <- catdf %>% filter(category == "Middle") %>% 
    mutate(d = contribution.scaled_PNI / contribution.scaled_LR)
  
  # Get ordering for each group
  pniOrder <- pnidf %>% arrange(desc(contribution.scaled_PNI)) %>% pull(name)
  lrOrder <- lrdf %>% arrange(contribution.scaled_LR) %>% pull(name)
  midOrder <- middf %>% 
    arrange(desc(d)) %>%
    pull(name)
  o <- c(pniOrder, midOrder, lrOrder)
  return(o)
}


data_dir <- file.path("data")
fig_dir <- file.path("figures", "manuscript")

combo <- readRDS(file.path(data_dir, "merged_cellchat.rds"))


df.kc_out <- rankNet(combo, mode = "comparison", stacked = T, 
                     do.stat = TRUE, 
                     sources.use = levels(combo@idents$joint)[str_detect(levels(combo@idents$joint), "KC")],
                     return.data = T)[["signaling.contribution"]]
df.kc_in <- rankNet(combo, mode = "comparison", stacked = T, 
                    do.stat = TRUE, 
                    targets.use = levels(combo@idents$joint)[str_detect(levels(combo@idents$joint), "KC")],
                    return.data = T)[["signaling.contribution"]]
df.kc_in$type <- "KC_IN"
df.kc_out$type <- "KC_OUT"
# Combine outgoing and incoming KC signaling
# Blacklist - pathways don't match literature/prior knowledge so exclude
blacklist <- c("TGFb", "WNT", "BMP")
sumsigdf <- df.kc_in %>%
  bind_rows(df.kc_out) %>%
  group_by(name, group) %>%
  summarize(
    contribution = mean(contribution),
    contribution.scaled = mean(contribution.scaled),
    pvalues = mean(pvalues)
  ) %>%
  ungroup() %>%
  filter(!(name %in% blacklist))
p <- plot_comp(sumsigdf, pthresh = 1e-4)
p2 <- plot_comp(sumsigdf, pthresh = 1e-8)

pathwaysForTable <- c("PTN", "KIT", "SAA", "AGRN", "CD226", "LIGHT", "PVR", "SEMA5", 
                      "CSF", "CHEMERIN", "ACTIVIN", "NT", "SLURP", "CD137", "NCAM", 
                      "SEMA6", "FLT3")
sumsigdf %>%
  filter(name %in% pathwaysForTable, group != "LR") %>%
  select(name, contribution, contribution.scaled) %>%
  arrange(desc(contribution.scaled)) %>%
  write_csv("data/CellChat_PNI_Hits_Supplemental_Table.csv")

## Save Plots
ggsave(p2, filename = file.path(fig_dir, "Fig2_CellChat_Compare_Conditions.eps"), width = 5, height = 8)

## Alternative plot for Ankit
test_plot <- function(res, pthresh = 1e-3, text_size = 18, title = "") {
  # Filter signaling to only include pathways above threshold
  # Need to use this approach as each signaling pathway has 2 rows (LR + PNI) with different p-values
  # Filtering just by p-value may exclude one row but not the other
  pvaldf <- res %>% 
    select(name, group, pvalues) %>% 
    pivot_wider(names_from = "group", values_from = "pvalues") %>%
    mutate(min_pval = min(LR, PNI)) %>%
    filter(min_pval < pthresh)
  res <- res %>% filter(name %in% pvaldf$name)
  # For some reason my filtering messes up the ordering so I need
  # to determine the ordering myself, ranking from most PNI to most LR 
  name_order <- get_order(res)
  x <- res %>%
    select(-contribution, -pvalues) %>%
    pivot_wider(names_from = "group", values_from = c("contribution.scaled")) %>%
    mutate(total_contrib = PNI  + LR) %>%
    mutate(PNI = PNI  / total_contrib) %>%
    mutate(LR = LR / total_contrib) %>%
    select(-total_contrib) %>%
    pivot_longer(
      !name,
      names_to = "group",
      values_to = "contribution"
    ) %>%
    mutate(name = factor(name, levels = name_order)) %>%
    filter(name %in% pvaldf$name)
  finalOrder <- x %>% filter(group == "PNI") %>%
    arrange(desc(contribution)) %>%
    pull(name)
  x %>%
    mutate(name = factor(name, levels = finalOrder)) %>%
    ggplot(aes(name, contribution)) + 
    geom_point(aes(color=group)) +
    theme_classic() + 
    labs(x = "Pathway", y = "Relative Information Flow") +
    scale_color_manual(values = c("LR" = "#BEAED4", "PNI" = "#7FC97F")) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1),
      text = element_text(size = 14)
    )
}

altPlot <- test_plot(sumsigdf, pthresh = 1e-8)

## See if other LR-specific pathways show up if we don't filter by in/out KC interactions
x <- rankNet(combo, mode="comparison", stacked = T, 
             do.stat = TRUE)

