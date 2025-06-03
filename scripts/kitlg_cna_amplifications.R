## Description: Draw bar plot showing KITLG CNA amplification frequency in LR vs PNI cSCC

library(dplyr)
library(ggplot2)
library(RColorBrewer)


fig_dir <- file.path("figures", "manuscript")
outfp <- file.path(fig_dir, "Fig2_KITLG_CNA_Proportion_Plot.eps")

# Colors
diploidBlue <- "#a6cee3"
diploidGray <- "gray"
ampLightRed <- "#f57f76"
ampRed <- "#ef3b2c"

# Define data
df <- data.frame(
  condition = c(rep("PNI+", 5), rep("PNI-", 10)),
  copyNumber = c(rep("Amplified", 3), rep("Diploid", 12))
)
p <- df %>%
  ggplot(aes(condition)) +
  geom_bar(aes(fill = copyNumber), position = 'fill', width = .7) +
  labs(x = "", y = "Percentage of Tumors", fill = "Copy Number", title = "KITLG") +
  theme_classic() +
  # scale_fill_manual(values = c("Diploid" = diploidBlue, "Amplified" = ampLightRed)) +
  scale_fill_manual(values = c("Diploid" = diploidGray, "Amplified" = ampLightRed)) +
  # scale_fill_manual(values = c("Diploid" = diploidBlue, "Amplified" = ampRed)) +
  theme(
    plot.title = element_text(hjust = .5),
    text = element_text(size = 14)
  ) +
  scale_y_continuous(labels = scales::percent)
print(p)

# Run statistical test
m <- matrix(c(0, 10, 3, 2), nrow = 2)

fisher.test(m)
chisq.test(m)

ggsave(p, filename = outfp, width = 4, height = 6)

## Old barplot code
# df <- data.frame(
#   condition = c("PNI-", "PNI+"),
#   freq = c(0, 3/5)
# )
# p <- df %>%
#   ggplot(aes(condition, freq)) +
#   geom_col(aes(fill = condition), width = .6) +
#   labs(x = "", y = "KITLG CNA Gain Frequency") +
#   scale_fill_brewer(palette = "Accent", direction = -1) +
#   scale_y_continuous(limit = c(0, 1)) +
#   guides(fill = "none") +
#   theme_classic()
# print(p)
