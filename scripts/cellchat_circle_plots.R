## Description: Draw cellchat circle plots for KIT, SEMA5, PDGF, and PTN
## Probably going in Figure 2


library(CellChat)

data_dir <- file.path("data")
fig_dir <- file.path("figures", "manuscript")

pni <- readRDS(file.path(data_dir, "pni.rds"))
lr <- readRDS(file.path(data_dir, "lowrisk_short_tme.rds"))
kitPlot <- netVisual_aggregate(pni, signaling = "KIT", layout = "circle")


svg(file.path(fig_dir, "Fig2_KIT_Circle.svg"))
kitPlot
dev.off()


