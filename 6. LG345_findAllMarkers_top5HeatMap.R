
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(MAST)

setwd("/athena/ganlab/scratch/lif4001/ChaoWang_IKKb/IKKb_CA/LG345_singlets")

LG345_integrated <- readRDS("LG345_integrated_PCA_0.1.rds")

DefaultAssay(LG345_integrated) <- 'RNA'

LG345_markers <- FindAllMarkers(LG345_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, test.use = "MAST")
write.csv(LG345_markers, "LG345_markers.csv")

LG345_markers <- read.csv(file = "LG345_markers.csv", header=T,row.names =1)
top5 <- LG345_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)
top5$gene <- as.character(top5$gene)
pdf("LG345_HeatMapTop5_0.1_new.pdf", width=24, height=16)
DoHeatmap(LG345_integrated, features = top5$gene) + NoLegend()
dev.off()






