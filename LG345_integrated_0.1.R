library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(MAST)

setwd("/athena/ganlab/scratch/lif4001/ChaoWang_IKKb/IKKb_CA/LG345_singlets")

#setwd("~/Desktop/data_analysis/Chao_IKKb")

LG345_integrated <- readRDS("LG345_integrated.rds")
DefaultAssay(LG345_integrated) <- 'integrated'

# LG345_integrated <- NormalizeData(LG345_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
# LG345_integrated <- FindVariableFeatures(LG345_integrated, selection.method = "vst", nfeatures = 3000)

LG345_integrated <- ScaleData(LG345_integrated, verbose = FALSE)
LG345_integrated <- RunPCA(LG345_integrated, features = VariableFeatures(object = LG345_integrated), verbose = FALSE)

LG345_integrated <- FindNeighbors(LG345_integrated, dims = 1:20)
LG345_integrated <- FindClusters(LG345_integrated, resolution = 0.1)
LG345_integrated <- RunUMAP(LG345_integrated, dims = 1: 20)

str(LG345_integrated)

DefaultAssay(LG345_integrated) <- 'RNA'
LG345_integrated <- NormalizeData(LG345_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
LG345_integrated <- ScaleData(LG345_integrated, features = rownames(LG345_integrated))

pdf("LG345_integrated_umap.pdf", width=5, height=4)
DimPlot(LG345_integrated, reduction = 'umap', label = T)
dev.off()
pdf("LG345_integrated_umap_split_individual.pdf", width=16, height=9)
DimPlot(LG345_integrated, reduction = "umap", split.by = "orig.ident", label = T, ncol = 4)
dev.off()
pdf("LG345_integrated_umap_split_Condition.pdf", width=16, height=4)
DimPlot(LG345_integrated, reduction = "umap", split.by = "Condition", label = T)
dev.off()

saveRDS(LG345_integrated, file = 'LG345_integrated_PCA_0.1.rds')

