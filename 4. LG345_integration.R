
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)

#load in data from Cell Ranger or other counts data ====

#for loading Cell Ranger counts:
setwd("/athena/ganlab/scratch/lif4001/ChaoWang_IKKb/IKKb_CA/LG345_singlets")
LG345_83 <- readRDS(file = "LG345_83_singlets_PCA.rds")
LG345_131 <- readRDS(file = "LG345_131_singlets_PCA.rds")
LG345_82 <- readRDS(file = "LG345_82_singlets_PCA.rds")
LG345_132 <- readRDS(file = "LG345_132_singlets_PCA.rds")
LG345_98 <- readRDS(file = "LG345_98_singlets_PCA.rds")
LG345_135 <- readRDS(file = "LG345_135_singlets_PCA.rds")
LG345_84 <- readRDS(file = "LG345_84_singlets_PCA.rds")
LG345_133 <- readRDS(file = "LG345_133_singlets_PCA.rds")

LG345_IKKbWT <- c(LG345_83, LG345_131)
LG345_IKKbCA <- c(LG345_82, LG345_132)
LG345_IKKbWT_P301S <- c(LG345_98, LG345_135)
LG345_IKKbCA_P301S <- c(LG345_84, LG345_133)
anchors_IKKbWT <- FindIntegrationAnchors(object.list = LG345_IKKbWT, dims = 1:30)
LG345_IKKbWT_integrated <- IntegrateData(anchorset = anchors_IKKbWT, dims = 1:30)
anchors_IKKbCA <- FindIntegrationAnchors(object.list = LG345_IKKbCA, dims = 1:30)
LG345_IKKbCA_integrated <- IntegrateData(anchorset = anchors_IKKbCA, dims = 1:30)
anchors_IKKbWT_P301S <- FindIntegrationAnchors(object.list = LG345_IKKbWT_P301S, dims = 1:30)
LG345_IKKbWT_P301S_integrated <- IntegrateData(anchorset = anchors_IKKbWT_P301S, dims = 1:30)
anchors_IKKbCA_P301S <- FindIntegrationAnchors(object.list = LG345_IKKbCA_P301S, dims = 1:30)
LG345_IKKbCA_P301S_integrated <- IntegrateData(anchorset = anchors_IKKbCA_P301S, dims = 1:30)

LG345_all <- c(LG345_IKKbWT_integrated,LG345_IKKbCA_integrated,LG345_IKKbWT_P301S_integrated,LG345_IKKbCA_P301S_integrated)
anchors_all <- FindIntegrationAnchors(object.list = LG345_all, dims = 1:30)
LG345_integrated <- IntegrateData(anchorset = anchors_all, dims = 1:30)

#rm(LG345_83, LG345_131, LG345_82, LG345_132, LG345_98, LG345_135, LG345_84,LG345_133,
#   LG345_IKKbWT, LG345_IKKbCA, LG345_IKKbWT_P301S, LG345_IKKbCA_P301S,
#   LG345_IKKbWT_integrated, LG345_IKKbCA_integrated, LG345_IKKbWT_P301S_integrated, LG345_IKKbCA_P301S_integrated)

#LG343_IKKbKO_integrated[["percent.mt"]] <- PercentageFeatureSet(object = LG343_IKKbKO_integrated, pattern = "^mt-")

pdf("LG345_QC.pdf", width=12, height=4)
Idents(LG345_integrated) <- "orig.ident"
VlnPlot(object = LG345_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
Idents(LG345_integrated) <- "Condition"
VlnPlot(object = LG345_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

saveRDS(LG345_integrated, file = "LG345_integrated.rds")

Idents(LG345_IKKbWT) <- "seurat_clusters"
pdf("LG345_IKKbWT.pdf", width=15, height=6)
p1 <- DimPlot(LG345_IKKbWT, reduction = "umap", group.by = "Condition")
p2 <- DimPlot(LG345_IKKbWT, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
dev.off()

Idents(LG345_IKKbWT_integrated) <- "seurat_clusters"
pdf("LG345_IKKbWT_integrated.pdf", width=15, height=6)
p1 <- DimPlot(LG345_IKKbWT_integrated, reduction = "umap", group.by = "Condition")
p2 <- DimPlot(LG345_IKKbWT_integrated, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
dev.off()

Idents(LG345_IKKbCA) <- "seurat_clusters"
pdf("LG345_IKKbCA.pdf", width=15, height=6)
p1 <- DimPlot(LG345_IKKbCA, reduction = "umap", group.by = "Condition")
p2 <- DimPlot(LG345_IKKbCA, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
dev.off()

Idents(LG345_IKKbCA_integrated) <- "seurat_clusters"
pdf("LG345_IKKbCA_integrated.pdf", width=15, height=6)
p1 <- DimPlot(LG345_IKKbCA_integrated, reduction = "umap", group.by = "Condition")
p2 <- DimPlot(LG345_IKKbCA_integrated, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
dev.off()

Idents(LG345_IKKbWT_P301S) <- "seurat_clusters"
pdf("LG345_IKKbWT_P301S.pdf", width=15, height=6)
p1 <- DimPlot(LG345_IKKbWT_P301S, reduction = "umap", group.by = "Condition")
p2 <- DimPlot(LG345_IKKbWT_P301S, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
dev.off()

Idents(LG345_IKKbWT_P301S_integrated) <- "seurat_clusters"
pdf("LG345_IKKbWT_P301S_integrated.pdf", width=15, height=6)
p1 <- DimPlot(LG345_IKKbWT_P301S_integrated, reduction = "umap", group.by = "Condition")
p2 <- DimPlot(LG345_IKKbWT_P301S_integrated, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
dev.off()

Idents(LG345_IKKbCA_P301S) <- "seurat_clusters"
pdf("LG345_IKKbCA_P301S.pdf", width=15, height=6)
p1 <- DimPlot(LG345_IKKbCA_P301S, reduction = "umap", group.by = "Condition")
p2 <- DimPlot(LG345_IKKbCA_P301S, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
dev.off()

Idents(LG345_IKKbCA_P301S_integrated) <- "seurat_clusters"
pdf("LG345_IKKbCA_P301S_integrated.pdf", width=15, height=6)
p1 <- DimPlot(LG345_IKKbCA_P301S_integrated, reduction = "umap", group.by = "Condition")
p2 <- DimPlot(LG345_IKKbCA_P301S_integrated, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
dev.off()

Idents(LG345_integrated) <- "seurat_clusters"
pdf("LG345_integrated.pdf", width=15, height=6)
p1 <- DimPlot(LG345_integrated, reduction = "umap", group.by = "Condition")
p2 <- DimPlot(LG345_integrated, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
dev.off()
