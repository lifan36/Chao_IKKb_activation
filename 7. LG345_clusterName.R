library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(MAST)

setwd("/athena/ganlab/scratch/lif4001/ChaoWang_IKKb/IKKb_CA/LG345_singlets")

LG345_integrated <- readRDS("LG345_integrated_PCA_0.1.rds")

DefaultAssay(LG345_integrated) <- 'RNA'
pdf("LG345_umap_test.pdf", width=8, height=6)
DimPlot(LG345_integrated, reduction = 'umap', label = T)
dev.off()

#Add marker genes

#Neurons
DotPlot(object = LG345_integrated, features = c("Map2")) + RotatedAxis()
#excitatory neurons
sig_EN<-c("Slc17a7", "Camk2a", "Nrgn")
markers.to.plot <- as.matrix(sig_EN)
pdf("LG345_annotation.pdf", width=3.5, height=6)
DotPlot(object = LG345_integrated, features = "Map2") + RotatedAxis()
DotPlot(object = LG345_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()


#inhibitory neurons
sig_IN<-c("Gad1", "Gad2")
markers.to.plot <- as.matrix(sig_IN)

DotPlot(object = LG345_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()


#oligodendrocytes
sig_OL<-c("Plp1", "Mbp", "Mobp")
markers.to.plot <- as.matrix(sig_OL)

DotPlot(object = LG345_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()


#OPCs
sig_OPC<-c("Scrg1", "Pdgfra", "Olig1")
markers.to.plot <- as.matrix(sig_OPC)

DotPlot(object = LG345_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()


#endothelial cells
sig_ENC<-c("Vtn", "Mgp", "Igfbp7")
markers.to.plot <- as.matrix(sig_ENC)
DotPlot(object = LG345_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()

#microglia
sig_MICROG<-c("Cx3cr1", "P2ry12", "Csf1r")
markers.to.plot <- as.matrix(sig_MICROG)
DotPlot(object = LG345_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()

#astrocytes
sig_AST<-c("Clu", "Aldoc", "Pla2g7")
markers.to.plot <- as.matrix(sig_AST)
DotPlot(object = LG345_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()

#T lymphocytes
#sig_TL<-c("CD3D", "CD3G", "CCL5")
sig_TL<-c("Cd3d", "Cd3g", "Ccl5")
markers.to.plot <- as.matrix(sig_TL)
DotPlot(object = LG345_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()
dev.off()
