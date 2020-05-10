library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(MAST)

setwd("/athena/ganlab/scratch/lif4001/ChaoWang_IKKb/IKKb_CA/LG345_singlets")

LG345_integrated <- readRDS("LG345_integrated_ready_4_DEGs.rds")
LG345_integrated$celltype.orig.ident <- paste(Idents(LG345_integrated), LG345_integrated$orig.ident, sep = "_")
LG345_integrated$celltype <- Idents(LG345_integrated)

setwd("/athena/ganlab/scratch/lif4001/ChaoWang_IKKb/IKKb_CA/DEGs")

#pdf("UMAP_0.1.pdf", width=8, height=6)
#DimPlot(LG345_integrated, label = TRUE)
#dev.off()
#pdf("UMAP_0.1_by_condition.pdf", width=16, height=4)
#DimPlot(LG345_integrated, split.by = "condition", label = TRUE)
#dev.off()

#Subset Neurons: EN and IN
Cluster_EN_IN <- subset(LG345_integrated, idents = c("excitatory neurons", "inhibitory neurons"))
#Subset Microglia
Cluster_MG <- subset(LG345_integrated, idents = "microglia")


#Test Cluster_EN_IN and Cluster_MG
markers.to.plot <- c("Plp1", "Mbp", "Mobp","Slc17a7", "Nrgn", "Gad1", "Gad2", "Clu", "Aldoc", 
                     "Pla2g7", "Cx3cr1", "P2ry12", "Csf1r","Scrg1", "Pdgfra", "Vtn", "Igfbp7",
                     "Bnc2", "Slc47a1", "Ttr")

pdf("Verify_Subset_neurons_MG.pdf", width=12, height=8)
DotPlot(Cluster_EN_IN, features = rev(markers.to.plot), cols = c("orange", "plum", "pink", "grey"), dot.scale = 8, 
        split.by = "condition") + RotatedAxis()
DotPlot(Cluster_MG, features = rev(markers.to.plot), cols = c("orange", "plum", "pink", "grey"), dot.scale = 8, 
        split.by = "condition") + RotatedAxis()
dev.off()

Idents(Cluster_EN_IN) <- "Condition"
Idents(Cluster_MG) <- "Condition"

#IKKbWT vs IKKbCA DEGs
IKKbWT_CA_EN_IN_DEGs <- FindMarkers(Cluster_EN_IN, ident.1 = "IKKbWT", ident.2 = "IKKbCA", logfc.threshold = 0,
                                      test.use = "MAST")
write.csv(IKKbWT_CA_EN_IN_DEGs, "IKKbWT_CA_EN_IN_DEGs.csv")

IKKbWT_CA_MG_DEGs <- FindMarkers(Cluster_MG, ident.1 = "IKKbWT", ident.2 = "IKKbCA", logfc.threshold = 0,
                                   test.use = "MAST")
write.csv(IKKbWT_CA_MG_DEGs, "IKKbWT_CA_MG_DEGs.csv")

#IKKbWT vs IKKbWT;P301S+  DEGs
IKKbWT_P301S_EN_IN_DEGs <- FindMarkers(Cluster_EN_IN, ident.1 = "IKKbWT", ident.2 = "IKKbWT;P301S+", logfc.threshold = 0,
                                          test.use = "MAST")
write.csv(IKKbWT_P301S_EN_IN_DEGs, "IKKbWT_P301S_EN_IN_DEGs.csv")

IKKbWT_P301S_MG_DEGs <- FindMarkers(Cluster_MG, ident.1 = "IKKbWT", ident.2 = "IKKbWT;P301S+", logfc.threshold = 0,
                                       test.use = "MAST")
write.csv(IKKbWT_P301S_MG_DEGs, "IKKbWT_P301S_MG_DEGs.csv")

#IKKbWT vs IKKbCA;PS301S+  DEGs
IKKbWT_CAP301S_EN_IN_DEGs <- FindMarkers(Cluster_EN_IN, ident.1 = "IKKbWT", ident.2 = "IKKbCA;P301S+", logfc.threshold = 0,
                                            test.use = "MAST")
write.csv(IKKbWT_CAP301S_EN_IN_DEGs, "IKKbWT_CAP301S_EN_IN_DEGs.csv")

IKKbWT_CAP301S_MG_DEGs <- FindMarkers(Cluster_MG, ident.1 = "IKKbWT", ident.2 = "IKKbCA;P301S+", logfc.threshold = 0,
                                         test.use = "MAST")
write.csv(IKKbWT_CAP301S_MG_DEGs, "IKKbWT_CAP301S_MG_DEGs.csv")














