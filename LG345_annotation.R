library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(MAST)

setwd("/athena/ganlab/scratch/lif4001/ChaoWang_IKKb/IKKb_CA/LG345_singlets")

LG345_integrated <- readRDS("LG345_integrated_PCA_0.1.rds")

DefaultAssay(LG345_integrated) <- 'RNA'


LG345_integrated <- RenameIdents(LG345_integrated,
                                        `0` = "oligodendrocytes", `1`="excitatory neurons", `2`="excitatory neurons", `3`="inhibitory neurons",
                                        `4`="astrocytes", `5`="excitatory neurons", `6`="mixed cell type", `7`="microglia",
                                        `8`="inhibitory neurons", `9`="inhibitory neurons", `10`="OPCs", `11`="inhibitory neurons",
                                        `12`="inhibitory neurons", `13`="excitatory neurons", `14`="excitatory neurons", `15`="vascular cells",
                                        `16`="excitatory neurons", `17`="excitatory neurons", `18`="endothelial cells", `19`="vascular cells",
                                        `20`="excitatory neurons", `21`="vascular cells" , `22`="ependymal cells"
)

pdf("LG345_integrated_umap_annotation.pdf", width=12, height=8)
DimPlot(LG345_integrated, reduction = 'umap', label = TRUE)
dev.off()
pdf("LG345_integrated_umap_annotation_split_Condition.pdf", width=16, height=4)
DimPlot(LG345_integrated, reduction = 'umap', split.by = "Condition", label = TRUE)
dev.off()

Idents(LG345_integrated) <- factor(Idents(LG345_integrated), levels = c(names(table(LG345_integrated@active.ident))))
markers.to.plot <- c("Plp1", "Mbp", "Mobp","Slc17a7", "Nrgn", "Gad1", "Gad2", "Clu", "Aldoc", 
                     "Pla2g7", "Cx3cr1", "P2ry12", "Csf1r","Scrg1", "Pdgfra", "Vtn", "Igfbp7",
                     "Bnc2", "Slc47a1", "Ttr")
pdf("LG345_integrated_DotPlot_top2-3.pdf", width=12, height=10)
DotPlot(LG345_integrated, features = rev(markers.to.plot), cols = c("orange", "plum", "pink", "grey"), dot.scale = 8, 
        split.by = "Condition") + RotatedAxis()
dev.off()

saveRDS(LG345_integrated, file = "LG345_integrated_ready_4_DEGs.rds")


