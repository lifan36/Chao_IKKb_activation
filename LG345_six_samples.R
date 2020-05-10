### example Seurat code for initial single cell RNA seq analysis ###
### written by Lay Kodama, Bang Liu, and Li Fan ###
### please reference https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html for details on the Seurat package parameters ###

#set working directory ====
setwd("/athena/ganlab/scratch/lif4001/ChaoWang_IKKb/IKKb_CA/LG345_82")

#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)

#load in data from Cell Ranger or other counts data ====

#for loading Cell Ranger counts:
LG345_82.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/ChaoWang_IKKb/cellranger/cellranger_count_LG345_82/outs/filtered_feature_bc_matrix")
LG345_82 <- CreateSeuratObject(counts = LG345_82.counts, project = "IKKbCA_1", min.cells = 3, min.features = 200)
LG345_82[["Condition"]] = c('IKKbCA')
rm(LG345_82.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
LG345_82[["percent.mt"]] <- PercentageFeatureSet(object = LG345_82, pattern = "^mt-") #recognize mitochondrial transcripts

all <- LG345_82
pdf("LG345_82_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG345_82_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
all
##An object of class Seurat 
##22387 features across 12315 samples within 1 assay 
##Active assay: RNA (22387 features)
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 8000 & nCount_RNA < 40000 & percent.mt < 5)
all
##An object of class Seurat 
##22092 features across 7500 samples within 1 assay 
##Active assay: RNA (22092 features)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG345_82_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("LG345_82_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"LG345_82_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("LG345_82_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG345_82_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.054*7873) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_425", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("LG345_82_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("LG345_82_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_425" #visualizing the singlet vs doublet cells
pdf("LG345_82_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

#Idents(object = all) <- "DF.classifications_0.25_0.005_552" #visualizing the singlet vs doublet cells
#pdf("LG345_82_3_UMAP_singlets_doublets_3.pdf", width=8, height=6)
#DimPlot(object = all, reduction = 'umap', label = F)
#dev.off()

saveRDS(all,"LG345_82_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"LG345_82_singlets.rds")

singlets<-readRDS("LG345_82_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("LG345_82_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("LG345_82_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("LG345_82_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("LG345_82_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("LG345_82_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("LG345_82_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"LG345_82_singlets_PCA.rds")

#annotate cell types based on outside databases - this step is tricky and arbitrary, but crucial to get as correct as possible:
#singlets <- RenameIdents(singlets, `0` = "exN", `1` = "Oligodendrocytes", `2` = "exN", 
#                         `3` = "ihN", `4` = "Astrocytes", `5` = "exN", `6` = "exN", `7` = "exN", '8' = "OPC", '9' = "Microglia", '10' = "Pericyte/EC")

#DimPlot(singlets, reduction = "tsne")

#load in meta data information for your genotypes etc - this step is very specific to your dataset:====
#Column_names_pbmc<-as.data.frame(singlets@meta.data$orig.ident)
#colnames(Column_names_pbmc) <- "Sample"
#sample_info <- read.csv("Sample_key.csv", header=T)
#library_info <- merge(Column_names_pbmc, sample_info, by="Sample")
#write.csv(library_info,"library_info.csv")

#singlets <- AddMetaData(object = singlets, metadata = library_info$Genotype, col.name = "Genotype")
#singlets <- AddMetaData(object =singlets, metadata = library_info$Sample, col.name = "Individual")
#saveRDS(singlets,"singlets_clustered.Rds")
#singlets<-readRDS("singlets_clustered.Rds")

#other functions you can use to analyze your data include:
#proportion of cells by genotype
#plotting tsnes by genotype
#calculating DEGs based on genotype for specific cell types
#subclustering specific cell types for finer clustering information
#please reference the different vignettes on Seurat's lab website: https://satijalab.org/seurat/vignettes.html

### example Seurat code for initial single cell RNA seq analysis ###
### written by Lay Kodama, Bang Liu, and Li Fan ###
### please reference https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html for details on the Seurat package parameters ###

#set working directory ====

#load in data from Cell Ranger or other counts data ====

#for loading Cell Ranger counts:
LG345_132.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/ChaoWang_IKKb/cellranger/cellranger_count_LG345_132/outs/filtered_feature_bc_matrix")
LG345_132 <- CreateSeuratObject(counts = LG345_132.counts, project = "IKKbCA_2", min.cells = 3, min.features = 200)
LG345_132[["Condition"]] = c('IKKbCA')
rm(LG345_132.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
LG345_132[["percent.mt"]] <- PercentageFeatureSet(object = LG345_132, pattern = "^mt-") #recognize mitochondrial transcripts

all <- LG345_132
pdf("LG345_132_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG345_132_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
all
##An object of class Seurat 
##22387 features across 12315 samples within 1 assay 
##Active assay: RNA (22387 features)
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 8000 & nCount_RNA < 40000 & percent.mt < 5)
all
##An object of class Seurat 
##22092 features across 7500 samples within 1 assay 
##Active assay: RNA (22092 features)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG345_132_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("LG345_132_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"LG345_132_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("LG345_132_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG345_132_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.091*12316) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_1121", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("LG345_132_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("LG345_132_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_1121" #visualizing the singlet vs doublet cells
pdf("LG345_132_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

#Idents(object = all) <- "DF.classifications_0.25_0.005_552" #visualizing the singlet vs doublet cells
#pdf("LG345_132_3_UMAP_singlets_doublets_3.pdf", width=8, height=6)
#DimPlot(object = all, reduction = 'umap', label = F)
#dev.off()

saveRDS(all,"LG345_132_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"LG345_132_singlets.rds")

singlets<-readRDS("LG345_132_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("LG345_132_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("LG345_132_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("LG345_132_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("LG345_132_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("LG345_132_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("LG345_132_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"LG345_132_singlets_PCA.rds")

#annotate cell types based on outside databases - this step is tricky and arbitrary, but crucial to get as correct as possible:
#singlets <- RenameIdents(singlets, `0` = "exN", `1` = "Oligodendrocytes", `2` = "exN", 
#                         `3` = "ihN", `4` = "Astrocytes", `5` = "exN", `6` = "exN", `7` = "exN", '8' = "OPC", '9' = "Microglia", '10' = "Pericyte/EC")

#DimPlot(singlets, reduction = "tsne")

#load in meta data information for your genotypes etc - this step is very specific to your dataset:====
#Column_names_pbmc<-as.data.frame(singlets@meta.data$orig.ident)
#colnames(Column_names_pbmc) <- "Sample"
#sample_info <- read.csv("Sample_key.csv", header=T)
#library_info <- merge(Column_names_pbmc, sample_info, by="Sample")
#write.csv(library_info,"library_info.csv")

#singlets <- AddMetaData(object = singlets, metadata = library_info$Genotype, col.name = "Genotype")
#singlets <- AddMetaData(object =singlets, metadata = library_info$Sample, col.name = "Individual")
#saveRDS(singlets,"singlets_clustered.Rds")
#singlets<-readRDS("singlets_clustered.Rds")

#other functions you can use to analyze your data include:
#proportion of cells by genotype
#plotting tsnes by genotype
#calculating DEGs based on genotype for specific cell types
#subclustering specific cell types for finer clustering information
#please reference the different vignettes on Seurat's lab website: https://satijalab.org/seurat/vignettes.html

### example Seurat code for initial single cell RNA seq analysis ###
### written by Lay Kodama, Bang Liu, and Li Fan ###
### please reference https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html for details on the Seurat package parameters ###


#load in data from Cell Ranger or other counts data ====

#for loading Cell Ranger counts:
LG345_98.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/ChaoWang_IKKb/cellranger/cellranger_count_LG345_98/outs/filtered_feature_bc_matrix")
LG345_98 <- CreateSeuratObject(counts = LG345_98.counts, project = "IKKbWT;P301S+_1", min.cells = 3, min.features = 200)
LG345_98[["Condition"]] = c('IKKbWT;P301S+')
rm(LG345_98.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
LG345_98[["percent.mt"]] <- PercentageFeatureSet(object = LG345_98, pattern = "^mt-") #recognize mitochondrial transcripts

all <- LG345_98
pdf("LG345_98_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG345_98_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
all
##An object of class Seurat 
##22387 features across 12315 samples within 1 assay 
##Active assay: RNA (22387 features)
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 8000 & nCount_RNA < 40000 & percent.mt < 5)
all
##An object of class Seurat 
##22092 features across 7500 samples within 1 assay 
##Active assay: RNA (22092 features)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG345_98_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("LG345_98_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"LG345_98_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("LG345_98_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG345_98_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.084*11145) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_936", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("LG345_98_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("LG345_98_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_936" #visualizing the singlet vs doublet cells
pdf("LG345_98_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

#Idents(object = all) <- "DF.classifications_0.25_0.005_552" #visualizing the singlet vs doublet cells
#pdf("LG345_98_3_UMAP_singlets_doublets_3.pdf", width=8, height=6)
#DimPlot(object = all, reduction = 'umap', label = F)
#dev.off()

saveRDS(all,"LG345_98_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"LG345_98_singlets.rds")

singlets<-readRDS("LG345_98_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("LG345_98_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("LG345_98_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("LG345_98_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("LG345_98_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("LG345_98_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("LG345_98_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"LG345_98_singlets_PCA.rds")

#annotate cell types based on outside databases - this step is tricky and arbitrary, but crucial to get as correct as possible:
#singlets <- RenameIdents(singlets, `0` = "exN", `1` = "Oligodendrocytes", `2` = "exN", 
#                         `3` = "ihN", `4` = "Astrocytes", `5` = "exN", `6` = "exN", `7` = "exN", '8' = "OPC", '9' = "Microglia", '10' = "Pericyte/EC")

#DimPlot(singlets, reduction = "tsne")

#load in meta data information for your genotypes etc - this step is very specific to your dataset:====
#Column_names_pbmc<-as.data.frame(singlets@meta.data$orig.ident)
#colnames(Column_names_pbmc) <- "Sample"
#sample_info <- read.csv("Sample_key.csv", header=T)
#library_info <- merge(Column_names_pbmc, sample_info, by="Sample")
#write.csv(library_info,"library_info.csv")

#singlets <- AddMetaData(object = singlets, metadata = library_info$Genotype, col.name = "Genotype")
#singlets <- AddMetaData(object =singlets, metadata = library_info$Sample, col.name = "Individual")
#saveRDS(singlets,"singlets_clustered.Rds")
#singlets<-readRDS("singlets_clustered.Rds")

#other functions you can use to analyze your data include:
#proportion of cells by genotype
#plotting tsnes by genotype
#calculating DEGs based on genotype for specific cell types
#subclustering specific cell types for finer clustering information
#please reference the different vignettes on Seurat's lab website: https://satijalab.org/seurat/vignettes.html

### example Seurat code for initial single cell RNA seq analysis ###
### written by Lay Kodama, Bang Liu, and Li Fan ###
### please reference https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html for details on the Seurat package parameters ###


#load in data from Cell Ranger or other counts data ====

#for loading Cell Ranger counts:
LG345_135.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/ChaoWang_IKKb/cellranger/cellranger_count_LG345_135/outs/filtered_feature_bc_matrix")
LG345_135 <- CreateSeuratObject(counts = LG345_135.counts, project = "IKKbWT;P301S+_2", min.cells = 3, min.features = 200)
LG345_135[["Condition"]] = c('IKKbWT;P301S+')
rm(LG345_135.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
LG345_135[["percent.mt"]] <- PercentageFeatureSet(object = LG345_135, pattern = "^mt-") #recognize mitochondrial transcripts

all <- LG345_135
pdf("LG345_135_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG345_135_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
all
##An object of class Seurat 
##22387 features across 12315 samples within 1 assay 
##Active assay: RNA (22387 features)
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 8000 & nCount_RNA < 40000 & percent.mt < 5)
all
##An object of class Seurat 
##22092 features across 7500 samples within 1 assay 
##Active assay: RNA (22092 features)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG345_135_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("LG345_135_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"LG345_135_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("LG345_135_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG345_135_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.084*11401) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_958", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("LG345_135_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("LG345_135_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_958" #visualizing the singlet vs doublet cells
pdf("LG345_135_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

#Idents(object = all) <- "DF.classifications_0.25_0.005_552" #visualizing the singlet vs doublet cells
#pdf("LG345_135_3_UMAP_singlets_doublets_3.pdf", width=8, height=6)
#DimPlot(object = all, reduction = 'umap', label = F)
#dev.off()

saveRDS(all,"LG345_135_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"LG345_135_singlets.rds")

singlets<-readRDS("LG345_135_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("LG345_135_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("LG345_135_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("LG345_135_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("LG345_135_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("LG345_135_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("LG345_135_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"LG345_135_singlets_PCA.rds")

#annotate cell types based on outside databases - this step is tricky and arbitrary, but crucial to get as correct as possible:
#singlets <- RenameIdents(singlets, `0` = "exN", `1` = "Oligodendrocytes", `2` = "exN", 
#                         `3` = "ihN", `4` = "Astrocytes", `5` = "exN", `6` = "exN", `7` = "exN", '8' = "OPC", '9' = "Microglia", '10' = "Pericyte/EC")

#DimPlot(singlets, reduction = "tsne")

#load in meta data information for your genotypes etc - this step is very specific to your dataset:====
#Column_names_pbmc<-as.data.frame(singlets@meta.data$orig.ident)
#colnames(Column_names_pbmc) <- "Sample"
#sample_info <- read.csv("Sample_key.csv", header=T)
#library_info <- merge(Column_names_pbmc, sample_info, by="Sample")
#write.csv(library_info,"library_info.csv")

#singlets <- AddMetaData(object = singlets, metadata = library_info$Genotype, col.name = "Genotype")
#singlets <- AddMetaData(object =singlets, metadata = library_info$Sample, col.name = "Individual")
#saveRDS(singlets,"singlets_clustered.Rds")
#singlets<-readRDS("singlets_clustered.Rds")

#other functions you can use to analyze your data include:
#proportion of cells by genotype
#plotting tsnes by genotype
#calculating DEGs based on genotype for specific cell types
#subclustering specific cell types for finer clustering information
#please reference the different vignettes on Seurat's lab website: https://satijalab.org/seurat/vignettes.html

### example Seurat code for initial single cell RNA seq analysis ###
### written by Lay Kodama, Bang Liu, and Li Fan ###
### please reference https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html for details on the Seurat package parameters ###



#load in data from Cell Ranger or other counts data ====

#for loading Cell Ranger counts:
LG345_84.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/ChaoWang_IKKb/cellranger/cellranger_count_LG345_84/outs/filtered_feature_bc_matrix")
LG345_84 <- CreateSeuratObject(counts = LG345_84.counts, project = "IKKbCA;P301S+_1", min.cells = 3, min.features = 200)
LG345_84[["Condition"]] = c('IKKbCA;P301S+')
rm(LG345_84.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
LG345_84[["percent.mt"]] <- PercentageFeatureSet(object = LG345_84, pattern = "^mt-") #recognize mitochondrial transcripts

all <- LG345_84
pdf("LG345_84_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG345_84_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
all
##An object of class Seurat 
##22387 features across 12315 samples within 1 assay 
##Active assay: RNA (22387 features)
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 8000 & nCount_RNA < 40000 & percent.mt < 5)
all
##An object of class Seurat 
##22092 features across 7500 samples within 1 assay 
##Active assay: RNA (22092 features)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG345_84_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("LG345_84_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"LG345_84_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("LG345_84_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG345_84_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.076*10976) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_834", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("LG345_84_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("LG345_84_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_834" #visualizing the singlet vs doublet cells
pdf("LG345_84_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

#Idents(object = all) <- "DF.classifications_0.25_0.005_552" #visualizing the singlet vs doublet cells
#pdf("LG345_84_3_UMAP_singlets_doublets_3.pdf", width=8, height=6)
#DimPlot(object = all, reduction = 'umap', label = F)
#dev.off()

saveRDS(all,"LG345_84_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"LG345_84_singlets.rds")

singlets<-readRDS("LG345_84_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("LG345_84_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("LG345_84_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("LG345_84_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("LG345_84_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("LG345_84_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("LG345_84_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"LG345_84_singlets_PCA.rds")

#annotate cell types based on outside databases - this step is tricky and arbitrary, but crucial to get as correct as possible:
#singlets <- RenameIdents(singlets, `0` = "exN", `1` = "Oligodendrocytes", `2` = "exN", 
#                         `3` = "ihN", `4` = "Astrocytes", `5` = "exN", `6` = "exN", `7` = "exN", '8' = "OPC", '9' = "Microglia", '10' = "Pericyte/EC")

#DimPlot(singlets, reduction = "tsne")

#load in meta data information for your genotypes etc - this step is very specific to your dataset:====
#Column_names_pbmc<-as.data.frame(singlets@meta.data$orig.ident)
#colnames(Column_names_pbmc) <- "Sample"
#sample_info <- read.csv("Sample_key.csv", header=T)
#library_info <- merge(Column_names_pbmc, sample_info, by="Sample")
#write.csv(library_info,"library_info.csv")

#singlets <- AddMetaData(object = singlets, metadata = library_info$Genotype, col.name = "Genotype")
#singlets <- AddMetaData(object =singlets, metadata = library_info$Sample, col.name = "Individual")
#saveRDS(singlets,"singlets_clustered.Rds")
#singlets<-readRDS("singlets_clustered.Rds")

#other functions you can use to analyze your data include:
#proportion of cells by genotype
#plotting tsnes by genotype
#calculating DEGs based on genotype for specific cell types
#subclustering specific cell types for finer clustering information
#please reference the different vignettes on Seurat's lab website: https://satijalab.org/seurat/vignettes.html

### example Seurat code for initial single cell RNA seq analysis ###
### written by Lay Kodama, Bang Liu, and Li Fan ###
### please reference https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html for details on the Seurat package parameters ###



#load in data from Cell Ranger or other counts data ====

#for loading Cell Ranger counts:
LG345_133.counts <- Read10X(data.dir = "/athena/ganlab/scratch/lif4001/ChaoWang_IKKb/cellranger/cellranger_count_LG345_133/outs/filtered_feature_bc_matrix")
LG345_133 <- CreateSeuratObject(counts = LG345_133.counts, project = "IKKbCA;P301S+_2", min.cells = 3, min.features = 200)
LG345_133[["Condition"]] = c('IKKbCA;P301S+')
rm(LG345_133.counts)
#vizualize QC metrics and filtering====
#mitochondrial transcripts - if the cell has high mitochondrial transcripts, it may signal a cell under stress/unhealthy
LG345_133[["percent.mt"]] <- PercentageFeatureSet(object = LG345_133, pattern = "^mt-") #recognize mitochondrial transcripts

all <- LG345_133
pdf("LG345_133_QC.pdf", width=12, height=4)
VlnPlot(object = all, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size=0)
dev.off()

#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("LG345_133_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
all
##An object of class Seurat 
##22387 features across 12315 samples within 1 assay 
##Active assay: RNA (22387 features)
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 8000 & nCount_RNA < 40000 & percent.mt < 5)
all
##An object of class Seurat 
##22092 features across 7500 samples within 1 assay 
##Active assay: RNA (22092 features)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
all <- ScaleData(object = all)
#perform and visualize PCA
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 1000)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("LG345_133_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)

pdf("LG345_133_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()

saveRDS(all,"LG345_133_PCA.rds") #it's good to save your R object periodically so you can start from this object without having to go through the processing steps again.
all<-readRDS("LG345_133_PCA.rds")

#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("LG345_133_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()

length(all@meta.data$seurat_clusters)
#homotypic doublet proportion estimate
annotations <- all@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.084*11293) #estimate the number of multiplets you expect from the kit you are using - should give you percent expected based on number of nuclei inputs
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#doublet finder with different classification stringencies
all <- doubletFinder_v3(all, PCs=1:15, pN=0.25, pK=0.005, nExp=nExp_poi, reuse.pANN = FALSE,sct=FALSE)
all <- doubletFinder_v3(all, PCs = 1:15, pN = 0.25, pK = 0.005, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.005_949", sct = FALSE)
#visualizing clusters and multiplet cells====
pdf("LG345_133_Elbow_2.pdf", width=8, height=6)
ElbowPlot(all,ndims=50) #only run clustering on PCs that explain the variations, ie at the bend of the elbow (not the entire dataset, otherwise will take too long)
dev.off()
all <- FindNeighbors(object = all, dims = 1:15)
all <- FindClusters(object = all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("LG345_133_UMAP_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = T)
dev.off()

Idents(object = all) <- "DF.classifications_0.25_0.005_949" #visualizing the singlet vs doublet cells
pdf("LG345_133_3_UMAP_singlets_doublets_2.pdf", width=8, height=6)
DimPlot(object = all, reduction = 'umap', label = F)
dev.off()

#Idents(object = all) <- "DF.classifications_0.25_0.005_552" #visualizing the singlet vs doublet cells
#pdf("LG345_133_3_UMAP_singlets_doublets_3.pdf", width=8, height=6)
#DimPlot(object = all, reduction = 'umap', label = F)
#dev.off()

saveRDS(all,"LG345_133_after_doublet_detection.rds")

#all@active.ident

#processing singlets ====
#remove doublets
singlets <- subset(all, idents=c("Singlet"))
rm(all)

saveRDS(singlets,"LG345_133_singlets.rds")

singlets<-readRDS("LG345_133_singlets.rds")

Idents(singlets) <- "seurat_clusters"
pdf("LG345_133_UMAP_singlets_before_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap',label = T)
dev.off()

pdf("LG345_133_Elbow_middle.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("LG345_133_UMAP_singlets_middle.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

#Reload singlets
singlets<-readRDS("LG345_133_singlets.rds")
Idents(singlets) <- "seurat_clusters"

#normalization
singlets <- NormalizeData(singlets, normalization.method = "LogNormalize", scale.factor = 10000)
#scaling the data
singlets <- ScaleData(object = singlets)
#perform and visualize PCA
singlets <- RunPCA(object = singlets, features = rownames(x = singlets), verbose = FALSE)

#PC capture
pdf("LG345_133_Elbow_after_processing.pdf", width=8, height=6)
ElbowPlot(singlets,ndims=50)
dev.off()
singlets <- FindNeighbors(object = singlets, dims = 1:15)
singlets <- FindClusters(object = singlets, resolution = 0.1)
singlets <- RunUMAP(object = singlets, dims = 1:15)
pdf("LG345_133_UMAP_singlets_after_processing.pdf", width=8, height=6)
DimPlot(object = singlets, reduction = 'umap', label = T)
dev.off()

saveRDS(singlets,"LG345_133_singlets_PCA.rds")

#annotate cell types based on outside databases - this step is tricky and arbitrary, but crucial to get as correct as possible:
#singlets <- RenameIdents(singlets, `0` = "exN", `1` = "Oligodendrocytes", `2` = "exN", 
#                         `3` = "ihN", `4` = "Astrocytes", `5` = "exN", `6` = "exN", `7` = "exN", '8' = "OPC", '9' = "Microglia", '10' = "Pericyte/EC")

#DimPlot(singlets, reduction = "tsne")

#load in meta data information for your genotypes etc - this step is very specific to your dataset:====
#Column_names_pbmc<-as.data.frame(singlets@meta.data$orig.ident)
#colnames(Column_names_pbmc) <- "Sample"
#sample_info <- read.csv("Sample_key.csv", header=T)
#library_info <- merge(Column_names_pbmc, sample_info, by="Sample")
#write.csv(library_info,"library_info.csv")

#singlets <- AddMetaData(object = singlets, metadata = library_info$Genotype, col.name = "Genotype")
#singlets <- AddMetaData(object =singlets, metadata = library_info$Sample, col.name = "Individual")
#saveRDS(singlets,"singlets_clustered.Rds")
#singlets<-readRDS("singlets_clustered.Rds")

#other functions you can use to analyze your data include:
#proportion of cells by genotype
#plotting tsnes by genotype
#calculating DEGs based on genotype for specific cell types
#subclustering specific cell types for finer clustering information
#please reference the different vignettes on Seurat's lab website: https://satijalab.org/seurat/vignettes.html