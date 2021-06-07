#Script for Analysis of 2 groups  with 3 biological replicates each
library(dplyr)
library(Seurat)
library(Matrix)
library(cowplot)

# Set up control objects
ctrl1 <- Read10X(data.dir = 'filtered_feature_bc_matrix_WT151_premRNA/')
ctrl1 <- CreateSeuratObject(counts = ctrl1, min.cells = 3, min.features = 200, project = "WT_151")
ctrl1$KO <- "CTRL"
ctrl2 <- Read10X(data.dir = 'filtered_feature_bc_matrix_WT152_premRNA/')
ctrl2 <- CreateSeuratObject(counts = ctrl2, min.cells = 3, min.features = 200, project = "WT_152")
ctrl2$KO <- "CTRL"
ctrl3 <- Read10X(data.dir = 'filtered_feature_bc_matrix_WT156_premRNA/')
ctrl3 <- CreateSeuratObject(counts = ctrl3, min.cells = 3, min.features = 200, project = "WT_156")
ctrl3$KO <- "CTRL"

# Set up KO objects
KO1 <- Read10X(data.dir = 'filtered_feature_bc_matrix_KO161_premRNA/')
KO1 <- CreateSeuratObject(counts = KO1, min.cells = 3, min.features = 200, project = "KO_161")
KO1$KO <- "KO"
KO2 <- Read10X(data.dir = 'filtered_feature_bc_matrix_KO691_premRNA/')
KO2 <- CreateSeuratObject(counts = KO2, min.cells = 3, min.features = 200, project = "KO_691")
KO2$KO <- "KO"
KO3 <- Read10X(data.dir = 'filtered_feature_bc_matrix_KO692_premRNA/')
KO3 <- CreateSeuratObject(counts = KO3, min.cells = 3, min.features = 200, project = "KO_692")
KO3$KO <- "KO"

#Merge Controls into a seurat object
CtrlCombined <- merge(ctrl1, y = c(ctrl2, ctrl3), add.cell.ids = c("WT_151", "WT_152", "WT_156"), project = "Combined")
Ctrl <- CtrlCombined

#Merge Klf4Kos into a seurat object
KOCombined <- merge(KO1, y = c(KO2, Ko3), add.cell.ids = c("KO_161", "KO_691", "KO_692"), project = "Combined")
KO <- KOCombined

#Merge Controls and KOs into a seurat object
DataCombined <- merge(Ctrl, y = KO, add.cell.ids = c("Ctrl", "KO"), project = "Combined")

#Assess mt content
DataCombined[["percent.mt"]] <- PercentageFeatureSet(object = DataCombined, pattern = "^mt-")
VlnPlot(Data.Combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

# Filter cells, only keep cells with more than 200 and less than 5000 features, and less than 10% MT-genes
DataCombined <- subset(x= DataCombined, subset= nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt <10)
DataCombined <- NormalizeData(object = DataCombined, verbose = FALSE)
DataCombined <- FindVariableFeatures(object = DataCombined, selection.method = "vst", nfeatures = 2000)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plot1 <- FeatureScatter(DataCombined, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(DataCombined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))

#Identify Variable Features 
obj.list <- SplitObject(DataCombined, split.by = "orig.ident")
obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Perform Integration and Identify Anchors
Data.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:30)
Data.Combined <- IntegrateData(anchorset = Data.anchors, dims = 1:30)

# Perform an integrated analysis
DefaultAssay(object = Data.Combined) <- "integrated"

# Run the standard workflow for visualization and clustering
Data.Combined <- ScaleData(object = Data.Combined, verbose = FALSE)
Data.Combined <- RunPCA(object = Data.Combined, npcs = 30, verbose = FALSE) 

# t-SNE and Clustering
Data.Combined <- RunUMAP(object = Data.Combined, reduction = "pca", dims = 1:30) 
Data.Combined <- FindNeighbors(object = Data.Combined, reduction = "pca", dims = 1:30)
Data.Combined <- FindClusters(Data.Combined, resolution = 0.5)

# Visualization of UMAP
DimPlot(Data.Combined, reduction = "umap", group.by = "orig.ident")
DimPlot(Data.Combined, reduction = "umap", group.by = "KO")
DimPlot(Data.Combined, reduction = "umap", label = TRUE)

#show WT and KO side by side
DimPlot(Data.Combined, reduction = "umap", split.by = "orig.ident", ncol = 2, label = FALSE)
DimPlot(Data.Combined, reduction = "umap", split.by = "KO", label = TRUE)

###############
## Save File ##
###############
saveRDS(Data.Combined, file = "6Samples2.rds")


#####################
## Load Saved File ##
#####################
Data.Combined <- readRDS("6Samples_subclustered.rds")

#To output average expression split by group
cluster.averages <- AverageExpression(Data.Combined, assays = "RNA", add.ident = "KO")
write.csv(cluster.averages[["RNA"]], file = "6samples_ClusterAverages_Updated.csv")

# Get number of cells per cluster and per sample/group of origin
table(Data.Combined@meta.data[["seurat_clusters"]], Data.Combined@meta.data$orig.ident)
table(Data.Combined@meta.data[["seurat_clusters"]], Data.Combined@meta.data$KO)

#Feature and Violin Plots
DefaultAssay(object = Data.Combined) <- "RNA"
VlnPlot(Data.Combined, features = c("Nphs1", "Nphs2", "Synpo", "Wt1"), split.by = "KO", ncol = 1, pt.size = 0)
VlnPlot(Data.Combined, features = c("Akap12"), split.by = "KO", pt.size = 1)
FeaturePlot(Data.Combined, features = c("Nphs1"), ncol = 2, split.by = "orig.ident")
FeaturePlot(Data.Combined, features = c("Cd44"), ncol = 2, split.by = "KO", pt.size = 0.01)
DotPlot(Data.Combined, features = c("Nphs1", "Nphs2"), cols = c("yellow", "red"), col.min = 0, dot.scale = 8) + RotatedAxis()

#For Dot plot (combined)
markers.to.plot <- c("Nphs1", "Emcn", "Csmd1", "Cd44", "Rad51b", "Slc5a12", "Slc7a7", "Slc7a13", "Slc12a1", "Slc12a3", "Slc8a1", "Egfem1", "Aqp2", "Kit", "Slc26a4", "Ctnna3", "Akap12", "Ptprc")
pdf("DotPlot.pdf", width=8, height=6)
DotPlot(Data.Combined, features = rev(markers.to.plot), cols = c("yellow", "red"), col.min = 0, dot.scale = 8) + RotatedAxis()
dev.off()


#Assigning new cluster names
new.cluster.ids <- c( "PT(S3)-1", "PT(S1-S2)-1", "PT(S3)-2", "PEC/Prolif.PT", "Mes/Fib", "CNT/DCT", "PT(S1-S2)-2", "LoH-1", "LoH-2", "CNT", "Novel1", "DCT", "Endo", "IC-B", "Prolif.PT", "Pod", "IC-A", "CD-PC", "Mø", "Fib", "Novel2", "JGA", "LoH-3")
names(new.cluster.ids) <- levels(Data.Combined)
Data.Combined <- RenameIdents(Data.Combined, new.cluster.ids)
DimPlot(Data.Combined, reduction = "umap", label = TRUE, pt.size = 0.5)

#Reorder new cluster names
levels(x = Data.Combined) <- c( "Novel2", "Novel1", "Mø", "Fib", "JGA", "IC-B", "IC-A", "CD-PC", "CNT", "CNT/DCT", "DCT", "LoH-3", "LoH-2", "LoH-1", "PT(S3)-2", "PT(S3)-1", "PT(S1-S2)-2", "PT(S1-S2)-1", "Prolif.PT", "PEC/Prolif.PT", "Mes/Fib", "Endo", "Pod")
new.cluster.ids <- c( "Novel2", "Novel1", "Mø", "Fib", "JGA", "IC-B", "IC-A", "CD-PC", "CNT", "CNT/DCT", "DCT", "LoH-3", "LoH-2", "LoH-1", "PT(S3)-2", "PT(S3)-1", "PT(S1-S2)-2", "PT(S1-S2)-1", "Prolif.PT", "PEC/Prolif.PT", "Mes/Fib", "Endo", "Pod")
#use this order for everything that isn't Dot Plot
levels(x = Data.Combined) <- rev(c( "Novel2", "Novel1", "Mø", "Fib", "JGA", "IC-B", "IC-A", "CD-PC", "CNT", "CNT/DCT", "DCT", "LoH-3", "LoH-2", "LoH-1", "PT(S3)-2", "PT(S3)-1", "PT(S1-S2)-2", "PT(S1-S2)-1", "Prolif.PT", "PEC/Prolif.PT", "Mes/Fib", "Endo", "Pod"))

#Receptors on Dot Plot
DefaultAssay(Data.Combined) <- "RNA"
receptors.to.plot <- c("Cd44", "Cd74", "Traf2", "Bsg", "Itga5", "Itgam", "Itgb2", "Ddr1", "Itga1", "Itga10", "Itga2", "Dip2a", "Bmp4", "Pappa", "Fzd8", "Igf1R", "Acvr1", "Eng", "Tgfbr1", "Tgfbr2", "Tgfbr3", "Cd47", "Cd93", "Itgav", "Itgb8", "Itgb1", "Fzd1", "Acvrl1", "Cav1", "Cd109", "Cxcr4", "Itgb6", "Sdc2", "Itgb3", "Acvr1b", "Flt4", "Il17cr", "Itga6", "Itga9", "Itgb7", "Plaur", "Ddr2", "Lrp4", "Sdc1")
DotPlot(Data.Combined, features = rev(receptors.to.plot), cols = c("yellow", "red"), col.min = 0, dot.scale = 8, split.by = "KO") + RotatedAxis()
p = DotPlot(Data.Combined, features = rev(receptors.to.plot), cols = c("yellow", "red"), col.min = 0, dot.scale = 8) 
p + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(axis.text.y = element_text(angle = 150, hjust = 0, vjust=0.5))

levels(x = immune.combined) <- c( "0.CTRL", "0.KO", "1.CTRL", "1.KO", "2.CTRL", "2.KO", "3.CTRL", "3.KO", "4.CTRL", "4.KO", "5.CTRL", "5.KO", "6.CTRL", "6.KO", "7.CTRL", "7.KO", "8.CTRL", "8.KO")
DefaultAssay(immune.combined) <- "RNA"
receptors.to.plot <- c("Lrp4", "Cd44", "Cd93", "Ddr1", "Ddr2", "Flt4", "Itga1", "Itga2", "Itga5", "Itgav", "Itga10", "Cd47", "Itgb8", "Itgam", "Itgb2", "Itgb3", "Itgb6", "Il17rc", "Itga6", "Itga9", "Itgb7", "Plaur", "Traf2", "Fzd8", "Acvr1b", "Tmem222", "SDC1", "Cd74")
p = DotPlot(immune.combined, features = rev(receptors.to.plot), cols = c("lightgray", "limegreen"), col.min = 0, dot.scale = 8) 
p + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(axis.text.y = element_text(angle = 150, hjust = 0, vjust=0.5))


DefaultAssay(Data.Combined) <- "RNA"
receptors.to.plot <- c("Lrp4", "Cd44", "Cd93", "Ddr1", "Ddr2", "Flt4", "Itga1", "Itga2", "Itga5", "Itgav", "Itga10", "Cd47", "Itgb8", "Itgam", "Itgb2", "Itgb3", "Itgb6", "Il17rc", "Itga6", "Itga9", "Itgb7", "Plaur", "Traf2", "Fzd8", "Acvr1b", "Tmem222", "SDC1", "Cd74")
p = DotPlot(Data.Combined, features = rev(receptors.to.plot), cols = c("lightgray", "limegreen"), col.min = 0, dot.scale = 8) 
p + theme(axis.text.x = element_text(angle = 90, hjust = 0.5)) + theme(axis.text.y = element_text(angle = 150, hjust = 0, vjust=0.5))

#Ligands on Dot Plot
DefaultAssay(Data.Combined) <- "RNA"
ligands.to.plot <- c("Agrn", "Anxa5", "Bmp1", "Cnn2", "Col1a1", "Col2a1", "Col4a1", "Col4a2", "Col5a2", "Cpe", "Ctgf", "Fbn1", "Fn1", "Fstl1", "Gfra1", "Gstp1", "Hapln1", "Igfbp4", "Igfbp5", "Igfbp7", "Inhba", "Lgals1", "Loxl2", "Ltbp2", "Mif", "Nucb1", "P4hb", "Ppia", "Prss23", "Sparc", "Stc2", "Tgfb2", "Tgfbi", "Thbs1", "Tkt")
DotPlot(Data.Combined, features = rev(ligands.to.plot), cols = c("yellow", "red"), col.min = 0, dot.scale = 8, split.by = "KO") + RotatedAxis()

FeatureScatter(Data.Combined, feature1 = "Cd44", feature2 = "Itgav")

#For Differential Expression between Control and KO of the same cluster (ie Cluster6_CTRL vs Cluster6_KO)
immune.combined <- Data.Combined
immune.combined$celltype.KO <- paste(Idents(immune.combined), immune.combined$KO, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.KO"

#Cluster 0
C0 <- FindMarkers(immune.combined, ident.1 = "0_KO", ident.2 = "0_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C0, n = 15)
#Cluster 1
C1 <- FindMarkers(immune.combined, ident.1 = "1_KO", ident.2 = "1_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C1, n = 15)
#Cluster 2
C2 <- FindMarkers(immune.combined, ident.1 = "2_KO", ident.2 = "2_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C2, n = 15)
#Cluster 3
C3 <- FindMarkers(immune.combined, ident.1 = "3_KO", ident.2 = "3_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C3, n = 15)
#Cluster 4
C4 <- FindMarkers(immune.combined, ident.1 = "4_KO", ident.2 = "4_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C4, n = 15)
#Cluster 5
C5 <- FindMarkers(immune.combined, ident.1 = "5_KO", ident.2 = "5_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C5, n = 15)
#Cluster 6
C6 <- FindMarkers(immune.combined, ident.1 = "6_KO", ident.2 = "6_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C6, n = 15)
#Cluster 7
C7 <- FindMarkers(immune.combined, ident.1 = "7_KO", ident.2 = "7_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C7, n = 15)
#Cluster 8
C8 <- FindMarkers(immune.combined, ident.1 = "8_KO", ident.2 = "8_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C8, n = 15)
#Cluster 9
C9 <- FindMarkers(immune.combined, ident.1 = "9_KO", ident.2 = "9_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C9, n = 15)
#Cluster 10
C10 <- FindMarkers(immune.combined, ident.1 = "10_KO", ident.2 = "10_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C10, n = 15)
#Cluster 11
C11 <- FindMarkers(immune.combined, ident.1 = "11_KO", ident.2 = "11_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C11, n = 15)
#Cluster 12
C12 <- FindMarkers(immune.combined, ident.1 = "12_KO", ident.2 = "12_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C12, n = 15)
#Cluster 13
C13 <- FindMarkers(immune.combined, ident.1 = "13_KO", ident.2 = "13_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C13, n = 15)
#Cluster 14
C14 <- FindMarkers(immune.combined, ident.1 = "14_KO", ident.2 = "14_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C14, n = 15)
#Cluster 15
C15 <- FindMarkers(immune.combined, ident.1 = "15_KO", ident.2 = "15_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C15, n = 15)
#Cluster 16
C16 <- FindMarkers(immune.combined, ident.1 = "16_KO", ident.2 = "16_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C16, n = 15)
#Cluster 17
C17 <- FindMarkers(immune.combined, ident.1 = "17_KO", ident.2 = "17_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C17, n = 15)
#Cluster 18
C18 <- FindMarkers(immune.combined, ident.1 = "18_KO", ident.2 = "18_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C18, n = 15)
#Cluster 19
C19 <- FindMarkers(immune.combined, ident.1 = "19_KO", ident.2 = "19_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C19, n = 15)
#Cluster 20
C20 <- FindMarkers(immune.combined, ident.1 = "20_KO", ident.2 = "20_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C20, n = 15)
#Cluster 21
C21 <- FindMarkers(immune.combined, ident.1 = "21_KO", ident.2 = "21_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C21, n = 15)
#Cluster 22
C22 <- FindMarkers(immune.combined, ident.1 = "22_KO", ident.2 = "22_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
head(C22, n = 15)

#add row names (Genes) as its own column 
C0 <- data.frame(Gene = row.names(C0), C0)
C1 <- data.frame(Gene = row.names(C1), C1)
C2 <- data.frame(Gene = row.names(C2), C2)
C3 <- data.frame(Gene = row.names(C3), C3)
C4 <- data.frame(Gene = row.names(C4), C4)
C5 <- data.frame(Gene = row.names(C5), C5)
C6 <- data.frame(Gene = row.names(C6), C6)
C7 <- data.frame(Gene = row.names(C7), C7)
C8 <- data.frame(Gene = row.names(C8), C8)
C9 <- data.frame(Gene = row.names(C9), C9)
C10 <- data.frame(Gene = row.names(C10), C10)
C11 <- data.frame(Gene = row.names(C11), C11)
C12 <- data.frame(Gene = row.names(C12), C12)
C13 <- data.frame(Gene = row.names(C13), C13)
C14 <- data.frame(Gene = row.names(C14), C14)
C15 <- data.frame(Gene = row.names(C15), C15)
C16 <- data.frame(Gene = row.names(C16), C16)
C17 <- data.frame(Gene = row.names(C17), C17)
C18 <- data.frame(Gene = row.names(C18), C18)
C19 <- data.frame(Gene = row.names(C19), C19)
C20 <- data.frame(Gene = row.names(C20), C20)
C21 <- data.frame(Gene = row.names(C21), C21)
C22 <- data.frame(Gene = row.names(C22), C22)

#Save
require(openxlsx)
list_of_datasets <- list("C0" = C0, "C1" = C1, "C2" = C2, "C3" = C3, "C4" = C4, "C5" = C5, "C6" = C6, "C7" = C7, "C8" = C8, "C9" = C9, "C10" = C10, "C11" = C11, "C12" = C12, "C13" = C13, "C14" = C14, "C15" = C15, "C16" = C16, "C17" = C17, "C18" = C18, "C19" = C19, "C20" =C20, "C21" =C21, "C22" = C22)
write.xlsx(list_of_datasets, file = "6Samp_DiffExprGenes_CTRLvsKO_Updated2.xlsx")

#For Differential Expression between a single cluster and all other cells (Cluster 1 vs all other clusters)
immune.combined <- Data.Combined
immune.combined$celltype.KO <- paste(Idents(immune.combined), immune.combined$KO, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype"

C0 <- FindMarkers(immune.combined, ident.1 = "0", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C1 <- FindMarkers(immune.combined, ident.1 = "1", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C2 <- FindMarkers(immune.combined, ident.1 = "2", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C3 <- FindMarkers(immune.combined, ident.1 = "3", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C4 <- FindMarkers(immune.combined, ident.1 = "4", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C5 <- FindMarkers(immune.combined, ident.1 = "5", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C6 <- FindMarkers(immune.combined, ident.1 = "6", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C7 <- FindMarkers(immune.combined, ident.1 = "7", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C8 <- FindMarkers(immune.combined, ident.1 = "8", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C9 <- FindMarkers(immune.combined, ident.1 = "9", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C10 <- FindMarkers(immune.combined, ident.1 = "10", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C11 <- FindMarkers(immune.combined, ident.1 = "11", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C12 <- FindMarkers(immune.combined, ident.1 = "12", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C13 <- FindMarkers(immune.combined, ident.1 = "13", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C14 <- FindMarkers(immune.combined, ident.1 = "14", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C15 <- FindMarkers(immune.combined, ident.1 = "15", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C16 <- FindMarkers(immune.combined, ident.1 = "16", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C17 <- FindMarkers(immune.combined, ident.1 = "17", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C18 <- FindMarkers(immune.combined, ident.1 = "18", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C19 <- FindMarkers(immune.combined, ident.1 = "19", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C20 <- FindMarkers(immune.combined, ident.1 = "20", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C21 <- FindMarkers(immune.combined, ident.1 = "21", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C22 <- FindMarkers(immune.combined, ident.1 = "22", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)


#add row names (Genes) as its own column 
C0 <- data.frame(Gene = row.names(C0), C0)
C1 <- data.frame(Gene = row.names(C1), C1)
C2 <- data.frame(Gene = row.names(C2), C2)
C3 <- data.frame(Gene = row.names(C3), C3)
C4 <- data.frame(Gene = row.names(C4), C4)
C5 <- data.frame(Gene = row.names(C5), C5)
C6 <- data.frame(Gene = row.names(C6), C6)
C7 <- data.frame(Gene = row.names(C7), C7)
C8 <- data.frame(Gene = row.names(C8), C8)
C9 <- data.frame(Gene = row.names(C9), C9)
C10 <- data.frame(Gene = row.names(C10), C10)
C11 <- data.frame(Gene = row.names(C11), C11)
C12 <- data.frame(Gene = row.names(C12), C12)
C13 <- data.frame(Gene = row.names(C13), C13)
C14 <- data.frame(Gene = row.names(C14), C14)
C15 <- data.frame(Gene = row.names(C15), C15)
C16 <- data.frame(Gene = row.names(C16), C16)
C17 <- data.frame(Gene = row.names(C17), C17)
C18 <- data.frame(Gene = row.names(C18), C18)
C19 <- data.frame(Gene = row.names(C19), C19)
C20 <- data.frame(Gene = row.names(C20), C20)
C21 <- data.frame(Gene = row.names(C21), C21)
C22 <- data.frame(Gene = row.names(C22), C22)

#Save
require(openxlsx)
list_of_datasets <- list("C0" = C0, "C1" = C1, "C2" = C2, "C3" = C3, "C4" = C4, "C5" = C5, "C6" = C6, "C7" = C7, "C8" = C8, "C9" = C9, "C10" = C10, "C11" = C11, "C12" = C12, "C13" = C13, "C14" = C14, "C15" = C15, "C16" = C16, "C17" = C17, "C18" = C18, "C19" = C19, "C20" = C20, "C21" = C21, "C22" = C22)
write.xlsx(list_of_datasets, file = "6Sample_ClusterMarkers_Updated2.xlsx")



#load saved file
Data.Combined <- readRDS("6Samples_subclustered.rds")

# Cell-Cycle scoring
# read in a list of cell cycle markers, from Tirosh et al, 2015
cellcyclegenes <- "regev_lab_cell_cycle_genes.lowercase.txt"
cc.genes <- readLines(con=cellcyclegenes)
length(cc.genes)
# segregate this list into markers of G2/M phase and markers of S
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:98]
DefaultAssay(object = Data.Combined) <- "RNA"
# assign cell-cycle scores
Data.Combined <- CellCycleScoring(Data.Combined, s.features = s.genes, g2m.features = g2m.genes, set.ident= FALSE)
head(Data.Combined@meta.data)
#Visualize results over umap
DimPlot(Data.Combined, reduction = "umap", group.by = "Phase", split.by = "KO")
#Get number of cells in each phase separated by cluster and group
table(Data.Combined@meta.data[["Phase"]], Data.Combined@meta.data$sub_cluster, Data.Combined@meta.data$KO)



Sub2 <- subset(immune.combined, idents = c("2.CTRL", "2.KO"))
saveRDS(Sub2, file = "PECSubcluster2_Only.rds")

Sub3 <- subset(immune.combined, idents = c("3.CTRL", "3.KO"))
saveRDS(Sub3, file = "PECSubcluster3_Only.rds")


###########################################################
###Subclustering of PEC cluster (Cluster 3 in 6 Sample)####
Data.Combined.c3 <- subset(Data.Combined, idents = "3")
DefaultAssay(object = Data.Combined.c3) <- "integrated"
Data.Combined.c3 <- FindNeighbors(Data.Combined.c3, dims = 1:10)
Data.Combined.c3 <- FindClusters(Data.Combined.c3, resolution = 0.5)

# Generate a new column called sub_cluster in the metadata
Data.Combined$sub_cluster <- as.character(Idents(Data.Combined))

# Change the information of cells containing sub-cluster information
Data.Combined$sub_cluster[Cells(Data.Combined.c3)] <- paste("c3",Idents(Data.Combined.c3))
DimPlot(Data.Combined, group.by = "sub_cluster")

# Visualization of UMAP
DimPlot(Data.Combined, reduction = "umap", group.by = "sub_cluster")
DimPlot(Data.Combined, reduction = "umap", split.by = "KO")
DimPlot(Data.Combined, reduction = "umap", label = FALSE)
DimPlot(Data.Combined.c3, reduction = "umap", label = FALSE)
DimPlot(immune.combined, reduction = "umap", label = FALSE)

DefaultAssay(Data.Combined) <- "RNA"
Idents(Data.Combined) <- "seurat_clusters"
FeaturePlot(Data.Combined, features = "Cd44", split.by = "KO")
VlnPlot(Data.Combined, features = c("Itgav"), ncol = 1, split.by = "KO", pt.size = 0)+facet_wrap(~ "Receptors")

VlnPlot(Data.Combined.c3, features = c("Cd44", "Slc5a12", "Havcr1", "Nphs1", "Sytl2"), ncol = 1, split.by = "KO", pt.size = 0)
RidgePlot(Data.Combined.c3, features = "Cd44", sort = "decreasing")
FeatureScatter(Data.Combined.c3, feature1= "Cd44", feature2 = "Havcr1")

VlnPlot(Data.Combined.c3, features = c("Cd44", "Slc5a12", "Slc7a7", "Havcr1"), ncol = 1, pt.size = 0)

#For Differential Expression between Control and KO of the same cluster
immune.combined <- Data.Combined.c3
immune.combined$celltype.KO <- paste(Idents(immune.combined), immune.combined$KO, sep = ".")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.KO"

C0 <- FindMarkers(immune.combined, ident.1 = "0_KO", ident.2 = "0_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C1 <- FindMarkers(immune.combined, ident.1 = "1_KO", ident.2 = "1_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C2 <- FindMarkers(immune.combined, ident.1 = "2_KO", ident.2 = "2_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C3 <- FindMarkers(immune.combined, ident.1 = "3_KO", ident.2 = "3_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C4 <- FindMarkers(immune.combined, ident.1 = "4_KO", ident.2 = "4_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C5 <- FindMarkers(immune.combined, ident.1 = "5_KO", ident.2 = "5_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C6 <- FindMarkers(immune.combined, ident.1 = "6_KO", ident.2 = "6_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C7 <- FindMarkers(immune.combined, ident.1 = "7_KO", ident.2 = "7_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C8 <- FindMarkers(immune.combined, ident.1 = "8_KO", ident.2 = "8_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)

#add row names (Genes) as its own column 
C0 <- data.frame(Gene = row.names(C0), C0)
C1 <- data.frame(Gene = row.names(C1), C1)
C2 <- data.frame(Gene = row.names(C2), C2)
C3 <- data.frame(Gene = row.names(C3), C3)
C4 <- data.frame(Gene = row.names(C4), C4)
C5 <- data.frame(Gene = row.names(C5), C5)
C6 <- data.frame(Gene = row.names(C6), C6)
C7 <- data.frame(Gene = row.names(C7), C7)
C8 <- data.frame(Gene = row.names(C8), C8)

#For Differential Expression between a single cluster and all other cells (Cluster 1 vs all other clusters)
immune.combined <- Data.Combined.c3
immune.combined$celltype.KO <- paste(Idents(immune.combined), immune.combined$KO, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype"

C0 <- FindMarkers(immune.combined, ident.1 = "0", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C1 <- FindMarkers(immune.combined, ident.1 = "1", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C2 <- FindMarkers(immune.combined, ident.1 = "2", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C3 <- FindMarkers(immune.combined, ident.1 = "3", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C4 <- FindMarkers(immune.combined, ident.1 = "4", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C5 <- FindMarkers(immune.combined, ident.1 = "5", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C6 <- FindMarkers(immune.combined, ident.1 = "6", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C7 <- FindMarkers(immune.combined, ident.1 = "7", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C8 <- FindMarkers(immune.combined, ident.1 = "8", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)

#add row names (Genes) as its own column 
C0 <- data.frame(Gene = row.names(C0), C0)
C1 <- data.frame(Gene = row.names(C1), C1)
C2 <- data.frame(Gene = row.names(C2), C2)
C3 <- data.frame(Gene = row.names(C3), C3)
C4 <- data.frame(Gene = row.names(C4), C4)
C5 <- data.frame(Gene = row.names(C5), C5)
C6 <- data.frame(Gene = row.names(C6), C6)
C7 <- data.frame(Gene = row.names(C6), C6)
C8 <- data.frame(Gene = row.names(C6), C6)


#Save
require(openxlsx)
list_of_datasets <- list("C0" = C0, "C1" = C1, "C2" = C2, "C3" = C3, "C4" = C4, "C5" = C5, "C6" = C6, "C7" = C7, "C8" = C8)
write.xlsx(list_of_datasets, file = "PEC_subcluster_markers.xlsx")

#To output average expression split by group
cluster.averages <- AverageExpression(immune.combined, assays = "RNA", add.ident = "KO")
write.csv(cluster.averages[["RNA"]], file = "Podocyte_subclusters_ClusterAverages.csv")

# Get number of cells per cluster and per sample/group of origin
table(immune.combined@meta.data[["seurat_clusters"]], immune.combined@meta.data$orig.ident)
table(immune.combined@meta.data[["seurat_clusters"]], Data.Combined@meta.data$KO)



#Save
require(openxlsx)
list_of_datasets <- list("C0" = C0, "C1" = C1, "C2" = C2, "C3" = C3, "C4" = C4, "C5" = C5, "C6" = C6, "C7" = C7)
write.xlsx(list_of_datasets, file = "PEC_subclusters_CTRLvsKO.xlsx")



###########################################################
###Subclustering of Podocytes cluster (Cluster 15 in 6 Sample)####
Data.Combined.c15 <- subset(Data.Combined, idents = "15")
DefaultAssay(object = Data.Combined.c15) <- "integrated"
Data.Combined.c15 <- FindNeighbors(Data.Combined.c15, dims = 1:10)
Data.Combined.c15 <- FindClusters(Data.Combined.c15, resolution = 0.5)

DimPlot(Data.Combined.c15, reduction = "umap", group.by = "sub_cluster")

# Generate a new column called sub_cluster in the metadata
Data.Combined$sub_cluster <- as.character(Idents(Data.Combined))

# Change the information of cells containing sub-cluster information
Data.Combined$sub_cluster[Cells(Data.Combined.c15)] <- paste("c15",Idents(Data.Combined.c15))
DimPlot(Data.Combined.c15, group.by = "sub_cluster")

DefaultAssay(Data.Combined.c15) <- "RNA"
FeaturePlot(Data.Combined.c15, features = "Nphs1")
VlnPlot(Data.Combined.c15, features = c("Nphs1", "Nphs2", "Fn1", "Stat3"), ncol = 1, split.by = "KO", pt.size = 0)
RidgePlot(Data.Combined.c15, features = "Klf4", sort = "decreasing")
FeatureScatter(Data.Combined.c15, feature1= "Nphs1", feature2 = "Nphs2")

#For Differential Expression between Control and KO of the same cluster
immune.combined <- Data.Combined.c15
immune.combined$celltype.KO <- paste(Idents(immune.combined), immune.combined$KO, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype.KO"

C0 <- FindMarkers(immune.combined, ident.1 = "0_KO", ident.2 = "0_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C1 <- FindMarkers(immune.combined, ident.1 = "1_KO", ident.2 = "1_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C2 <- FindMarkers(immune.combined, ident.1 = "2_KO", ident.2 = "2_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C3 <- FindMarkers(immune.combined, ident.1 = "3_KO", ident.2 = "3_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C4 <- FindMarkers(immune.combined, ident.1 = "4_KO", ident.2 = "4_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)
C5 <- FindMarkers(immune.combined, ident.1 = "5_KO", ident.2 = "5_CTRL", min.pct = 0.1, logfc.threshold = 0.20, verbose = FALSE)

#add row names (Genes) as its own column 
C0 <- data.frame(Gene = row.names(C0), C0)
C1 <- data.frame(Gene = row.names(C1), C1)
C2 <- data.frame(Gene = row.names(C2), C2)
C3 <- data.frame(Gene = row.names(C3), C3)
C4 <- data.frame(Gene = row.names(C4), C4)
C5 <- data.frame(Gene = row.names(C5), C5)

#Save
require(openxlsx)
list_of_datasets <- list("C0" = C0, "C1" = C1, "C2" = C2, "C3" = C3, "C4" = C4, "C5" = C5)
write.xlsx(list_of_datasets, file = "Podocyte_subclusters_CTRLvsKO.xlsx")

#For Differential Expression between a single cluster and all other cells (Cluster 1 vs all other clusters)
immune.combined <- Data.Combined.c15
immune.combined$celltype.KO <- paste(Idents(immune.combined), immune.combined$KO, sep = "_")
immune.combined$celltype <- Idents(immune.combined)
Idents(immune.combined) <- "celltype"

C0 <- FindMarkers(immune.combined, ident.1 = "0", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C1 <- FindMarkers(immune.combined, ident.1 = "1", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C2 <- FindMarkers(immune.combined, ident.1 = "2", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C3 <- FindMarkers(immune.combined, ident.1 = "3", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C4 <- FindMarkers(immune.combined, ident.1 = "4", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C5 <- FindMarkers(immune.combined, ident.1 = "5", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)
C6 <- FindMarkers(immune.combined, ident.1 = "6", ident.2 = NULL, min.pct = 0.2, logfc.threshold = 0.20, verbose = FALSE)

#add row names (Genes) as its own column 
C0 <- data.frame(Gene = row.names(C0), C0)
C1 <- data.frame(Gene = row.names(C1), C1)
C2 <- data.frame(Gene = row.names(C2), C2)
C3 <- data.frame(Gene = row.names(C3), C3)
C4 <- data.frame(Gene = row.names(C4), C4)
C5 <- data.frame(Gene = row.names(C5), C5)
C6 <- data.frame(Gene = row.names(C6), C6)


#Save
require(openxlsx)
list_of_datasets <- list("C0" = C0, "C1" = C1, "C2" = C2, "C3" = C3, "C4" = C4, "C5" = C5)
write.xlsx(list_of_datasets, file = "Podocyte_subcluster_markers.xlsx")

#To output average expression split by group
cluster.averages <- AverageExpression(immune.combined, assays = "RNA", add.ident = "KO")
write.csv(cluster.averages[["RNA"]], file = "Podocyte_subclusters_ClusterAverages.csv")

# Get number of cells per cluster and per sample/group of origin
table(immune.combined@meta.data[["seurat_clusters"]], immune.combined@meta.data$orig.ident)
table(immune.combined@meta.data[["seurat_clusters"]], Data.Combined@meta.data$KO)


