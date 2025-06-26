installed.packages("Seurat")
library(Seurat)
library(ggplot2)

xenium_path <- "/P0_P6_Adult_mousecerebellum_5K/"
All.obj <- LoadXenium(data.dir = xenium_path,fov = "fov")
P0_Whole_Cell_ID <- read.csv("/P0_whole_cerebellum_cells_ID.csv", row.names = 1)
cell_id <- rownames(P0_Whole_Cell_ID)
P0_Whole.obj <- All.obj[,cell_id]
P6_Whole_Cell_ID <- read.csv("/P6_whole_cerebellum_cells_ID.csv", row.names = 1)
cell_id_P6w <- rownames(P6_Whole_Cell_ID)
P6_Whole.obj <- All.obj[,cell_id_P6w]
P0_Whole.obj$orig.ident <- "P0"
P6_Whole.obj$orig.ident <- "P6"
Whole.merge <- merge(P0_Whole.obj , y= P6_Whole.obj)
Whole1 <- subset(Whole.merge, nCount_Xenium > 300 & nFeature_Xenium >300 )
all.genes <- row.names(Whole1)
Whole1<- NormalizeData(Whole1, assay = "Xenium")
Whole1<- ScaleData(Whole1, features = all.genes ,assay = "Xenium")
Whole1 <- FindVariableFeatures(Whole1)
Whole1 <- RunPCA(Whole1, npcs = 30, features = all.genes)
Whole1 <- RunUMAP(Whole1, dims = 1:30)
Whole1 <- FindNeighbors(Whole1, reduction = "pca", dims = 1:30)
Whole1 <- FindClusters(Whole1, resolution = 1.0,random.seed = 777)