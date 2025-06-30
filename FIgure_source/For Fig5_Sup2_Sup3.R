#for annotation
VlnPlot(Whole1,features = c("Mki67","Pcna","Nes","Ccnd1","Hey1","Sfrp1","Reln","Cntn2","Neurod1","Sox2","Sox9","Fgfr1","Gdf10","Pecam1","Vwf","Calb1","Lhx1","Foxp2","Grid2","Itpr1","Slc17a6","Tbr1","Pax2","Gad1","Pdgfra","Olig1","Olig2","Aif1","Cx3cr1","Eomes","Rsph1","Foxj1"),stack = T, split.by = "seurat_clusters",flip = T)+
  coord_cartesian(clip = "off")
#set colors
library(RColorBrewer)
my_colors <- c(brewer.pal(12, "Paired"), brewer.pal(12, "Set3"))
#named each cluster
new.cluster.ids1 <- c("P0P6 GCP", "P6 IGL-GC", "P0P6 PWM cell", "P6 IGL-GC", "P6 iEGL-GC", "P6 iEGL-GC","P6 BG-like cell", "P0P6 Endothelial cell", "P0 Differentiating exNeuron","P6 Purkinje cell","P0 Purkinje cell","P0 Astroglial lineage cell","P0P6 Deep Cerebellar Nuclei cell","P6 Astrocyte-like cell","P0 GCP","P0P6 Inhibitory Neuron","P6 Undetermined","P0P6 Oligodendrocyte lineage cell","P6 Inhibitory Neuron","P6 Microglia","P0P6 Deep Cerebellar Nuclei cell","P6 Unipolar Blush cell", "P0P6 Ciliated cell","P6 Choroid plexus cell")
names(new.cluster.ids1) <- levels(Whole1)
Whole2 <- RenameIdents(Whole1, new.cluster.ids1)
VlnPlot(Whole2,features = c("Mki67","Pcna","Nes","Ccnd1","Hey1","Sfrp1","Reln","Cntn2","Neurod1","Sox2","Sox9","Fgfr1","Gdf10","Pecam1","Vwf","Calb1","Lhx1","Foxp2","Grid2","Itpr1","Slc17a6","Tbr1","Pax2","Gad1","Pdgfra","Olig1","Olig2","Aif1","Cx3cr1","Eomes","Rsph1","Foxj1"),stack = T, split.by = "seurat_clusters",flip = T)+
  scale_fill_manual(values=my_colors)+theme(strip.text.y = element_text(size=10))+ coord_cartesian(clip = "off")

new.cluster.ids2 <- c("P0P6 GCP", "P6 IGL-GC-1", "P0P6 PWM cell", "P6 IGL-GC-2", "P6 iEGL-GC-1", "P6 iEGL-GC-2","P6 BG-like cell", "P0P6 Endothelial cell", "P0 Differentiating exNeuron","P6 Purkinje cell","P0 Purkinje cell","P0 Astroglial lineage cell","P0P6 Deep Cerebellar Nuclei cell-1","P6 Astrocyte-like cell","P0 GCP","P0P6 Inhibitory Neuron","P6 Undetermined","P0P6 Oligodendrocyte lineage cell","P6 Inhibitory Neuron","P6 Microglia","P0P6 Deep Cerebellar Nuclei cell-2","P6 Unipolar Blush cell", "P0P6 Ciliated cell","P6 Choroid plexus cell")
names(new.cluster.ids2) <- levels(Whole1)
Whole3 <- RenameIdents(Whole1, new.cluster.ids2)
saveRDS(Whole3,"~/Desktop/annotatedcluster.rds")


DimPlot(Whole3, repel = FALSE, label.size = 8,pt.size = 1) +
  scale_color_manual(values = my_colors) +theme(strip.text =element_text(size=18)) +NoAxes()
q1 <- DimPlot(Whole3, repel = FALSE, label.size = 8,pt.size = 1) +
  scale_color_manual(values = my_colors) +theme(strip.text =element_text(size=18)) +NoAxes()
ggsave("~/Desktop/UMAPc24-1.pdf", plot=q1, width = 15, height = 10)


DimPlot(Whole3, repel = FALSE,split.by = "orig.ident", label.size = 8) +
  scale_color_manual(values = my_colors) +theme(strip.text =element_text(size=18)) +NoAxes()
q2 <- DimPlot(Whole3, repel = FALSE,split.by = "orig.ident", label.size = 8) +
  scale_color_manual(values = my_colors) +theme(strip.text =element_text(size=18)) +NoAxes()
ggsave("~/Desktop/UMAPc24-2.pdf", plot=q2, width = 15, height = 10)

#extract data for spacial information
c0 <- subset(Whole1,idents = 0)
write.csv(Idents(c0),"/Users/kyokasuyama/Desktop/Xenium/csv/c0forXenium.csv")

c1 <- subset(Whole1,idents = 1)
write.csv(Idents(c1),"/Users/kyokasuyama/Desktop/Xenium/csv/c1forXenium.csv")

c2 <- subset(Whole1,idents = 2)
write.csv(Idents(c2),"/Users/kyokasuyama/Desktop/Xenium/csv/c2forXenium.csv")

c3 <- subset(Whole1,idents = 3)
write.csv(Idents(c3),"/Users/kyokasuyama/Desktop/Xenium/csv/c3forXenium.csv")

c4 <- subset(Whole1,idents = 4)
write.csv(Idents(c4),"/Users/kyokasuyama/Desktop/Xenium/csv/c4forXenium.csv")

c5 <- subset(Whole1,idents = 5)
write.csv(Idents(c5),"/Users/kyokasuyama/Desktop/Xenium/csv/c5forXenium.csv")

c6 <- subset(Whole1,idents = 6)
write.csv(Idents(c6),"/Users/kyokasuyama/Desktop/Xenium/csv/c6forXenium.csv")

c7 <- subset(Whole1,idents = 7)
write.csv(Idents(c7),"/Users/kyokasuyama/Desktop/Xenium/csv/c7forXenium.csv")

c8 <- subset(Whole1,idents = 8)
write.csv(Idents(c8),"/Users/kyokasuyama/Desktop/Xenium/csv/c8forXenium.csv")

c9 <- subset(Whole1,idents = 9)
write.csv(Idents(c9),"/Users/kyokasuyama/Desktop/Xenium/csv/c9forXenium.csv")

c10 <- subset(Whole1,idents = 10)
write.csv(Idents(c10),"/Users/kyokasuyama/Desktop/Xenium/csv/c10forXenium.csv")

c11 <- subset(Whole1,idents = 11)
write.csv(Idents(c11),"/Users/kyokasuyama/Desktop/Xenium/csv/c11forXenium.csv")

c12 <- subset(Whole1,idents = 12)
write.csv(Idents(c12),"/Users/kyokasuyama/Desktop/Xenium/csv/c12forXenium.csv")

c13 <- subset(Whole1,idents = 13)
write.csv(Idents(c13),"/Users/kyokasuyama/Desktop/Xenium/csv/c13forXenium.csv")

c14 <- subset(Whole1,idents = 14)
write.csv(Idents(c14),"/Users/kyokasuyama/Desktop/Xenium/csv/c14forXenium.csv")

c15 <- subset(Whole1,idents = 15)
write.csv(Idents(c15),"/Users/kyokasuyama/Desktop/Xenium/csv/c15forXenium.csv")

c16 <- subset(Whole1,idents = 16)
write.csv(Idents(c16),"/Users/kyokasuyama/Desktop/Xenium/csv/c16forXenium.csv")

c17 <- subset(Whole1,idents = 17)
write.csv(Idents(c17),"/Users/kyokasuyama/Desktop/Xenium/csv/c17forXenium.csv")

c18 <- subset(Whole1,idents = 18)
write.csv(Idents(c18),"/Users/kyokasuyama/Desktop/Xenium/csv/c18forXenium.csv")

c19 <- subset(Whole1,idents = 19)
write.csv(Idents(c19),"/Users/kyokasuyama/Desktop/Xenium/csv/c19forXenium.csv")

c20 <- subset(Whole1,idents = 20)
write.csv(Idents(c20),"/Users/kyokasuyama/Desktop/Xenium/csv/c20forXenium.csv")

c21 <- subset(Whole1,idents = 21)
write.csv(Idents(c21),"/Users/kyokasuyama/Desktop/Xenium/csv/c21forXenium.csv")

c22 <- subset(Whole1,idents = 22)
write.csv(Idents(c22),"/Users/kyokasuyama/Desktop/Xenium/csv/c22forXenium.csv")

c23 <- subset(Whole1,idents = 23)
write.csv(Idents(c23),"/Users/kyokasuyama/Desktop/Xenium/csv/c23forXenium.csv")

#extract Astroglial lineage cell
P0P6Astglineage <- subset(Whole1,idents=c(6,11,13))
#read cell-cycle Scoring and Regression files from web https://satijalab.org/seurat/articles/cell_cycle_vignette.html
exp.mat <- read.table(file = "/Users/kyokasuyama/Desktop/Xenium/Xenium data from Web/cell_cycle_vignette_files/nestorawa_forcellcycle_expressionMatrix.txt")
library(stringr)
#cell cycle scoring
s.genes <- cc.genes$s.genes
s.genes_mouse<-str_to_sentence(s.genes)
g2m.genes <- cc.genes$g2m.genes
g2m.genes_mouse <- str_to_sentence(g2m.genes)
P0P6Astglineage <- JoinLayers(P0P6Astglineage)
P0P6Astglineage <- CellCycleScoring(P0P6Astglineage,s.features = s.genes_mouse,g2m.features = g2m.genes_mouse,set.ident = T)
P0P6Astglineage <- NormalizeData(P0P6Astglineage,assay = "Xenium")
P0P6Astglineage <- ScaleData(P0P6Astglineage, features = rownames(P0P6Astglineage),assay = "Xenium")
P0P6Astglineage <- RunPCA(P0P6Astglineage)
P0P6Astglineage <- RunUMAP(P0P6Astglineage,dims = 1:30)
P0P6Astglineage<-FindNeighbors(P0P6Astglineage,reduction = "pca",dims = 1:30)
P0P6Astglineage <- FindClusters(P0P6Astglineage,resolution = 1.0) 
#cell cyle regression
install.packages("dplyr")
library("dplyr")
P0P6Astglineage <- ScaleData(P0P6Astglineage, vars.to.regress = c("S.Score","G2M.Score"), features = rownames(P0P6Astglineage))
P0P6Astglineage <- FindVariableFeatures(P0P6Astglineage)
P0P6Astglineage <- RunPCA(P0P6Astglineage)
P0P6Astglineage <- RunUMAP(P0P6Astglineage,dims = 1:30)
P0P6Astglineage <- FindNeighbors(P0P6Astglineage,reduction = "pca",dims = 1:30)
P0P6Astglineage <- FindClusters(P0P6Astglineage,resolution = 1.0) 

DimPlot(P0P6Astglineage)
VlnPlot(P0P6Astglineage, features = c("Gdf10","Gria1","Sept4","Sox2","Sox9","Hopx","Pcna"),stack = T,flip = T,split.by = "ident")+scale_fill_manual(values=my_colors)+theme(
  strip.text.y = element_text(size = 18))


#extract spacial information
sc0 <- subset(P0P6Astglineage,idents = 0)
write.csv(Idents(sc0),"/Users/kyokasuyama/Desktop/Xenium/csv/sc0forXenium.csv")
sc1 <- subset(P0P6Astglineage,idents = 1)
write.csv(Idents(sc1),"/Users/kyokasuyama/Desktop/Xenium/csv/sc1forXenium.csv")
sc2 <- subset(P0P6Astglineage,idents = 2)
write.csv(Idents(sc2),"/Users/kyokasuyama/Desktop/Xenium/csv/sc2forXenium.csv")
sc3 <- subset(P0P6Astglineage,idents = 3)
write.csv(Idents(sc3),"/Users/kyokasuyama/Desktop/Xenium/csv/sc3forXenium.csv")
sc4 <- subset(P0P6Astglineage,idents = 4)
write.csv(Idents(sc4),"/Users/kyokasuyama/Desktop/Xenium/csv/sc4forXenium.csv")
sc5 <- subset(P0P6Astglineage,idents = 5)
write.csv(Idents(sc5),"/Users/kyokasuyama/Desktop/Xenium/csv/sc5forXenium.csv")
sc6 <- subset(P0P6Astglineage,idents = 6)
write.csv(Idents(sc6),"/Users/kyokasuyama/Desktop/Xenium/csv/sc6forXenium.csv")
sc7 <- subset(P0P6Astglineage,idents = 7)
write.csv(Idents(sc7),"/Users/kyokasuyama/Desktop/Xenium/csv/sc7forXenium.csv")
sc8 <- subset(P0P6Astglineage,idents = 8)
write.csv(Idents(sc8),"/Users/kyokasuyama/Desktop/Xenium/csv/sc8forXenium.csv")
sc9 <- subset(P0P6Astglineage,idents = 9)
write.csv(Idents(sc9),"/Users/kyokasuyama/Desktop/Xenium/csv/sc9forXenium.csv")

new.cluster.ids3 <- c("P0 BG + IGL astrocyte", "P6 BGLP-1", "P6 WM astrocyte", "P6 BG","P6 BGLP-2", "P6 IGL astrocyte", "P0 AsLP + WM astrocyte","P0 BGLP", "Undetermined-1", "Undetermined-2")
names(new.cluster.ids3) <- levels(P0P6Astglineage)
P0P6Astglineage1 <- RenameIdents(P0P6Astglineage, new.cluster.ids3)
DimPlot(P0P6Astglineage1, repel = FALSE, label.size = 8,pt.size = 3) +
  scale_color_manual(values = my_colors) +theme(strip.text =element_text(size=18)) +NoAxes()
q3 <- DimPlot(P0P6Astglineage1, repel = FALSE, label.size = 8, pt.size = 3) +
  scale_color_manual(values = my_colors) +theme(strip.text =element_text(size=18)) +NoAxes()
ggsave("~/Desktop/UMAP0P6Ast-1.pdf", plot=q3, width = 15, height = 10)

DimPlot(P0P6Astglineage1, repel = FALSE, label.size = 8,split.by = "orig.ident",pt.size = 2) +
  scale_color_manual(values = my_colors) +theme(strip.text =element_text(size=18)) +NoAxes()
q4 <-DimPlot(P0P6Astglineage1, repel = FALSE, label.size = 8,split.by = "orig.ident",pt.size = 2) +
  scale_color_manual(values = my_colors) +theme(strip.text =element_text(size=18)) +NoAxes()
ggsave("~/Desktop/UMAP0P6Ast-2.pdf", plot=q4, width = 15, height = 10)

VlnPlot(P0P6Astglineage1, features = c("Gdf10","Gria1","Sept4","Sox2","Sox9","Hopx","Pcna"),stack = T,flip = T,split.by = "ident")+scale_fill_manual(values=my_colors)+theme(
  strip.text.y = element_text(size = 18))
q5<-VlnPlot(P0P6Astglineage1, features = c("Gdf10","Gria1","Sept4","Sox2","Sox9","Hopx","Pcna"),stack = T,flip = T,split.by = "ident")+scale_fill_manual(values=my_colors)+theme(
  strip.text.y = element_text(size = 18))
ggsave("~/Desktop/UMAP0P6Ast-VLN.pdf", plot=q5, width = 15, height = 10)


#extract P0 Astroglial lineage clusters
P0Astglineage <- subset(P0P6Astglineage,orig.ident=="P0")

library(viridis)
FeaturePlot(P0Astglineage,features = c("Gdf10","Gria1","Pcna","Mki67") ,cols = viridis(100))
q8 <-FeaturePlot (P0Astglineage,features = c("Gdf10","Gria1","Pcna","Mki67") ,cols = viridis(100))
ggsave("~/Desktop/FeatureP0Astglineage.pdf", plot=q8, width = 15, height = 10)

P0Astglineage_c067<-subset(P0Astglineage,ident=c(0,6,7))
my_colors2 <- setNames( my_colors[c(1, 7, 8)],levels(c(0,6,7))  ) 


new.cluster.ids4 <- c("P0 BG + IGL astrocyte",  "P0 AsLP + WM astrocyte","P0 BGLP" )
names(new.cluster.ids4) <- levels(P0Astglineage_c067)
P0Astglineage_c067_1 <- RenameIdents(P0Astglineage_c067, new.cluster.ids4)
VlnPlot(P0Astglineage_c067_1, features = c("Gdf10","Gria1","Sox2","Sox9","Hopx","Pcna"),stack = T,flip = T,split.by = "ident")+scale_fill_manual(values=my_colors2)+theme(strip.text.y = element_text(size=15))
q6<-VlnPlot(P0Astglineage_c067_1, features = c("Gdf10","Gria1","Sox2","Sox9","Hopx","Pcna"),stack = T,flip = T,split.by = "ident")+scale_fill_manual(values=my_colors2)+theme(strip.text.y = element_text(size=15))
ggsave("~/Desktop/VLNP0Astglineage_c067.pdf", plot=q6, width = 15, height = 10)

#extract P6 Astroglial lineage clusters
P6Astglineage <- subset(P0P6Astglineage,orig.ident=="P6")

FeaturePlot(P6Astglineage,features = c("Gdf10","Gria1","Pcna","Mki67") ,cols = viridis(100))
q9 <-FeaturePlot (P6Astglineage,features = c("Gdf10","Gria1","Pcna","Mki67") ,cols = viridis(100))
ggsave("~/Desktop/FeatureP6Astglineage.pdf", plot=q9, width = 15, height = 10)

P6Astglineage_c12345<-subset(P6Astglineage,ident=c(1,2,3,4,5))
my_colors3 <- setNames( my_colors[c(2,3,4,5,6)],levels(c(1,2,3,4,5))  ) 
new.cluster.ids5 <- c("P6 BGLP-1",  "P6 WM astrocyte","P6 BG","P6 BGLP-2","P6 IGL astrocyte" )
names(new.cluster.ids5) <- levels(P6Astglineage_c12345)
P6Astglineage_c12345_1 <- RenameIdents(P6Astglineage_c12345, new.cluster.ids5)
VlnPlot(P6Astglineage_c12345_1, features = c("Gdf10","Gria1","Sox2","Sox9","Hopx","Pcna"),stack = T,flip = T,split.by = "ident")+scale_fill_manual(values=my_colors3)+theme(strip.text.y = element_text(size=15))
q7<-VlnPlot(P6Astglineage_c12345_1, features = c("Gdf10","Gria1","Sox2","Sox9","Hopx","Pcna"),stack = T,flip = T,split.by = "ident")+scale_fill_manual(values=my_colors3)+theme(strip.text.y = element_text(size=15))
ggsave("~/Desktop/VLNP6Astglineage_c12345.pdf", plot=q7, width = 15, height = 10)
