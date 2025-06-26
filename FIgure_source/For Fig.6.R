#extract P0BGLP
P0BGLP <- subset(P0Astglineage,idents = 7)

#extract P6BGLP and combine  c1 and c4
P6c1c4 <- subset(P6Astglineage,idents = c(1,4))
P6BGLP <-P6c1c4
Idents(P6BGLP) <- "1&4"
P6BGLP

#mergeBGLP
merge_P0P6BGLP <- merge(P0BGLP,y=P6BGLP)

#scaling
merge_P0P6BGLP <-NormalizeData(merge_P0P6BGLP)
merge_P0P6BGLP<- ScaleData(merge_P0P6BGLP)
merge_P0P6BGLP <- JoinLayers(merge_P0P6BGLP)

#extract DEG
BGLP.marker <- FindAllMarkers(merge_P0P6BGLP,only.pos = TRUE)#positive only
write.csv(BGLP.marker,"/Users/kyokasuyama/Desktop/BGLPmarker.csv")
View(BGLP.marker)

#preparation for GO and heatmap

install.packages("clusterProfiler")
library(clusterProfiler)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Mm.eg.db")
library (org.Mm.eg.db)


#heatmap
install.packages("ggplot2")
library(ggplot2)

BGLP.marker %>%dplyr::filter(avg_log2FC >1 & p_val_adj<0.05 &cluster == 7) %>% arrange(avg_log2FC) -> up_7
BGLP.marker %>%dplyr::filter(avg_log2FC >1 & p_val_adj<0.05 &cluster == "1&4")%>% arrange(avg_log2FC) -> up_14
heatmap_gene <- c(up_7$gene,rev(up_14$gene))
df_heat <- data.frame(heatmap_gene)
r<-DoHeatmap(merge_P0P6BGLP, features = df_heat$heatmap_gene,raster = FALSE) + scale_fill_gradientn(colors = c("blue", "white", "red"))　+ theme(axis.text.y = element_text(size=4))
ggsave("~/Desktop/Doheatmap.pdf", plot = r, width = 8, height = 14)

df_heat$heatmap_gene %>% head(n=224)->top224
write.csv(top224,"~/Desktop/top224.csv")
df_heat$heatmap_gene %>% tail(n=87) ->tail87
write.csv(tail87,"~/Desktop/tail87.csv")
#GO analysis
DEG_entrez <- bitr(top224, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Mm.eg.db)
all_genes <- keys(org.Mm.eg.db)
res.go1<- enrichGO(DEG_entrez$ENTREZID,OrgDb =  org.Mm.eg.db, keyType = "ENTREZID", 
                   universe = all_genes,ont = "All",pvalueCutoff = 0.05, 
                   pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = TRUE)
res.go1%>% filter(ONTOLOGY == "BP") -> res.go.bp1
res.go.bp1%>% as.data.frame() %>% write.csv("~/Desktop/top224BP.csv")
barplot(res.go.bp1,showCategory = 20)+
  theme(
    axis.text.y = element_text(size = 33))->s
ggsave("~/Desktop/barplot224.pdf",plot = s, width = 18,height = 24)

DEG_entrez <- bitr(tail87, fromType = "SYMBOL",toType = "ENTREZID",OrgDb = org.Mm.eg.db)
all_genes <- keys(org.Mm.eg.db)
res.go1<- enrichGO(DEG_entrez$ENTREZID,OrgDb =  org.Mm.eg.db, keyType = "ENTREZID", 
                   universe = all_genes,ont = "All",pvalueCutoff = 0.05, 
                   pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10, maxGSSize = 500, readable = TRUE)
res.go1%>% filter(ONTOLOGY == "BP") -> res.go.bp2
res.go.bp2%>% as.data.frame() %>% write.csv("~/Desktop/tail87BP.csv")
barplot(res.go.bp2,showCategory = 20)+
  theme(
    axis.text.y = element_text(size = 33))->t
ggsave("~/Desktop/barplot187.pdf",plot = t, width = 18,height = 24)

#BoxPlot
genes <- c("Gdf10","Gria1","Tnc","Sept4","Nes","Sox4","Sox11","Foxm1","Nfia")
df <- FetchData(merge_P0P6BGLP, vars = genes)
df$cluster <- Idents(merge_P0P6BGLP)
library(tidyr)
df_long <- pivot_longer(df, cols = all_of(genes), names_to = "Gene", values_to = "Expression")
df_long$Gene <- factor(df_long$Gene, levels = genes)
ggplot(df_long, aes(x = cluster, y = Expression, fill = as.factor(cluster))) +
  geom_boxplot() +
  facet_wrap(~Gene, scales = "free_y") + 
  labs(title = "Gene Expression", x = "Cluster", y = "Expression Level") +
  theme_minimal() +
  scale_fill_brewer(palette = "Set3")
