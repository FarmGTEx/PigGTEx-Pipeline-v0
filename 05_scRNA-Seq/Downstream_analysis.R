# Step 4. Downstream analysis using Seurat
library(Seurat)
library(DoubletFinder)
library(ggplot2)
library(dplyr)

brainf.data <- Read10X(data.dir = "Brain/rawdata/CNS0127223/")
brainh.data <- Read10X(data.dir = "Brain/rawdata/CNS0127224/")
braino.data <- Read10X(data.dir = "Brain/rawdata/CNS0127226/")
brainp.data <- Read10X(data.dir = "Brain/rawdata/CNS0127227/")
braint.data <- Read10X(data.dir = "Brain/rawdata/CNS0127228/")

brainf <- CreateSeuratObject(counts = brainf.data, project = "Frontal", min.cells = 3, min.features = 200)
brainf[["percent.mt"]] <- PercentageFeatureSet(brainf, pattern = "^MT-")
brainf <- subset(brainf, subset = nFeature_RNA > 200 & percent.mt < 5)
brainf <- NormalizeData(brainf, normalization.method = "LogNormalize", scale.factor = 10000)
brainf <- FindVariableFeatures(brainf, selection.method = "vst", nfeatures = 2000)

brainh <- CreateSeuratObject(counts = brainh.data, project = "Hypothalamus", min.cells = 3, min.features = 200)
brainh[["percent.mt"]] <- PercentageFeatureSet(brainh, pattern = "^MT-")
brainh <- subset(brainh, subset = nFeature_RNA > 200 & percent.mt < 5)
brainh <- NormalizeData(brainh, normalization.method = "LogNormalize", scale.factor = 10000)
brainh <- FindVariableFeatures(brainh, selection.method = "vst", nfeatures = 2000)

braino <- CreateSeuratObject(counts = braino.data, project = "Occipital", min.cells = 3, min.features = 200)
braino[["percent.mt"]] <- PercentageFeatureSet(braino, pattern = "^MT-")
braino <- subset(braino, subset = nFeature_RNA > 200 & percent.mt < 5)
braino <- NormalizeData(braino, normalization.method = "LogNormalize", scale.factor = 10000)
braino <- FindVariableFeatures(braino, selection.method = "vst", nfeatures = 2000)

brainp <- CreateSeuratObject(counts = brainp.data, project = "Parietal", min.cells = 3, min.features = 200)
brainp[["percent.mt"]] <- PercentageFeatureSet(brainp, pattern = "^MT-")
brainp <- subset(brainp, subset = nFeature_RNA > 200 & percent.mt < 5)
brainp <- NormalizeData(brainp, normalization.method = "LogNormalize", scale.factor = 10000)
brainp <- FindVariableFeatures(brainp, selection.method = "vst", nfeatures = 2000)

braint <- CreateSeuratObject(counts = braint.data, project = "Temporal", min.cells = 3, min.features = 200)
braint[["percent.mt"]] <- PercentageFeatureSet(braint, pattern = "^MT-")
braint <- subset(braint, subset = nFeature_RNA > 200 & percent.mt < 5)
braint <- NormalizeData(braint, normalization.method = "LogNormalize", scale.factor = 10000)
braint <- FindVariableFeatures(braint, selection.method = "vst", nfeatures = 2000)

brain <- merge(brainf, y = c(brainh, braino, brainp, braint), 
               add.cell.ids = c("Frontal","Hypothalamus","Occipital","Parietal","Temporal"), 
               project = "Brain")
brain
head(colnames(brain))
table(brain$orig.ident)

brain.list <- SplitObject(brain, split.by = "orig.ident")
brain.list <- lapply(X = brain.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


# Perform integration
brain.anchors <- FindIntegrationAnchors(object.list = brain.list, dims = 1:20)
brain.combined <- IntegrateData(anchorset = brain.anchors, dims = 1:20)
DefaultAssay(brain.combined) <- "integrated"
brain.combined <- ScaleData(brain.combined, verbose = FALSE)
brain.combined <- RunPCA(brain.combined, npcs = 20, verbose = FALSE)
brain.combined <- JackStraw(brain.combined, num.replicate = 100)
brain.combined <- ScoreJackStraw(brain.combined, dims = 1:20)

JackStrawPlot(brain.combined, dims = 1:15) +
  theme(legend.text = element_text(size=8),
        legend.key.size = unit(0.5, 'lines'),
        legend.key.width = unit(0.01,'cm'),
        legend.title = element_text(size = 9),
        plot.title = element_text(hjust = 0.5, size = 10),
        axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size=9,face="bold"),
        axis.title.y = element_text(size=9,face="bold"))+
  guides(colour = guide_legend(override.aes = list(size=0.9)))

ElbowPlot(brain.combined) 

# t-SNE and Clustering
brain.combined <- RunUMAP(brain.combined, reduction = "pca", dims = 1:20)
brain.combined <- FindNeighbors(brain.combined, reduction = "pca", dims = 1:20)
brain.combined <- FindClusters(brain.combined, resolution = 0.4)


# Identify conserved cell type markers
DefaultAssay(brain.combined) <- "RNA"
cluster0.markers <- FindConservedMarkers(brain.combined, ident.1 = 0, grouping.var = "orig.ident", verbose = FALSE)
cluster1.markers <- FindConservedMarkers(brain.combined, ident.1 = 1, grouping.var = "orig.ident", verbose = FALSE)
cluster2.markers <- FindConservedMarkers(brain.combined, ident.1 = 2, grouping.var = "orig.ident", verbose = FALSE)
cluster3.markers <- FindConservedMarkers(brain.combined, ident.1 = 3, grouping.var = "orig.ident", verbose = FALSE)
cluster4.markers <- FindConservedMarkers(brain.combined, ident.1 = 4, grouping.var = "orig.ident", verbose = FALSE)
cluster5.markers <- FindConservedMarkers(brain.combined, ident.1 = 5, grouping.var = "orig.ident", verbose = FALSE)
cluster6.markers <- FindConservedMarkers(brain.combined, ident.1 = 6, grouping.var = "orig.ident", verbose = FALSE)
cluster7.markers <- FindConservedMarkers(brain.combined, ident.1 = 7, grouping.var = "orig.ident", verbose = FALSE)
cluster8.markers <- FindConservedMarkers(brain.combined, ident.1 = 8, grouping.var = "orig.ident", verbose = FALSE)
saveRDS(brain.combined, "brain_tutorial_new.rds")


# Summarize the cell type
cell_cluster <- as.data.frame(cbind(colnames(brain.combined), brain.combined$seurat_clusters))
colnames(cell_cluster) <- c("cell", "cluster")
write.table(cell_cluster, "cell_cluster.txt", quote = F, sep = "\t", row.names = F)


# Replot UMAP using Azimuth results
data <- readRDS("azimuth_umap.Rds")
type <- read.table("azimuth_pred.tsv", header = T, sep = "\t")
umap <- as.data.frame(data@cell.embeddings)
umap$cell <- rownames(umap)
umap <- merge(type, umap)
umap$predicted.subclass <- as.factor(umap$predicted.subclass)

allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00",
            "#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0",
            "#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE",
            "#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")

class_avg <- umap %>%
  group_by(predicted.subclass) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2))

pdf("brain_azimuth_umap.pdf",width = 6, height=3.5)
ggplot(umap ,aes(x=UMAP_1,y=UMAP_2))+
  geom_point(aes(color=predicted.subclass), size = 3) +
  labs(x="", y="", title = "Brain") +
  scale_color_manual(values = allcolour, 
                     labels=c("Astro(768)","Endo(80)","L2/3 IT(7306)","L5 IT(873)",
                              "Micro-PVM(541)","Oligo(3960)","OPC(3)","VLMC(438)")) +
  theme(text=element_text(size=18)) +
  labs(title = "", x = "UMAP1", y="UMAP2") +
  theme_classic() +
  theme(#panel.background = element_rect(fill='white', colour='black'), 
        panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=18),
        axis.title.x=element_text(colour='black', size=18),
        axis.title.y=element_text(colour='black', size=18),
        axis.text=element_text(colour='black',size=18),
        legend.title=element_blank(),
        legend.text=element_text(size=18),
        legend.key=element_blank())+
  theme(plot.title = element_text(size=22,colour = "black",face = "bold", hjust = 0.5)) +
  guides(colour = guide_legend(override.aes = list(size=5)))
dev.off()
