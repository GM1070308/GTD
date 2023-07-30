library(Seurat)
library(tidyverse)
library(data.table)
library(tibble)
library(dplyr)
library(ggplot2)
library(ggsci)
setwd('~/EPIC/')

#Read single cell count file and metadata
raw_data <- fread('raw_data_10x.txt')
raw_data <- column_to_rownames(raw_data, 'Gene')

metadata <- read.csv('meta_10x.txt', sep = '\t')
meta_tcells <- filter(metadata, annotation == 'Tcells')

############# Blood T cells #############
raw_bTcells <- raw_data[, colnames(raw_data) %in% rownames(filter(meta_tcells, final_cluster !=8 ))]

btcells <- CreateSeuratObject(counts = raw_bTcells)

str(btcells@meta.data)
btcells <- NormalizeData(btcells)
btcells <- FindVariableFeatures(btcells)
btcells <- ScaleData(btcells)
btcells <- RunPCA(btcells, features = VariableFeatures(object = btcells))
ElbowPlot(btcells)
btcells <- FindNeighbors(btcells, dims = 1:15)
btcells <- FindClusters(btcells, resolution = c(0.03))
DimPlot(btcells, group.by = 'RNA_snn_res.0.03', label = TRUE)
Idents(btcells) <- 'RNA_snn_res.0.03'

all_markers <-FindAllMarkers(btcells, 
                             min.pct =  0.25, 
                             min.diff.pct = 0.25)

#Retain marker genes that are overexpressed
all_markers <- filter(all_markers, avg_log2FC > 0)

gen_marker_table <- function(x){
  all_markers[all_markers$cluster == x, ] %>%
    head(n=10)
}

top10_blood <- map_dfr(0:6, gen_marker_table)

btcells <- RunUMAP(btcells,dims= 1:15)

umap_blood <- as.data.frame(btcells@reductions$umap@cell.embeddings)
b_clusters <- btcells@meta.data[5]
umap_blood <- merge(x=umap_blood, y=b_clusters, by = 0)
umap_blood$seurat_clusters <- sub('0', 'bTreg', umap_blood$seurat_clusters)
umap_blood$seurat_clusters <- sub('1', 'bMAIT', umap_blood$seurat_clusters)
umap_blood$seurat_clusters <- sub('2', 'bCD4', umap_blood$seurat_clusters)
umap_blood$seurat_clusters <- as.character(umap_blood$seurat_clusters)
umap_blood <- column_to_rownames(umap_blood, 'Row.names')

#UMAP plot

colours <- pal_npg("nrc")(10)

ggplot(umap_blood, aes(x = UMAP_1, y = UMAP_2, colour = seurat_clusters)) +
  geom_point() +
  theme_void() +
  labs(x = 'UMAP 1', y = 'UMAP 2') +
  scale_colour_manual(name = NULL,values = colours, labels = c('CD4+', 'MAIT','Treg')) +
  annotate("segment", 
           x = -15, xend = -15 + c(5, 0), 
           y = -6, yend = -6 + c(0, 2.5), 
           arrow = arrow(type = "closed", length = unit(10, 'pt'))) +
  annotate( "text", x=-12.5, y = -6.2, label = 'UMAP 1', colour = 'black') +
  annotate( "text", x=-15.5, y = -4.8, label = 'UMAP 2', colour = 'black', angle = 90) +
  guides(shape = guide_legend(order = 3)) +
  ggtitle('Blood') +
  theme(legend.position = c(.5,.02),
        legend.direction="horizontal",
        legend.text=element_text(size=15),
        legend.key.size = unit(1, 'line'), 
        plot.title = element_text(hjust = 0.5, size = 17)) +
  guides(color = guide_legend(override.aes = list(size=4, shape = 15)))
  
ggsave('blood_umap.png')

new_annot_btcells <- umap_blood[3]

#### Decidual T cells ########

raw_dTcells <- raw_data[, colnames(raw_data) %in% rownames(filter(meta_tcells, final_cluster == 8 ))]

dtcells <- CreateSeuratObject(counts = raw_dTcells)

dtcells <- NormalizeData(dtcells)
dtcells <- FindVariableFeatures(dtcells)
dtcells <- ScaleData(dtcells)
dtcells <- RunPCA(dtcells, features = VariableFeatures(object = dtcells))
ElbowPlot(dtcells)
dtcells <- FindNeighbors(dtcells, dims = 1:15)
dtcells <- FindClusters(dtcells, resolution = c(0.1))
DimPlot(dtcells, group.by = 'RNA_snn_res.0.1', label = TRUE)
Idents(tcells) <- 'RNA_snn_res.0.1'

all_markers <-FindAllMarkers(dtcells, 
                             min.pct =  0.25, 
                             min.diff.pct = 0.25)

gen_marker_table <- function(x){
  all_markers[all_markers$cluster == x, ] %>%
    head(n=10)
}

all_markers <- filter(all_markers, avg_log2FC > 0)
top10_decidual <- map_dfr(0:6, gen_marker_table)

dtcells <- RunUMAP(dtcells,dims= 1:15)

umap_deci <- as.data.frame(dtcells@reductions$umap@cell.embeddings)
d_clusters <- dtcells@meta.data[5]
umap_deci <- merge(x=umap_deci, y=d_clusters, by = 0)
umap_deci$seurat_clusters <- sub('0', 'dCD4', umap_deci$seurat_clusters)
umap_deci$seurat_clusters <- sub('1', 'dCD8', umap_deci$seurat_clusters)
umap_deci$seurat_clusters <- sub('2', 'dTreg', umap_deci$seurat_clusters)
umap_deci$seurat_clusters <- as.character(umap_deci$seurat_clusters)
umap_deci <- column_to_rownames(umap_deci, 'Row.names')

ggplot(umap_deci, aes(x = UMAP_1, y = UMAP_2, colour = seurat_clusters)) +
  geom_point() +
  theme_void() +
  labs(x = 'UMAP 1', y = 'UMAP 2') +
  scale_colour_manual(name = NULL,values = c("#E64B35FF","#8491B4FF","#00A087FF"), labels = c('CD4+', 'CD8+', 'Treg')) +
  annotate("segment", 
           x = -7, xend = -7 + c(2.5, 0), 
           y = -6, yend = -6 + c(0, 2.5), 
           arrow = arrow(type = "closed", length = unit(10, 'pt'))) +
  annotate( "text", x=-5.8, y = -6.2, label = 'UMAP 1', colour = 'black') +
  annotate( "text", x=-7.2, y = -4.8, label = 'UMAP 2', colour = 'black', angle = 90) +
  guides(shape = guide_legend(order = 3)) +
  ggtitle('Decidua') +
  theme(legend.position = c(.5,.02),
        legend.direction="horizontal",
        legend.text=element_text(size=15),
        legend.key.size = unit(1, 'line'), 
        plot.title = element_text(hjust = 0.5, size = 17)) +
  guides(color = guide_legend(override.aes = list(size=4, shape = 15)))

ggsave('deci_umap.png')

new_annot_dcells <- umap_deci[3]

new_annot_tcells <- rbind(new_annot_btcells, new_annot_dcells)
colnames(new_annot_tcells) <- 'annotation'

