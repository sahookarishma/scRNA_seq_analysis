setwd("E:/scrnaseq/GSEXXXXX_RAW/tumor")

library(Seurat)
library(SeuratDisk)

mtx_obj <- ReadMtx(mtx = "matrix.mtx.gz",
                   features = "features.tsv.gz",
                   cells = "barcodes.tsv.gz")
seurat_mtx <- CreateSeuratObject(counts = mtx_obj)
seurat_mtx 

dim(seurat_mtx)
write.table(as.matrix(GetAssayData(object = seurat_mtx, slot = "counts")), 
            'data.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)

setwd("E:/scrnaseq/GSEXXXXX_RAW/normal")
mtx_obj2 <- ReadMtx(mtx = "matrix.mtx.gz",
                   features = "features.tsv.gz",
                   cells = "barcodes.tsv.gz")

seurat_mtx2 <- CreateSeuratObject(counts = mtx_obj2)
dim(seurat_mtx2)

write.table(as.matrix(GetAssayData(object = seurat_mtx2, slot = "counts")), 
            'data_normal_1.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)



setwd("E:/scrnaseq/GSEXXXXX_RAW/tumor2")
mtx_obj3 <- ReadMtx(mtx = "matrix.mtx.gz",
                    features = "features.tsv.gz",
                    cells = "barcodes.tsv.gz")

seurat_mtx3 <- CreateSeuratObject(counts = mtx_obj3)
dim(seurat_mtx3)

write.table(as.matrix(GetAssayData(object = seurat_mtx3, slot = "counts")), 
            'data_tumor_2.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)

setwd("E:/scrnaseq/GSEXXXXX_RAW/normal2")
mtx_obj4 <- ReadMtx(mtx = "matrix.mtx.gz",
                    features = "features.tsv.gz",
                    cells = "barcodes.tsv.gz")

seurat_mtx4 <- CreateSeuratObject(counts = mtx_obj4)
dim(seurat_mtx4)

write.table(as.matrix(GetAssayData(object = seurat_mtx4, slot = "counts")), 
            'data_normal_2.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)

setwd("E:/scrnaseq/GSEXXXXX_RAW/metastasis")
mtx_obj5 <- ReadMtx(mtx = "matrix.mtx.gz",
                    features = "features.tsv.gz",
                    cells = "barcodes.tsv.gz")

seurat_mtx5 <- CreateSeuratObject(counts = mtx_obj5)
dim(seurat_mtx5)

write.table(as.matrix(GetAssayData(object = seurat_mtx5, slot = "counts")), 
            'data_metastasis_1.csv', 
            sep = ',', row.names = T, col.names = T, quote = F)


setwd("E:/scrnaseq/GSEXXXXX_RAW")
library(Seurat)
library(ggplot2)
library(tidyverse)
library(gridExtra)

# get data location
dirs <- list.dirs(path = 'E:/scrnaseq/GSEXXXXX_RAW/', recursive = F, full.names = F)

dirs


for(x in dirs){
  name <- gsub('_filtered_feature_cc_matrix','', x)
  
  cts <- ReadMtx(mtx = paste0('E:/scrnaseq/GSEXXXXX_RAW/',x,'/matrix.mtx.gz'),
                 features = paste0('E:/scrnaseq/GSEXXXXX_RAW/',x,'/features.tsv.gz'),
                 cells = paste0('E:/scrnaseq/GSEXXXXX_RAW/',x,'/barcodes.tsv.gz'))
  
  # create seurat objects
  assign(name, CreateSeuratObject(counts = cts))
}

merged_seurat <- merge(metastasis, y = c(normal1, normal2, tumor1, tumor2),
                       add.cell.ids = ls()[1:5],
                       project = 'CC')
merged_seurat

View(merged_seurat@meta.data)
# create a sample column
merged_seurat$sample <- rownames(merged_seurat@meta.data)
View(merged_seurat$sample)

merged_seurat@meta.data <- separate(merged_seurat@meta.data, col = 'sample', into = c('Patient', 'Type', 'Barcode'), 
                                    sep = '_')



View(merged_seurat@meta.data)
merged_seurat$mitoPercent <- PercentageFeatureSet(merged_seurat, pattern='^MT-')

merged_seurat_filtered <- subset(merged_seurat, subset = nCount_RNA > 800 &
                                   nFeature_RNA > 500 &
                                   mitoPercent < 10)
merged_seurat_filtered

merged_seurat


# Visualize QC metrics as a violin plot
VlnPlot(merged_seurat_filtered, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3)
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
merged_seurat_filtered[["percent.mt"]] <- PercentageFeatureSet(merged_seurat_filtered, pattern = "^MT-")

merged_seurat_filtered[["percent.mt"]]
# Visualize QC metrics as a violin plot
VlnPlot(merged_seurat_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(merged_seurat_filtered, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(merged_seurat_filtered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
merged_seurat_filtered <- NormalizeData(merged_seurat_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
merged_seurat_filtered <- NormalizeData(merged_seurat_filtered)


merged_seurat_filtered <- FindVariableFeatures(merged_seurat_filtered, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(merged_seurat_filtered), 10)
plot1 <- VariableFeaturePlot(merged_seurat_filtered)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(merged_seurat_filtered)
merged_seurat_filtered <- ScaleData(merged_seurat_filtered, features = all.genes)
merged_seurat_filtered <- RunPCA(merged_seurat_filtered, features = VariableFeatures(object = merged_seurat_filtered))
print(merged_seurat_filtered[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(merged_seurat_filtered, dims = 1:2, reduction = "pca")
DimPlot(merged_seurat_filtered, reduction = "pca") + NoLegend()
DimHeatmap(merged_seurat_filtered, dims = 1, cells = 500, balanced = TRUE)
ElbowPlot(merged_seurat_filtered)
merged_seurat_filtered <- FindNeighbors(merged_seurat_filtered, dims = 1:10)
merged_seurat_filtered <- FindClusters(merged_seurat_filtered, resolution = 0.5)
head(Idents(merged_seurat_filtered), 5)
merged_seurat_filtered <- RunUMAP(merged_seurat_filtered, dims = 1:10)

DimPlot(merged_seurat_filtered, reduction = "umap")
saveRDS(merged_seurat_filtered, file = "../output/scrnaseq_cc_tutorial.rds")

cluster2.markers <- FindMarkers(merged_seurat_filtered, ident.1 = 1)
head(cluster2.markers, n = 5)

merged_seurat_filtered.markers <- FindAllMarkers(merged_seurat_filtered, only.pos = TRUE)
merged_seurat_filtered.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)






