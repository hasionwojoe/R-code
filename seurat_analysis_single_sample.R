library(Seurat)
library(cowplot)
library(ggplot2)
library(tidyverse)
library(scCustomize)


#Read raw data
raw_data <- Read10X(data.dir = paste(basic_path,raw_data_path,sep=""))
rownames(raw_data) = genes
colnames(raw_data) = barcodes

seurat_object <- CreateSeuratObject(counts = raw_data, project = "combined", min.cells = 3, min.features = 200)

# calculate mitochondrial QC metrics
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")

seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 20)

########################################################

# Normalize data
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)

########################################################

# Identification of highly variable features
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 15000)

########################################################

#Scaling the data
all.genes <- rownames(seurat_object)
seurat_object <- ScaleData(seurat_object, features = all.genes)

########################################################

# Perform linear dimensional reduction (PCA)
# input can be defined using features argument when using a different subset

seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object))

########################################################

# Determine PC cutoff

ebp <- ElbowPlot(seurat_object)
ggsave(paste(plot_path,"elbow_plot.png",sep=""), plot = ebp)

pc_cutoff = 13

########################################################

# Cluster the cells
seurat_object <- FindNeighbors(seurat_object, dims = 1:pc_cutoff)
seurat_object <- FindClusters(seurat_object, resolution = 0.20)

########################################################

# Run non-linear dimensional reduction (UMAP)

#UMAP
seurat_object <- RunUMAP(seurat_object, dims = 1:pc_cutoff)
umap_plot <- DimPlot(seurat_object, reduction = "umap", label=TRUE, raster=FALSE)
ggsave(paste(plot_path,"umap_cluster_plot.png",sep=""), plot = umap_plot)

########################################################

saveRDS(seurat_object, file = paste(output_path,"seurat_object.rds",sep=""))

