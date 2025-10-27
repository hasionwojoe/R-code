library(Seurat)
library(cowplot)
library(ggplot2)
library(tidyverse)
library(scCustomize)

basic_path = "YOUR_PATH_HERE"

metadata_file = paste(basic_path,"metadata.csv",sep="")
metadata_df = read.csv2(metadata_file,header=TRUE,sep="\t",stringsAsFactors = FALSE)

seurat_object_list = list()


################################################################

#load seurat objects and add them to list

for(i in 1:length(rownames(metadata_df))){
  current_sample = as.character(metadata_df$Sample[i])
  current_basic_path = as.character(metadata_df$Basic_Path[i])
  current_seurat_object_file = as.character(metadata_df$Seurat_Object[i])
  current_seurat_object <- readRDS(current_seurat_object_file)
  current_seurat_object$sample = current_sample
  seurat_object_list = append(seurat_object_list,current_seurat_object)
}


################################################################
#perform integration in the classical way

min_cell_number <- min(sapply(seurat_object_list, ncol))

# Enable parallelization with future library
plan("multicore", workers = 10)
#set Max size for objects
options(future.globals.maxSize = 30000 * 1024^2)
options(future.seed=TRUE)


#First Analysis
combined_anchors <- FindIntegrationAnchors(object.list = seurat_object_list, dims = 1:20)

# create list of common genes to keep
to_integrate <- Reduce(intersect, lapply(combined_anchors@object.list, rownames))

# integrate data and keep full geneset
combined_data <- IntegrateData(anchorset = combined_anchors, dims = 1:20, features.to.integrate = to_integrate)

#save rds file
saveRDS(combined_data, file = paste(basic_output_dir,"combined_data.rds",sep=""))

#define default assay
DefaultAssay(combined_data) <- "integrated"

# Run the standard workflow for visualization and clustering
combined_data <- ScaleData(combined_data, verbose = FALSE)
combined_data <- RunPCA(combined_data, npcs = 20, verbose = FALSE)

# Clustering
combined_data <- RunUMAP(combined_data, reduction = "pca", dims = 1:20)
combined_data <- FindNeighbors(combined_data, reduction = "pca", dims = 1:20)
combined_data <- FindClusters(combined_data, resolution = 0.3)

######################################################
#Visualization

#by cluster
p2 <- DimPlot(combined_data, reduction = "umap", label = TRUE,raster=FALSE)
ggsave(paste(plot_dir,"umap_cluster_plot.png",sep=""), plot = p2, width=40,height=20,units="cm")

