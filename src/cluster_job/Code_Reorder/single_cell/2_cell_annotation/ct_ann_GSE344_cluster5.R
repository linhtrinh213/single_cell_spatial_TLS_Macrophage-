library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(celldex)
library(SingleR)
library(Seurat)
library (ggplot2)
library(BiocParallel)
set.seed (0)

GSE_integrated = readRDS("/omics/groups/OE0436/internal/Linh/results/GSE_integrated")

# test dataset 
norm_counts <- GSE_integrated[["SCT"]]$data 

# creating ref dataset: GSE_integrated subset for GSE344: 
GSE344_10k = readRDS("/omics/groups/OE0436/data/Linh/Datasets/Annotation/ct_ann_GSE344_10k")
GSE344_ran_value = GSE344_10k@assays[["SCT"]]@data # extract the dgCmatrix 
GSE344_ran_value = as.matrix(GSE344_ran_value)

#labels
celltype = GSE344_10k@meta.data[["x"]] #extracting the cell type 

# run SingleR
ct_ann_GSE344_10kSingleR <- SingleR(test=norm_counts, ref=GSE344_ran_value, labels=celltype, BPPARAM=MulticoreParam(24))

# saveRDS 
saveRDS (ct_ann_GSE344_10kSingleR, "/omics/groups/OE0436/data/Linh/Datasets/Annotation/ct_ann_GSE344_10kSingleR")

