
set.seed(0) #need to add in the sciprt to clutser!!! 
library(Seurat)
library (glmGamPoi)
library(ggplot2)
library(sctransform)
library(dplyr)
library(ggplot2)
library(patchwork)
library (hdf5r)


GSE_all_sub = readRDS(".../Datasets/Merged/GSE_all_sub") #object no scale. 4 datasets
 
GSE186344 <- subset(GSE_all_sub, subset = dataset_ID == "GSE186344") 
GSE174401 <- subset(GSE_all_sub, subset = dataset_ID == "GSE174401") 
GSE234832 <- subset(GSE_all_sub, subset = dataset_ID == "GSE234832") 

GSE_SCT <- merge(x =  GSE186344, y = c(GSE174401,GSE234832))
GSE_SCT <- SCTransform(GSE_SCT, return.only.var.genes = TRUE)


saveRDS (GSE_SCT, ".../results/GSE_SCT")

