library(tidyverse)
library(magrittr)
library(liana)
library(Seurat)
library(dplyr)
library(future.apply)
plan(multisession, workers = 4) # Adjust 'workers' to the number of cores you want to use
# make it much faster yeeee

############################################### data input
seurat_transformed = readRDS(".../Datasets/Spatial_transcriptomics/Seurat_obj/seurat_transformed")


# NOTE: HAVE TO CHANGE seurat_transformed[[12]]@meta.data[["TLS_anno"]] from "" to "TLS"

#################### Have to chnage "" to "NO_TLS" of sample 12
################################################  run LIANA

# Initialize a list to store raesults
liana_results_list <- list()

# Loop through each Seurat object and run LIANA
for (i in 1:16) {
  seurat_obj <- seurat_transformed[[i]]
  
  Idents(seurat_obj) = seurat_obj$TLS_anno
  # Run LIANA
  liana_result <- liana_wrap(seurat_obj)
  
  # Aggregate results
  liana_result <- liana_result %>%
    liana_aggregate()
  
  # Store the result in the list
  liana_results_list[[paste0("sample_", i)]] <- liana_result
}


# Assign names to the list based on the index of the sample
names(liana_results_list) <- names(seurat_transformed)

# Save the results
saveRDS(liana_results_list, ".../Datasets/Spatial_transcriptomics/Seurat_obj/liana_results_list")

#liana_results_list_ind = readRDS(".../Datasets/Spatial_transcriptomics/Seurat_obj/liana_results_list")
