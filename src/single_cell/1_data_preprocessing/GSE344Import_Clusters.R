set.seed(0)

# Load Libraries
library(Seurat)
#library(ggplot2)
#library(sctransform)
#library(scDblFinder)
#library(SingleCellExperiment)
#library("viridis")
#library(knitr)
library(dplyr)
#library(readxl)
#library(SingleR)
#library(celldex)
#library(ggplot2)
library(patchwork)
#library ("harmony")



folder_path = ".../Datasets/GSE186344/" #path to folder
files <- list.files(folder_path) #list all of the files in folder_path
samples = c(888:898, 900, 902, 904) #in the samples vector, specify which samples you wanna take (by specify the unique number of the samples)

detected_files_list <- list() # Create a list to store the detected files for each sample
import_files = list () #list to store all the imported file by the ReaddMtx func of Seurat

# Loop through each element of samples 
for (sample in samples) {
  
  # Detect files containing the current unique sample number
  detected_files <- files[stringr::str_detect(files, as.character(sample))] #in each loop, find the files that contain the unique number of the samples
  
  # Save the detected files to the list
  detected_files_list[[as.character(sample)]] <- detected_files
  
  #paste (string1, string2. sep = ""): combine 2 strings without any space in between 
  matrix = paste(folder_path,detected_files_list[[as.character(sample)]][stringr::str_detect(detected_files_list[[as.character(sample)]],"matrix")],sep = "")
  cells = paste(folder_path,detected_files_list[[as.character(sample)]][stringr::str_detect(detected_files_list[[as.character(sample)]],"barcodes")],sep = "") 
  features = paste(folder_path,detected_files_list[[as.character(sample)]][stringr::str_detect(detected_files_list[[as.character(sample)]],"features")],sep = "")
  
  # Read/Import the file into imported_list. In the detected_files_list, access the current unique samples, search for and take the files that have pattern matrix/barcodes or features  
  import_files [[as.character(sample)]] <- ReadMtx(mtx = matrix, cells = cells, features = features)
}

# Export import_files
# Loop through each matrix in the list and export them individually
#for (i in seq_along(import_files)) {
  # Construct the file path for exporting
#  file_name <- paste0(".../Datasets/GSE186344/matrix_", i, ".mtx")
  
  # Export the matrix using writeMM
#  writeMM(import_files[[i]], file_name)
#}


# create Seurat objects. After creating objects, delect the import_files to free memory

# Create a list to store Seurat objects for each sample
seurat_objects <- list()

# Loop through each sample and create a Seurat object
for (sample in samples) {
    seurat_objects[[as.character(sample)]] <- CreateSeuratObject(
    counts = import_files[[as.character(sample)]],
    project = "GSE344",
    min.cells = 3,
    min.features = 200
  )
} # creating Seurat objects are actually not that long


# Merge the Seurat objects together
merged_seurat <- merge(x = seurat_objects [[1]], y = c(seurat_objects [[2]], seurat_objects [[3]],seurat_objects [[4]], seurat_objects [[5]], seurat_objects [[6]], seurat_objects [[7]],
                                                       seurat_objects [[8]], seurat_objects [[9]], seurat_objects [[10]], seurat_objects [[11]], seurat_objects [[12]],
                                                       seurat_objects [[13]], seurat_objects [[14]]),  add.cell.ids = samples) #add.cells.id so that cell barcode don't overlap accidentally


# Export seurat object in 10X format so that it can be reimported
SaveSeuratRds(merged_seurat,".../Datasets/GSE186344/merged_seurat") #take a bit 

#imported_seurat_obj = readRDS (".../Datasets/GSE186344/merged_seurat")
