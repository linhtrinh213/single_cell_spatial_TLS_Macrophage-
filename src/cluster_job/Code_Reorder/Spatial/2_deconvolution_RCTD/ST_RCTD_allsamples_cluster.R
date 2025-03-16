
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(readxl)
library(spacexr) #RCTD
library(Matrix)
library(doParallel)



ref = readRDS("/omics/groups/OE0436/internal/Linh/Output/ST/sc_ref") #SCT, with annotation, joined layers (ccRCC)
seurat_transformed = readRDS("/omics/groups/OE0436/internal/Linh/Datasets/Spatial_transcriptomics/Seurat_obj/seurat_transformed")

################################################################## sc ref
counts = ref@assays[["SCT"]]@counts
counts@Dimnames[[1]] = Features(ref)
counts@Dimnames[[2]] = Cells(ref)
cell_types = ref@meta.data[["celltype1"]]
#cell_types <- gsub("/", "_", cell_types)
#any(is.na(cell_types))
#na_cells <- which(is.na(cell_types)) 
#counts <- counts[, -na_cells]
#cell_types <- cell_types[-na_cells]

names(cell_types) = counts@Dimnames[[2]]

#nUMI = GSE_SCT_sub@meta.data[["nCount_RNA"]] 
#names(nUMI) = colnames(GSE_SCT_sub)

### Create the Reference object
reference <- Reference(counts, as.factor(cell_types), colSums(counts)) #Warning: Reference: nUMI does not match colSums of counts. If this is unintended, please correct this discrepancy. If this is intended, there is no problem.



################################################################## 

# Define a function to process each Seurat object
return_RCTD <- function(seurat_obj, reference, max_cores = 8) {
  coords <- GetTissueCoordinates(seurat_obj)
  colnames(coords) <- c("x", "y")
  
  # Remove columns with NA in column names
  coords <- coords[, !is.na(colnames(coords))]
  
  counts <- seurat_obj[["Spatial"]]$counts 
  query <- SpatialRNA(coords, counts, colSums(counts))
  myRCTD <- create.RCTD(query, reference, max_cores = max_cores)
  myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')
  
  return(myRCTD)
}

myRCTD_list = lapply (seurat_transformed, function(x) return_RCTD(x,reference))

saveRDS(myRCTD_list,"/omics/groups/OE0436/internal/Linh/Output/ST/myRCTD_list")