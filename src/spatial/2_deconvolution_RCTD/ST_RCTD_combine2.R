
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(readxl)
library(spacexr)
library(Matrix)
library(doParallel)



ref = readRDS(".../Output/ST/sc_ref") #SCT, with annotation, joined layers (ccRCC)
seurat_transformed = readRDS(".../Datasets/Spatial_transcriptomics/Seurat_obj/seurat_transformed")

################################################################## offset the coords
# Offset for slice 2 (adjust the offset values as needed)
# Initialize with the first sample
slide1 = seurat_transformed[[1]]
slide2 = seurat_transformed[[2]]
slide3 = seurat_transformed[[3]]
slide4 = seurat_transformed[[4]]
slide5 = seurat_transformed[[5]]
slide6 = seurat_transformed[[6]]
slide7 = seurat_transformed[[7]]
slide8 = seurat_transformed[[8]]
slide9 = seurat_transformed[[9]]
slide10 = seurat_transformed[[10]]
slide11 = seurat_transformed[[11]]
slide12 = seurat_transformed[[12]]
slide13 = seurat_transformed[[13]]
slide14 = seurat_transformed[[14]]
slide15 = seurat_transformed[[15]]
slide16 = seurat_transformed[[16]]

slide2@images[["slice1"]]@coordinates[["imagerow"]] <- slide2@images[["slice1"]]@coordinates[["imagerow"]] + max(slide1@images[["slice1"]]@coordinates[["imagerow"]]) + 100 #slide2@images[["slice1"]] is correct no worry
slide2@images[["slice1"]]@coordinates[["imagecol"]] <- slide2@images[["slice1"]]@coordinates[["imagecol"]] + max(slide1@images[["slice1"]]@coordinates[["imagecol"]]) + 100 #slide2@images[["slice1"]] is correct no worry

slide3@images[["slice1"]]@coordinates[["imagerow"]] <- slide3@images[["slice1"]]@coordinates[["imagerow"]] + max(slide2@images[["slice1"]]@coordinates[["imagerow"]]) + 100 #slide2@images[["slice1"]] is correct no worry
slide3@images[["slice1"]]@coordinates[["imagecol"]] <- slide3@images[["slice1"]]@coordinates[["imagecol"]] + max(slide2@images[["slice1"]]@coordinates[["imagecol"]]) + 100 #slide2@images[["slice1"]] is correct no worry

slide4@images[["slice1"]]@coordinates[["imagerow"]] <- slide4@images[["slice1"]]@coordinates[["imagerow"]] + max(slide3@images[["slice1"]]@coordinates[["imagerow"]]) + 100 #slide2@images[["slice1"]] is correct no worry
slide4@images[["slice1"]]@coordinates[["imagecol"]] <- slide4@images[["slice1"]]@coordinates[["imagecol"]] + max(slide3@images[["slice1"]]@coordinates[["imagecol"]]) + 100 #slide2@images[["slice1"]] is correct no worry

slide5@images[["slice1"]]@coordinates[["imagerow"]] <- slide5@images[["slice1"]]@coordinates[["imagerow"]] + max(slide4@images[["slice1"]]@coordinates[["imagerow"]]) + 100 #slide2@images[["slice1"]] is correct no worry
slide5@images[["slice1"]]@coordinates[["imagecol"]] <- slide5@images[["slice1"]]@coordinates[["imagecol"]] + max(slide4@images[["slice1"]]@coordinates[["imagecol"]]) + 100 #slide2@images[["slice1"]] is correct no worry

slide6@images[["slice1"]]@coordinates[["imagerow"]] <- slide6@images[["slice1"]]@coordinates[["imagerow"]] + max(slide5@images[["slice1"]]@coordinates[["imagerow"]]) + 100 #slide2@images[["slice1"]] is correct no worry
slide6@images[["slice1"]]@coordinates[["imagecol"]] <- slide6@images[["slice1"]]@coordinates[["imagecol"]] + max(slide5@images[["slice1"]]@coordinates[["imagecol"]]) + 100 #slide2@images[["slice1"]] is correct no worry

slide7@images[["slice1"]]@coordinates[["imagerow"]] <- slide7@images[["slice1"]]@coordinates[["imagerow"]] + max(slide6@images[["slice1"]]@coordinates[["imagerow"]]) + 100 #slide2@images[["slice1"]] is correct no worry
slide7@images[["slice1"]]@coordinates[["imagecol"]] <- slide7@images[["slice1"]]@coordinates[["imagecol"]] + max(slide6@images[["slice1"]]@coordinates[["imagecol"]]) + 100 #slide2@images[["slice1"]] is correct no worry

slide8@images[["slice1"]]@coordinates[["imagerow"]] <- slide8@images[["slice1"]]@coordinates[["imagerow"]] + max(slide7@images[["slice1"]]@coordinates[["imagerow"]]) + 100 #slide2@images[["slice1"]] is correct no worry
slide8@images[["slice1"]]@coordinates[["imagecol"]] <- slide8@images[["slice1"]]@coordinates[["imagecol"]] + max(slide7@images[["slice1"]]@coordinates[["imagecol"]]) + 100 #slide2@images[["slice1"]] is correct no worry

slide9@images[["slice1"]]@coordinates[["imagerow"]] <- slide9@images[["slice1"]]@coordinates[["imagerow"]] + max(slide8@images[["slice1"]]@coordinates[["imagerow"]]) + 100 #slide2@images[["slice1"]] is correct no worry
slide9@images[["slice1"]]@coordinates[["imagecol"]] <- slide9@images[["slice1"]]@coordinates[["imagecol"]] + max(slide8@images[["slice1"]]@coordinates[["imagecol"]]) + 100 #slide2@images[["slice1"]] is correct no worry

slide10@images[["slice1"]]@coordinates[["imagerow"]] <- slide10@images[["slice1"]]@coordinates[["imagerow"]] + max(slide9@images[["slice1"]]@coordinates[["imagerow"]]) + 100 #slide2@images[["slice1"]] is correct no worry
slide10@images[["slice1"]]@coordinates[["imagecol"]] <- slide10@images[["slice1"]]@coordinates[["imagecol"]] + max(slide9@images[["slice1"]]@coordinates[["imagecol"]]) + 100 #slide2@images[["slice1"]] is correct no worry

slide11@images[["slice1"]]@coordinates[["imagerow"]] <- slide11@images[["slice1"]]@coordinates[["imagerow"]] + max(slide10@images[["slice1"]]@coordinates[["imagerow"]]) + 100 #slide2@images[["slice1"]] is correct no worry
slide11@images[["slice1"]]@coordinates[["imagecol"]] <- slide11@images[["slice1"]]@coordinates[["imagecol"]] + max(slide10@images[["slice1"]]@coordinates[["imagecol"]]) + 100 #slide2@images[["slice1"]] is correct no worry

slide12@images[["slice1"]]@coordinates[["imagerow"]] <- slide12@images[["slice1"]]@coordinates[["imagerow"]] + max(slide11@images[["slice1"]]@coordinates[["imagerow"]]) + 100 #slide2@images[["slice1"]] is correct no worry
slide12@images[["slice1"]]@coordinates[["imagecol"]] <- slide12@images[["slice1"]]@coordinates[["imagecol"]] + max(slide11@images[["slice1"]]@coordinates[["imagecol"]]) + 100 #slide2@images[["slice1"]] is correct no worry

slide13@images[["slice1"]]@coordinates[["imagerow"]] <- slide13@images[["slice1"]]@coordinates[["imagerow"]] + max(slide12@images[["slice1"]]@coordinates[["imagerow"]]) + 100 #slide2@images[["slice1"]] is correct no worry
slide13@images[["slice1"]]@coordinates[["imagecol"]] <- slide13@images[["slice1"]]@coordinates[["imagecol"]] + max(slide12@images[["slice1"]]@coordinates[["imagecol"]]) + 100 #slide2@images[["slice1"]] is correct no worry

slide14@images[["slice1"]]@coordinates[["imagerow"]] <- slide14@images[["slice1"]]@coordinates[["imagerow"]] + max(slide13@images[["slice1"]]@coordinates[["imagerow"]]) + 100 #slide2@images[["slice1"]] is correct no worry
slide14@images[["slice1"]]@coordinates[["imagecol"]] <- slide14@images[["slice1"]]@coordinates[["imagecol"]] + max(slide13@images[["slice1"]]@coordinates[["imagecol"]]) + 100 #slide2@images[["slice1"]] is correct no worry

slide15@images[["slice1"]]@coordinates[["imagerow"]] <- slide15@images[["slice1"]]@coordinates[["imagerow"]] + max(slide14@images[["slice1"]]@coordinates[["imagerow"]]) + 100 #slide2@images[["slice1"]] is correct no worry
slide15@images[["slice1"]]@coordinates[["imagecol"]] <- slide15@images[["slice1"]]@coordinates[["imagecol"]] + max(slide14@images[["slice1"]]@coordinates[["imagecol"]]) + 100 #slide2@images[["slice1"]] is correct no worry

slide16@images[["slice1"]]@coordinates[["imagerow"]] <- slide16@images[["slice1"]]@coordinates[["imagerow"]] + max(slide15@images[["slice1"]]@coordinates[["imagerow"]]) + 100 #slide2@images[["slice1"]] is correct no worry
slide16@images[["slice1"]]@coordinates[["imagecol"]] <- slide16@images[["slice1"]]@coordinates[["imagecol"]] + max(slide15@images[["slice1"]]@coordinates[["imagecol"]]) + 100 #slide2@images[["slice1"]] is correct no worry

merge_slice = merge(slide1,c(slide2,slide3,slide4,slide5,slide6,slide7,slide8,slide9,slide10,slide11,slide12,slide13,slide14,slide15,slide16))


################################################################## sc ref
counts = ref@assays[["SCT"]]@counts
counts@Dimnames[[1]] = Features(ref)
counts@Dimnames[[2]] = Cells(ref)
cell_types = ref@meta.data[["celltype1"]]
names(cell_types) = counts@Dimnames[[2]]
### Create the Reference object
reference <- Reference(counts, as.factor(cell_types), colSums(counts)) 


################################################################## RCTD
# coord of the first samples 
imagerow = c(merge_slice@images[["slice1"]]@coordinates[["imagerow"]]
             ,merge_slice@images[["slice1.2"]]@coordinates[["imagerow"]]
             ,merge_slice@images[["slice1.3"]]@coordinates[["imagerow"]]
             ,merge_slice@images[["slice1.4"]]@coordinates[["imagerow"]]
             ,merge_slice@images[["slice1.5"]]@coordinates[["imagerow"]]
             ,merge_slice@images[["slice1.6"]]@coordinates[["imagerow"]]
             ,merge_slice@images[["slice1.7"]]@coordinates[["imagerow"]]
             ,merge_slice@images[["slice1.8"]]@coordinates[["imagerow"]]
             ,merge_slice@images[["slice1.9"]]@coordinates[["imagerow"]]
             ,merge_slice@images[["slice1.10"]]@coordinates[["imagerow"]]
             ,merge_slice@images[["slice1.11"]]@coordinates[["imagerow"]]
             ,merge_slice@images[["slice1.12"]]@coordinates[["imagerow"]]
             ,merge_slice@images[["slice1.13"]]@coordinates[["imagerow"]]
             ,merge_slice@images[["slice1.14"]]@coordinates[["imagerow"]]
             ,merge_slice@images[["slice1.15"]]@coordinates[["imagerow"]]
             ,merge_slice@images[["slice1.16"]]@coordinates[["imagerow"]])

imagecol = c(merge_slice@images[["slice1"]]@coordinates[["imagecol"]]
             ,merge_slice@images[["slice1.2"]]@coordinates[["imagecol"]]
             ,merge_slice@images[["slice1.3"]]@coordinates[["imagecol"]]
             ,merge_slice@images[["slice1.4"]]@coordinates[["imagecol"]]
             ,merge_slice@images[["slice1.5"]]@coordinates[["imagecol"]]
             ,merge_slice@images[["slice1.6"]]@coordinates[["imagecol"]]
             ,merge_slice@images[["slice1.7"]]@coordinates[["imagecol"]]
             ,merge_slice@images[["slice1.8"]]@coordinates[["imagecol"]]
             ,merge_slice@images[["slice1.9"]]@coordinates[["imagecol"]]
             ,merge_slice@images[["slice1.10"]]@coordinates[["imagecol"]]
             ,merge_slice@images[["slice1.11"]]@coordinates[["imagecol"]]
             ,merge_slice@images[["slice1.12"]]@coordinates[["imagecol"]]
             ,merge_slice@images[["slice1.13"]]@coordinates[["imagecol"]]
             ,merge_slice@images[["slice1.14"]]@coordinates[["imagecol"]]
             ,merge_slice@images[["slice1.15"]]@coordinates[["imagecol"]]
             ,merge_slice@images[["slice1.16"]]@coordinates[["imagecol"]])

coords = data.frame(imagerow, imagecol)
colnames(coords) <- c("x", "y")
rownames(coords) = Cells(merge_slice)
# Remove columns with NA in column names
coords <- coords[, !is.na(colnames(coords))]
counts <- merge_slice[["SCT"]]$counts 
query <- SpatialRNA(coords, counts, colSums(counts)) 

myRCTD <- create.RCTD(query, reference, max_cores = 8)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

saveRDS(myRCTD,".../Output/ST/myRCTD_combine_try2")

