#ct_ann_GSE344

library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(celldex)
library(SingleR)
library(Seurat)
library (ggplot2)
library(stringr)
set.seed (0)


GSE_integrated = readRDS("/omics/groups/OE0436/internal/Linh/results/GSE_integrated")

# GSE_integrated subset for GSE344: 

GSE344 = subset (GSE_integrated, subset = dataset_ID == "GSE186344")
GSE344 =  GSE344@assays[["SCT"]]@data # extract the normalised count as well as genes and cell barcode 

GSE344_barcode = colnames (GSE344_matrix)

# import all the annotation files 
setwd("/omics/groups/OE0436/data/Linh/Datasets/GSE186344/")

folder_path = "/omics/groups/OE0436/data/Linh/Datasets/GSE186344/" #path to folder
files <- list.files(folder_path) #list all of the files in folder_path
celltype = data.frame ()

# create a dataframe that contains cell barcode and cell types infos
for (i  in 1:length (files)) {
  if (stringr::str_detect(files[[i]],"Cell_Types_Annotations") == TRUE) {
    celltype  <- rbind(celltype, read.csv(gzfile (files[[i]])))
  } 
}


###################### match cell barcode and cell type 

new_row2 = c(rep (NA, length (GSE344_barcode)))
extr_barcode = str_split(celltype [,1], "_")  #extract barcode info 
barcode= c()
barcode = sapply(extr_barcode, tail, n =1) #take the barcode info (last part of the list)


GSE344_split = str_split(GSE344_barcode, "_")  #extract barcode info 
GSE344_split = sapply(GSE344_split, tail, n =1) #take the barcode info (last part of the list)
GSE344_split = str_split(GSE344_split, "-") 
GSE344_split = sapply(GSE344_split, head, n =1) #take the barcode info (last part of the list)

# iterate over each barcode in vector A
for (i in 1:length(GSE344_split)) {
  col_index = which(grepl(GSE344_split[i],barcode) == TRUE) #extract col index
  if (length (col_index)>0) {
    new_row2[i] = celltype [col_index,2]
    celltype = celltype [-col_index[1],] #so that the cell barcode wont repeat
    barcode = barcode[-col_index[1]]
  }
}

write.table(new_row2,"/omics/groups/OE0436/data/Linh/Datasets/GSE186344/new_row2")
###########################################
new_row2 = read.table("/omics/groups/OE0436/data/Linh/Datasets/GSE186344/new_row2")

GSE344@Dimnames[[3]] = new_row2 [,1]

# check for same cell barcode

which(grepl("CAACTAGGTCCGAGTC",GSE344@Dimnames[[2]]) == TRUE) #437 43265 50624
GSE344@Dimnames[[3]]  [43265]
#[1] "MTC" : correct

# save GSE344 as ref dataset
writeMM (GSE344,"/omics/groups/OE0436/data/Linh/Datasets/Annotation/ref_GSE344")

