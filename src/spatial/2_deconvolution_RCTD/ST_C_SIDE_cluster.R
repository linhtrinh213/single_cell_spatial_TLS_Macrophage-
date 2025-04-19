
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(readxl)
library(spacexr)
library(Matrix)
library(doParallel)


myRCTD_list = readRDS(".../Output/ST/myRCTD_list")
seurat_transformed = readRDS(".../Datasets/Spatial_transcriptomics/Seurat_obj/seurat_transformed")



CSIDE_input = function (seu_obj) {
  # Identifying TLS positive regions (output: logical vector)
  TLS = seu_obj$TLS_anno == "TLS"
  # preparing the explanatory variable (#The explanatory variable itself is a vector of values, constrained between 0 and 1, with names matching the pixel names of the myRCTD object.)
  pixel_names <- Cells(seu_obj)
  explanatory_variable <- rep(0, length(pixel_names))
  explanatory_variable[pixel_names %in% rownames(GetTissueCoordinates(seu_obj)[TLS, ])] <- 1  # Assign 1 to TLS positive regions
  names(explanatory_variable) <- pixel_names
  return(explanatory_variable)
}

explanatory_variable_list <- lapply(seurat_transformed, CSIDE_input) #apply the function onto seurat object lists 


RunRCTD = function(myRCTD, explanatory_variable) {
  myRCTD@config$max_cores <- 2
  myRCTD <- run.CSIDE.single(myRCTD, explanatory_variable, cell_types = "Macro", 
                               doublet_mode = FALSE, weight_threshold = 0, 
                               cell_type_threshold = 0)
    return(myRCTD)

}

myRCTD_list = mapply (RunRCTD, myRCTD_list,explanatory_variable_list, SIMPLIFY = F) #mapply is used to apply RunRCTD to each pair of myRCTD_list and explanatory_variable_list. SIMPLIFY = FALSE ensures the output is a list.
#SIMPLIFY = FALSE ensures the output is a list. 

saveRDS(myRCTD_list, ".../Output/ST/myRCTD_DEA")

############################################################ try with 0.25 weight threshold
myRCTD_list_DEA_new  = list()

#myRCTD_list_DEA_new [[1]] =  run.CSIDE.single(myRCTD_list [[1]], explanatory_variable_list[[1]], cell_types = "Macro", 
                                           #       doublet_mode = FALSE, weight_threshold = 0.25, 
                                            #      cell_type_threshold = 0) #it was 30 (cel) # only 15 genes turned out to be significant 
myRCTD_list_DEA_new [[1]] =  run.CSIDE.single(myRCTD_list [[1]], explanatory_variable_list[[1]], cell_types = "Macro", 
                                              doublet_mode = FALSE, weight_threshold = 0.15, 
                                              cell_type_threshold = 30) #seems better -> try with 0.15
myRCTD_list_DEA_new [[2]] =  run.CSIDE.single(myRCTD_list [[2]], explanatory_variable_list[[2]], cell_types = "Macro", 
                                              doublet_mode = FALSE, weight_threshold = 0.15, 
                                              cell_type_threshold = 30)
myRCTD_list_DEA_new [[3]] =  run.CSIDE.single(myRCTD_list [[3]], explanatory_variable_list[[3]], cell_types = "Macro", 
                                              doublet_mode = FALSE, weight_threshold = 0.15, 
                                              cell_type_threshold = 0)
myRCTD_list_DEA_new [[4]] =  run.CSIDE.single(myRCTD_list [[4]], explanatory_variable_list[[4]], cell_types = "Macro", 
                                              doublet_mode = FALSE, weight_threshold = 0.15, 
                                              cell_type_threshold = 0)
myRCTD_list_DEA_new [[5]] =  run.CSIDE.single(myRCTD_list [[5]], explanatory_variable_list[[5]], cell_types = "Macro", 
                                              doublet_mode = FALSE, weight_threshold = 0.15, 
                                              cell_type_threshold = 0)
myRCTD_list_DEA_new [[6]] =  run.CSIDE.single(myRCTD_list [[6]], explanatory_variable_list[[6]], cell_types = "Macro", 
                                              doublet_mode = FALSE, weight_threshold = 0.15, 
                                              cell_type_threshold = 0)
myRCTD_list_DEA_new [[7]] =  run.CSIDE.single(myRCTD_list [[7]], explanatory_variable_list[[7]], cell_types = "Macro", 
                                              doublet_mode = FALSE, weight_threshold = 0.15, 
                                              cell_type_threshold = 0) ########################################
myRCTD_list_DEA_new [[8]] =  run.CSIDE.single(myRCTD_list [[8]], explanatory_variable_list[[8]], cell_types = "Macro", 
                                              doublet_mode = FALSE, weight_threshold = 0.15, 
                                              cell_type_threshold = 0) #Error in find_sig_genes_individual(cell_type, cell_types, gene_fits, gene_list_type,  : 
#find_sig_genes_individual: cell type Macro has not converged on any genes. Consider removing this cell type from the model using the cell_types option
myRCTD_list_DEA_new [[9]] =  run.CSIDE.single(myRCTD_list [[9]], explanatory_variable_list[[9]], cell_types = "Macro", 
                                              doublet_mode = FALSE, weight_threshold = 0.15, 
                                              cell_type_threshold = 0)
myRCTD_list_DEA_new [[10]] =  run.CSIDE.single(myRCTD_list [[10]], explanatory_variable_list[[10]], cell_types = "Macro", 
                                              doublet_mode = FALSE, weight_threshold = 0.15, 
                                              cell_type_threshold = 0)
myRCTD_list_DEA_new [[11]] =  run.CSIDE.single(myRCTD_list [[11]], explanatory_variable_list[[11]], cell_types = "Macro", 
                                              doublet_mode = FALSE, weight_threshold = 0.15, 
                                              cell_type_threshold = 0)
myRCTD_list_DEA_new [[12]] =  run.CSIDE.single(myRCTD_list [[12]], explanatory_variable_list[[12]], cell_types = "Macro", 
                                              doublet_mode = FALSE, weight_threshold = 0.15, 
                                              cell_type_threshold = 0)
myRCTD_list_DEA_new [[13]] =  run.CSIDE.single(myRCTD_list [[13]], explanatory_variable_list[[13]], cell_types = "Macro", 
                                              doublet_mode = FALSE, weight_threshold = 0.15, 
                                              cell_type_threshold = 0)
myRCTD_list_DEA_new [[14]] =  run.CSIDE.single(myRCTD_list [[14]], explanatory_variable_list[[14]], cell_types = "Macro", 
                                              doublet_mode = FALSE, weight_threshold = 0.15, 
                                              cell_type_threshold = 0)
myRCTD_list_DEA_new [[15]] =  run.CSIDE.single(myRCTD_list [[15]], explanatory_variable_list[[15]], cell_types = "Macro", 
                                              doublet_mode = FALSE, weight_threshold = 0.15, 
                                              cell_type_threshold = 0)
myRCTD_list_DEA_new [[16]] =  run.CSIDE.single(myRCTD_list [[16]], explanatory_variable_list[[16]], cell_types = "Macro", 
                                              doublet_mode = FALSE, weight_threshold = 0.15, 
                                              cell_type_threshold = 0)

saveRDS(myRCTD_list, ".../Output/ST/myRCTD_DEA_0.25")
############################################################################### Inspecting results 
newthr = myRCTD_list_DEA_new[[1]]@de_results[["all_gene_list"]][["Macro"]]
newthr$gene_name = rownames(newthr)
