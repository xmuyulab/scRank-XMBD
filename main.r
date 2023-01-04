# main

## library
if (T) {
  library(Seurat)
  library(clustree)
  library(ggplot2)
  library(dplyr)
  library(Biobase)
  library(CMScaller)
  library(GSVA)
  library(clusterProfiler)
  library(patchwork)
  library(ggpubr)
  library(ggtext)
  library(stringr)
  library(RColorBrewer)
  library(pheatmap)
  library(ggrepel)
  library(VennDiagram)
  library(cowplot)
  library(rstatix)
  library(glmnet)
  library(kableExtra)
  library(survival)
  library(singleCellNet)
  library(CMSclassifier)
  library(igraph)
  library(tableone)
  library(forestplot)
  library(tidyverse)
  library(scibet)
  library(viridis)
  library(survminer)
  library(maxstat)
  library(ggsci)
  library(ggthemes)
  library(org.Hs.eg.db)
  library(scibet)
  library(MetBrewer)
  source("./dataPreprocess.r")
  source("./utils.r")
  source("./visulize.r")
  source("./model.r")
  set.seed(619)
}
  
## part 1 - annotation scRNA-seq dataset

### pre-processing

GSE144735_annoFn()
GSE132465_preProcessFn()
GSE132257_preProcessFn()

### Scibet annotation query dataset (GSE144735 as reference)
GSE144735_seurat <- readRDS("./result/GSE144735_SeuratObj_anno.rds")
GSE144735_ref <- SciBetRef_transfrom(GSE144735_seurat)

load_and_predict_query(GSE144735_ref, "GSE132257")
load_and_predict_query(GSE144735_ref, "GSE132465")

readRDS("./result/GSE132257_SeuratObj_anno.rds")
rm(GSE144735_ref)

### cell subtype annotation
GSE132257_seurat <- readRDS("./result/GSE132257_SeuratObj_anno.rds")
GSE144735_seurat <- readRDS("./result/GSE144735_SeuratObj_anno.rds")
GSE132465_seurat <- readRDS("./result/GSE132465_SeuratObj_anno.rds")

rm(GSE144735_ref)

### cell subtype annotation
GSE132257_seurat <- readRDS("./result/GSE132257_SeuratObj.rds")
GSE144735_seurat <- readRDS("./result/GSE144735_SeuratObj.rds")
GSE132465_seurat <- readRDS("./result/GSE132465_SeuratObj.rds")

query_list <- list()
query_list[["GSE132257_seurat"]] <- GSE132257_seurat
query_list[["GSE132465_seurat"]] <- GSE132465_seurat

## view the celltype(major)
celltype_label <- GSE144735_seurat@meta.data$Cell_type
cat(table(celltype_label), "\n")
celltype_label_list <- c("B cells", "Epithelial cells", "Myeloid cells", "Fibroblasts", "T cells", "Endothelial cells")

predict_query_subtype(GSE144735_seurat, query_list, celltype_label_list)

## reload scRNA-seq dataset

datasets <- c("GSE144735", "GSE132257", "GSE132465")
dataset_list <- list()
for (i in datasets) {
  dataset_list[[i]] <- readRDS(paste0("./result/", i, "_SeuratObj_anno.rds"))
}
## plot the UMAP of major celltype (cell subtype) and cell-population fraction for each scRNA-seq data
plot_MarkersAndCellfraction(dataset_list)

## part 2 - Applying SingleCellNet to bulid GP classifier
datasets <- c("GSE144735", "GSE132257", "GSE132465")

### setting parameters
ClassList <- c("Tumor", "Border")

ncells <- 10
nTopGenes <- 100
nTopGenePairs <- 250

for (i in 1:length(dataset_list)) {
  data <- dataset_list[[i]]
  dataset_id <- datasets[i]

  ### the tumor and border tissue were considered to use this workflow
  exp_matrix <- Interface_Seurat_SCN(data, ClassList)[["exp_matrix"]]
  anno_matrix <- Interface_Seurat_SCN(data, ClassList)[["anno_matrix"]]
  training_SCN_classifier(exp_matrix,anno_matrix, dataset_id, nTopGenes, nTopGenePairs, ncells)

  ## After the above processing, we will get a classifer (under ./result) and
  ## a the PR curves for assessing the performance all cell-subpopulations (under ./figure)

  rm(exp_matrix,anno_matrix)
}


## part 3 - identifying C-GPs for each cell subpopulation
exp_matrix_list <- load_bulk_Exp()
clinic_list <- load_bulk_Clinical()

list_ <- combine_Datasets(exp_matrix_list, clinic_list)
training_exp <- list_[[1]]
training_clinical <- list_[[2]]
rm(list_)

scRNA_name_vector <- c("GSE144735", "GSE132465", "GSE132257")
ncells <- 10
deltaS <- 0.6
clinical_cutoff <- 0.1

for (scRNA_name in scRNA_name_vector) {
  path <- paste0("./result/", scRNA_name, "(ncells=", ncells, ")top_genepairs.rds")
  average_exp_of_top_genepairs <- readRDS(path)

  specific_genepairs_list <- get_C_Gps(average_exp_of_top_genepairs, scRNA_name, deltaS, ncells)
  ## visualize the celltype-specific gene pairs
  p <- visualize_celltype_specific_genepairs(
    average_exp_of_top_genepairs, specific_genepairs_list
  )
  png(filename = paste0("./figure/GPs average expression(", scRNA_name, ",ncells=", ncells, ",deltaS=", deltaS, ").png"), height = 800, width = 800)
  print(p)
  dev.off()

  prognostic_CGPsList <- get_C_Gps_with_prognosis(training_exp,training_clinical,specific_genepairs_list, scRNA_name, ncells, deltaS, clinical_cutoff)

  ## plotting the interesting cell sub-populations
  print(names(prognostic_CGPsList))
  CellSubTypeList <- c(
    "CD4+ Tfh", "CD8+ GZMK+ CTL", "Regulatory T cells", "IgA+ IGLC2+ Plasma B cell", "IgA+ IGLL5+ Plasma B cell",
    "Macro_SPP1", "Macro_C1QC", "eCAF", "myCAF_DES", "Fibro_SGK1"
  )
  
  major_celltype_df <- load_major_celltype_name()
  boxplotForprognotic_CGPs(prognostic_CGPsList, major_celltype_df, CellSubTypeList, scRNA_name, ncells, deltaS, clinical_cutoff)
}

## part 4 - the REO of some C-GPs
scRNA_name_vector <- c("GSE144735", "GSE132465")
scRNA_seurat_list <- list()
for (scRNA_name in scRNA_name_vector) {
  scRNA_seurat_list[[scRNA_name]] <- readRDS(paste0(
    "./result/", scRNA_name, "_SeuratObj_anno.rds"
  ))
  scRNA_seurat_list[[scRNA_name]] <- subset(scRNA_seurat_list[[scRNA_name]], (Class == "Tumor" | Class == "Border"))
}

list_ <- get_stable_C_Gps(scRNA_name_vector, scRNA_seurat_list, ncells, deltaS)
stable_CGPs <- list_[[1]]
consistance_list <- list_[[2]]
celltype_names <- list_[[3]]
rm(list_)

plot_n_CGPs(stable_CGPs, scRNA_seurat_list, consistance_list, celltype_names)

### meanwhile, you could plot certain C-GP in different scRNA-seq dataset
### visualize the certain gene pair rank in different cell subtype
dir.create("./figure/Certain_genepair")
setwd("./figure/Certain_genepair")

cellsubtypes <- c("IgA+ Plasma B cell", "Cycling B cell", "CD8+ KLRD1+ CTL", "Naive T cells", "CD8+ GZMK+ CTL")
visualize_certain_gene_pair("B2M", "IGKC", cellsubtypes, scRNA_seurat_list)
visualize_certain_gene_pair("CCL5", "GADD45B", cellsubtypes, scRNA_seurat_list)
setwd("../../")

## part 5 - build clinical signature
exp_matrix_list <- load_bulk_Exp()
clinic_list <- load_bulk_Clinical()

### validation sets
validate_id <- c(1, 3, 8)

list_ <- split_Datasets(exp_matrix_list, clinic_list, validate_id)
all_exp <- list_[[1]]
all_clinical <- list_[[2]]
training_size <- list_[[3]]
test_size <- list_[[4]]
rm(list_)

### setting parameters
scRNA_name_vector <- c("GSE144735", "GSE132465", "GSE132257")
iteration_times <- 1000
ncells <- 10
deltaS <- 0.6


stable_pairs_list <- combine_StableCGPs(scRNA_name_vector, ncells, deltaS)
celltype_list <- names(stable_pairs_list)

validation_sets <- c("GSE14333", "GSE17536", "GSE39582")
LassoCox_signature(stable_pairs_list, celltype_list, iteration_times, all_exp, all_clinical, training_size, test_size, validation_sets)

### to select the iteration successful
path <- "./model/"
iteration_result <- analysis_KMs(ncells, celltype_list, path)

### Lollipop chart
dir.create("./figure/model")
p <- Lollipop_chart(iteration_result)
pdf(file = paste0("./figure/model/Cell subtypes recurrence-risk prediction-value, ncells=", ncells, ".pdf"), width = 10, height = 15)
print(p)
dev.off()

### multiple variables cox
### for this analysis, we highly recommend to re-write this function to get the approciate results
exp_matrix_list <- load_bulk_Exp()
raw_clinical <- load_bulk_RawClinical()
clinic_list <- load_bulk_Clinical()

celltype_list <- list.files("./model/")
validation_sets <- c("GSE14333", "GSE17536", "GSE39582")

for (celltype in celltype_list) {
  path_ <- paste0("./model/", celltype)
  MultiCox_Data_Transform(path_, validation_sets, exp_matrix_list, clinic_list, raw_clinical)
}

### output the C-GPs of signatures
## you have to assigned the celltype and iteration id
if (F) {
  celltype <- c("IgA+ IGLC2+ Plasma B cell", "Fibro_SGK1", "CD8+ KLRD1+ CTL")
  iteration <- c(39, 327, 826)
  for (i in 1:3) {
    save_signature(celltype[i], iteration[i])
  }
}

