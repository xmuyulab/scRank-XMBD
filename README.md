# scRank_XMBD
## Introduction
We proposed a data analysis framework to prioritize prognostic-associated subpopulations based on relative expression orderings (REOs). Cell type specific gene pairs (C-GPs) were identified to evaluate prognostic value for each cell type. Individualized recurrence risk signatures at single-cell resolution were developed based on REOs. The results shown that REOs-based signatures could classify accurately among most cell subtypes. C-GPs achieves higher precision compared with four current methods. Moreover, we developed single-cell gene pair signatures (scGPSs) to predict recurrence risk for patients individually.<br>
<br>

![](https://github.com/Lin-Yux/scRank_XMBD/blob/non-master/workflow.png)
## Docker
We provide our R package conda enviroment images in [docker hub](https://hub.docker.com/)<br>
You can use docker in Linux to run the code directly. You can start with the following command.<br>
```
docker pull watermelontreesjs/scrank_xmbd
docker run --name scRank_XMBD -itd watermelontreesjs/scrank_xmbd
docker exec -it scRank_XMBD /bin/bash
cd ./root/code/
conda activate scRank_XMBD
R
```
You can exit docker environment by 
```
exit
```
## Usage
In ```main.r``` , we provide the main steps for scRank construction. Including:<br>
### data availability
The example data stored in ```data``` and ```result```.<br>
```
source("./dataPreprocess.r")
source("./utils.r")
source("./visulize.r")
source("./model.r")
set.seed(619)
```
### Identification of cell type specific gene pairs with cell subpopulation-classifying value.
#### Applying SingleCellNet to bulid GP classifier
```
data <- readRDS("./result/GSE144735_SeuratObj_anno.rds")

# parameter setting 
ncells <- 10
nTopGenes <- 100
nTopGenePairs <- 250

# use tumor and border tissue
ClassList <- c("Tumor", "Border")

exp_matrix <- Interface_Seurat_SCN(data, ClassList)[["exp_matrix"]]
anno_matrix <- Interface_Seurat_SCN(data, ClassList)[["anno_matrix"]]
training_SCN_classifier(exp_matrix,anno_matrix, "GSE144735", nTopGenes, nTopGenePairs, ncells)
```
#### get cell subtype specific gene pair(C-GPs)
```
# load exp and clinical data
exp_matrix_list <- load_bulk_Exp()
clinic_list <- load_bulk_Clinical()
list_ <- combine_Datasets(exp_matrix_list, clinic_list)
training_exp <- list_[[1]]
training_clinical <- list_[[2]]

# parameter setting 
ncells <- 10
deltaS <- 0.6

average_exp_of_top_genepairs <- readRDS("./result/GSE144735(ncells=10)top_genepairs.rds"")
specific_genepairs_list <- get_C_Gps(average_exp_of_top_genepairs, "GSE144735", deltaS, ncells)

# visualize the celltype-specific gene pairs
visualize_celltype_specific_genepairs(average_exp_of_top_genepairs, specific_genepairs_list)
```
### Evaluation of prognostic value for each cell type.
```
prognostic_CGPsList <- get_C_Gps_with_prognosis(training_exp,training_clinical,specific_genepairs_list, scRNA_name, ncells, deltaS, clinical_cutoff)
```
### Development of individualized recurrence risk signatures.
#### build clinical signature
```
# load exp and clinical data
exp_matrix_list <- load_bulk_Exp()
clinic_list <- load_bulk_Clinical()
raw_clinical <- load_bulk_RawClinical()
validate_id <- c(1, 3, 9)
list_ <- split_Datasets(exp_matrix_list, clinic_list, validate_id)

# prepare training and test list
all_exp <- list_[[1]]
all_clinical <- list_[[2]]
training_size <- list_[[3]]
test_size <- list_[[4]]

# parameter setting 
scRNA_name_vector <- c("GSE144735", "GSE132465", "GSE132257")
iteration_times <- 1000
ncells <- 10
deltaS <- 0.6

# combine stable C-GPs from different scRNA-seq dataset
stable_pairs_list <- combine_StableCGPs(scRNA_name_vector, ncells, deltaS)
celltype_list <- names(stable_pairs_list)
validation_sets <- c("GSE14333", "GSE17536", "GSE39582")

# applied lasso-cox model to bulid signature
LassoCox_signature(stable_pairs_list, celltype_list, iteration_times, all_exp, all_clinical, training_size, test_size, validation_sets)
```
#### select the iteration successful
```
path <- "./model/"
iteration_result <- analysis_KMs(ncells, celltype_list, path)

# Lollipop chart
Lollipop_chart(iteration_result)
```
#### multiple variables cox
```
celltype_list <- list.files("./model/")

for (celltype in celltype_list) {
  path_ <- paste0("./model/", celltype)
  vali_clinical <- MultiCox_Data_Transform(path_, validation_sets, exp_matrix_list, clinic_list, raw_clinical)
 
  if(vali_clinical=="next"){next;}

  MultiCox(vali_clinical)
}
```
In ```REO Stablity in Bulk``` , we evaluate the REO of some gene-pairs in bulk RNA-seq dataset with ground truth.<br>
<br>
In ```Methods Benchmark``` , we evaluate performance of five methods (scRank, CibersortX, GSEA, Scissor, Uni-Markers) used for prioritizing prognostic-associated subpopulations.<br>

## How to cite
我也不知道
