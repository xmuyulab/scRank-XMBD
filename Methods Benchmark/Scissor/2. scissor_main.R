## part2.1 scissor-main
library(Seurat)
library(Scissor)
library(clustree)
library(ggplot2)
library(dplyr)

source("../../dataPreprocess.r")
source("./0. functions.R")

set.seed(619)

## load bulk clinical and expression matrix
clinic_list=load_bulk_Clinical(path="different")
exp_matrix_list=load_bulk_Exp(path="different")
names=names(exp_matrix_list)

## Scissor
for (bulk_name in names) {
  cat(bulk_name,'\n')
  if(bulk_name=="GSE29621"){next;}
  
  ## load bulkRNA-seq datasets
  exp=exp_matrix_list[[bulk_name]]
  clinical=clinic_list[[bulk_name]]
  if(length(table(colnames(exp)%in%rownames(clinical)))==1){  
    cat("The bulk dataset: ",bulk_name,"\t","colnames of exp are equal to clinical rownames.",'\n')}
  
  ## set the order is identity between two matrices
  exp=exp[,match(rownames(clinical),colnames(exp))]
  exp=as.matrix(exp)
  
  ## phenotype
  phenotype=clinical
  colnames(phenotype)=c("time", "status")
  rm(clinical)
  
  scRNA_gseid_vector=c("GSE144735","GSE132257","GSE132465_2","GSE132465_1")
  for(scRNA_gseid in scRNA_gseid_vector){
    cat("The bulkRNA-seq dataset be processed in this iteration is ",bulk_name,"\n")
    cat("The scRNA-seq dataset be processed in this iteration is ",scRNA_gseid,"\n")
    ## load scRNA-seq dataset with cluster id
    seurat=load_scRNAseq_dataset(scRNA_gseid)
    
    print(names(seurat))  ## check the network information we need
    #DimPlot(seurat, reduction = 'tsne', label = T, label.size = 10)
    
    infos4=Scissor(exp, seurat, phenotype , alpha = 0.2,
                   family = "cox", Save_file = paste0("./Scissor_obj/","Scissor_",scRNA_gseid,"_",bulk_name,".RData"))
    
    ## visualize result
    Scissor_select <- rep(0, ncol(seurat))
    names(Scissor_select) <- colnames(seurat)
    Scissor_select[infos4$Scissor_pos] <- 1
    Scissor_select[infos4$Scissor_neg] <- 2
    seurat <- AddMetaData(seurat, metadata = Scissor_select, col.name = "scissor")
    if(length(table(seurat@meta.data$scissor))==1){next;}
    pdf(file = paste0("./figures/Scissor_Umap/","Scissor_",scRNA_gseid,"_",bulk_name,"(alpha = 0.2).pdf"),width = 10,height = 7.5)
    p=DimPlot(seurat, reduction = 'tsne', group.by = 'scissor', 
              cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))
    print(p)
    dev.off()
    
    saveRDS(seurat,file = paste0("./result/",scRNA_gseid,"_",bulk_name,"_scissor_seurat.rds"))
  }
}

