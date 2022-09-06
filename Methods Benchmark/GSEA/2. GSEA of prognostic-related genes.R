## part 3.2 GSEA and cox
## library
source("./0. functions.R")
library(Seurat)
library(GSVA)
library(limma) 
library(dplyr)
library(ggplot2)
library(ggsci)
library(gridExtra)
library(ggpubr)
library(pheatmap)
library(clusterProfiler)
library(enrichplot)
library(RColorBrewer)

## extract mean exp matrix
cox_result=readRDS("./CRC_bulk_RNAseq_cox_list.rds")

for(scRNA_seq in c("GSE144735","GSE132465","GSE132257")){
  tumor_exp_list=load_scRNA_seq_seurat(scRNA_seq)
  
  ## calculate the mean expression of each celltype and cell subtype
  mean_exp_list=calculate_mean_expression(tumor_exp_list)
  saveRDS(mean_exp_list,file=paste0(scRNA_seq,"_mean_exp_for_GSEA.rds"))
  
  ## GSEA-cox of each cellsubtype
  mean_exp_list=readRDS(paste0("./",scRNA_seq,"_mean_exp_for_GSEA.rds"))
  
  GSEA_cox_results_risk=calculate_GSEA_score_cox(cox_result,mean_exp_list[[2]],prognostic="risk")
  GSEA_cox_results_protect=calculate_GSEA_score_cox(cox_result,mean_exp_list[[2]],prognostic="protective")
  
  ## heatmap visualize
  pdf(file = paste0("./figures/",scRNA_seq," risk-related cox genes.pdf"),width = 20,height = 6)
  p=heatmap_visualization(GSEA_cox_results_risk)
  print(p)
  dev.off()
  
  pdf(file = paste0("./figures/",scRNA_seq," protect-related cox genes.pdf"),width = 20,height = 6)
  p=heatmap_visualization(GSEA_cox_results_protect)
  print(p)
  dev.off()
  
  delta_mat=GSEA_cox_results_risk-GSEA_cox_results_protect
  pdf(file = paste0("./figures/",scRNA_seq," risk-protect_mat cox genes.pdf"),width = 20,height = 6)
  p=heatmap_visualization(delta_mat)
  print(p)
  dev.off()
}

