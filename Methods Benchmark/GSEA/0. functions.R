## part3 functions

## load datasets
load_exp_list=function(){
  ## load exp matrix
  exp_matrix_files = list.files('../../data/bulkExp')
  names = sapply(strsplit(x = exp_matrix_files, split = '[.]'), function(x) x[1])
  names = sapply(strsplit(x = names, split = '_'), function(x) {x[1]})
  exp_matrix_list = list()
  for (i in 1:length(exp_matrix_files)) {
    exp_matrix_list[[names[i]]] = readRDS(paste0('../../data/bulkExp/',exp_matrix_files[i]))
  }
  return(exp_matrix_list)
}

load_clinical_list=function(){
  ##  load clincal data
  clinic_list = list()
  clinical_files = list.files('../../data/bulkClinical/')
  names=gsub('_','\\.',clinical_files)
  names=sapply(names,function(x){strsplit(x,'\\.')[[1]][2]})
  names=as.character(names)
  for (i in 1:length(clinical_files)) {
    clinic_list[[names[i]]] = read.table(paste0('../../data/bulkClinical/',
                                                clinical_files[i]), sep = "\t", header = T,stringsAsFactors = F,row.names = 1)}
  return(clinic_list)
}

cox_model=function(exp_df,clinical_df){
  if(as.numeric(table(colnames(exp_df)%in%rownames(clinical_df)))!=nrow(clinical_df)){
    cat("Patients' ID not match")
    return(NULL)
  }
  
  else{
    if(ncol(clinical_df)==3){clinical_df=clinical_df[,-1]}
    exp_df=exp_df[,match(rownames(clinical_df),colnames(exp_df))]
    if(class(exp_df)=="matrix"){exp_df=as.data.frame(exp_df)}
    res_df=as.data.frame(matrix(NA, nrow = 1, ncol = 5))
    colnames(res_df)=c("Gene_symbol","p value","HR","HR confint lower","HR confint upper")
    
    for(i in 1:nrow(exp_df)){
      ## generated dataframe
      gene_tem=exp_df[i,]
      tem_df=cbind(clinical_df,t(gene_tem))
      colnames(tem_df)=c("time","event","gene")
      
      ## cox
      tem_cox=tryCatch(coxph(Surv(time,event)~gene,tem_df),
                       error=function(x){
                         return(1)})
      if(length(tem_cox)==1){next;}
      x=summary(tem_cox)
      
      # get p-value
      p_value=signif(x$wald["pvalue"], digits=2)
      
      # get HR
      HR = signif(x$coef[2], digits=2);
      # get 95% confidence interval
      HR_confint_lower <- signif(x$conf.int[,"lower .95"],2)
      HR_confint_upper <- signif(x$conf.int[,"upper .95"],2)
      
      res<-c(rownames(exp_df)[i],p_value,HR,HR_confint_lower,HR_confint_upper)
      names(res)<-c("Gene_symbol","p value","HR","HR confint lower","HR confint upper")
      
      res_df=rbind(res_df,res)
    }
    res_df=res_df[-1,]
    rownames(res_df)=res_df$Gene_symbol
    res_df=res_df[,-1]
    for(j in ncol(res_df)){
      res_df[,j]=as.numeric(res_df[,j])
    }
    res_df$FDR=p.adjust(res_df$`p value`,method = "BH")
    return(res_df)
  }
}

## select PDR minor 0.1
fileter_cox_results=function(cox_result,cutoff){
  if(cutoff==0){
    index_significant=which(cox_result$`p value`<=0.05)
  }
  else{
    index_significant=which(cox_result$FDR<=cutoff)
  }
  cox_result=cox_result[index_significant,]
  return(cox_result)
}

## load seuart data and extract tumor tissue exp matrix
load_scRNA_seq_seurat=function(scRNA_data){
  seurat_obj=readRDS(paste0("../../result/",scRNA_data,"_SeuratObj_anno.rds"))
  
  seurat_obj=subset(seurat_obj,(Class=="Tumor" | Class=='Border'))
  
  exp=as.data.frame(seurat_obj@assays$RNA@counts) 
  print(dim(exp))
  
  Celltype=seurat_obj@meta.data$Cell_type
  Cellsubtype=seurat_obj@meta.data$Cell_subtype
  
  exp_anno_list=list()
  exp_anno_list[["exp"]]=exp
  exp_anno_list[["Celltype"]]=Celltype
  exp_anno_list[["Cellsubtype"]]=Cellsubtype
  
  return(exp_anno_list)
}

## calculate the mean expression of each celltype and cell subtype
calculate_mean_expression=function(tumor_exp_list){
  exp=tumor_exp_list[["exp"]]
  Celltype_df=as.data.frame(table(tumor_exp_list[["Celltype"]]))
  colnames(Celltype_df)=c("Celltype_name","Freq")
  
  Cellsubtype_df=as.data.frame(table(tumor_exp_list[["Cellsubtype"]]))
  colnames(Cellsubtype_df)=c("Cellsubtype_name","Freq")
  
  Cell_type_mean_exp=matrix(data = NA,nrow = nrow(exp),ncol = nrow(Celltype_df))
  for(i in 1:nrow(Celltype_df)){
    
    exp_tem=exp[,which(tumor_exp_list[["Celltype"]]==as.character(Celltype_df[i,1]))]
    cat("The cell type is ",as.character(Celltype_df[i,1])," ",'\n')
    
    ## for cell number less than 1
    if(class(exp_tem)=="numeric"){
      cat("This celltype's cell number less than 2",'\n')
      exp_tem_mean=exp_tem
      Cell_type_mean_exp[,i]=exp_tem_mean
    }
    
    ## for cell number more than 1
    else{ 
      cat("The dim of exp of this celltype is ",dim(exp_tem),'\n')
      exp_tem_mean=apply(exp_tem, 1, function(x){mean(x)})
      Cell_type_mean_exp[,i]=exp_tem_mean
    }
  }
  colnames(Cell_type_mean_exp)=Celltype_df[,1]
  rownames(Cell_type_mean_exp)=rownames(exp)
  
  print(head(Cell_type_mean_exp))
  
  Cell_subtype_mean_exp=matrix(data = NA,nrow = nrow(exp),ncol = nrow(Cellsubtype_df))
  for(i in 1:nrow(Cellsubtype_df)){
    
    exp_tem=exp[,which(tumor_exp_list[["Cellsubtype"]]==as.character(Cellsubtype_df[i,1]))]
    cat("The cell subtype is ",as.character(Cellsubtype_df[i,1])," ",'\n')
    
    ## for cell number less than 1
    if(class(exp_tem)=="numeric"){
      cat("This cellsubtype's cell number less than 2",'\n')
      exp_tem_mean=exp_tem
      Cell_subtype_mean_exp[,i]=exp_tem_mean
      
    }
    
    ## for cell number more than 1
    else{ 
      cat("The dim of exp of this cellsubtype is ",dim(exp_tem),'\n')
      exp_tem_mean=apply(exp_tem, 1, function(x){mean(x)})
      Cell_subtype_mean_exp[,i]=exp_tem_mean
    }
  }
  colnames(Cell_subtype_mean_exp)=Cellsubtype_df[,1]
  rownames(Cell_subtype_mean_exp)=rownames(exp)
  
  print(head(Cell_subtype_mean_exp))
  
  mean_exp_list=list()
  mean_exp_list[["Celltype_mean_exp"]]=Cell_type_mean_exp
  mean_exp_list[["Cell_subtype_mean_exp"]]=Cell_subtype_mean_exp
  
  return(mean_exp_list)
}

## GSEA for cox
calculate_GSEA_score_cox=function(cox_result,scRNA_mean_exp_list,prognostic){
  
    scRNA_mean_exp_tem=scRNA_mean_exp_list
    
    ## GSEA-Main
    ssGSEA_result=matrix(data = NA,nrow = 0,ncol = ncol(scRNA_mean_exp_list))
    
    for(j in 1:length(cox_result)){
      cox_df=cox_result[[j]]
      ## if cox_df==0，return Null
      if(nrow(cox_df)<5){
        next
      }
      
      ## get HR>1(risk-related) gene from COX
      if(prognostic=="risk"){
        risk_cox_genes=cox_df[which(cox_df$HR>1),]
      }
      if(prognostic=="protective"){
        risk_cox_genes=cox_df[which(cox_df$HR<1),]
      }
      
      ## make gene list
      risk_gene_list=list(rownames(risk_cox_genes))
      names(risk_gene_list)=names(cox_result)[j] ## bulk datasets'geo_accession as name
      
      
      fit=try(
        ssGSEA_res <-gsva(as.matrix(scRNA_mean_exp_tem),risk_gene_list,
                          method = "ssgsea",kcdf = "Gaussian",abs.ranking = T)
      )
      if(class(fit)=="try-error"){
        ssGSEA_result=rbind(ssGSEA_result,rep(NA,ncol(scRNA_mean_exp_list[[i]])))
      }
      
      else{
        ssGSEA_result=rbind(ssGSEA_result,ssGSEA_res)
      }
    }
  return(ssGSEA_result)
}

## heatmap for GSEA-cox
heatmap_visualization=function(GSEA_cox_results_df){
  ## make same major celltype together
  if(1){
    B_cells=c("CD19+ CD20+ Plasma B cell","Cycling B cell","IgA+ IGLC2+ Plasma B cell",
              "IgA+ IGLL5+ Plasma B cell","IgA+ Plasma B cell","IgG+ Plasma B cell",
              "Naive B cell") ## 7
    
    T_cells=c("CD4+ IL7R+ T cells","CD4+ Tfh","CD4+ Th17","CD8+ Exhausted-CTL",
              "CD8+ GZMK+ CTL","CD8+ KLRD1+ CTL","Gamma-delta T cells","Naive T cells",
              "Regulatory T cells","NK cells") ## 10
    
    Myeloid_cells=c("cDC_CD1C","Macro_C1QC","Macro_INHBA",
                    "Macro_DNASE1L3","Macro_LYVE1","Macro_SPP1",
                    "Mono_FCN1","pDC_GZMB") ## 8
    
    Fibroblast_cells=c("eCAF","Fibro_ADAMDEC1","Fibro_BMP5","Fibro_GPM6B",
                       "Fibro_ID1","Fibro_LY6H","Fibro_OGN",
                       "Fibro_SGK1","iCAF","mCAF","myCAF_DES",
                       "myCAF_NOTCH3","Pericyte") ## 13
    
    Mast_cells=c("Mast cells") ## 1
    
    Endothelial_cells=c("EC_ACKR1","EC_ESM1","EC_GPIHBP1","EC_IGFBP3",
                        "EC_PLVAP_CD36","EC_PLVAP_PRCP",
                        "EC_TFF3","EC_Unspecific") ## 8
    
    Epithelial_cells=c("CMS1","CMS2","CMS3","CMS4") ## 4
  }
  
  cellsubtype_order=c(B_cells,Endothelial_cells,Epithelial_cells,Fibroblast_cells,Mast_cells,Myeloid_cells,T_cells)
  
  GSEA_cox_results_df=GSEA_cox_results_df[,match(cellsubtype_order,colnames(GSEA_cox_results_df))]
  index=!(is.na(GSEA_cox_results_df[1,]))
  GSEA_cox_results_df=GSEA_cox_results_df[,index]
  
  Cell_type=c(rep("B cell",length(B_cells)), 
              rep("Endothelial cells",length(Endothelial_cells)), 
              rep("Epithelial cell",length(Epithelial_cells)),
              rep("Fibroblasts",length(Fibroblast_cells)), 
              rep("Mast",length(Mast_cells)), 
              rep("Myeloid cells",length(Myeloid_cells)), 
              rep("T cell",length(T_cells))
              )
  Cell_type=Cell_type[index]
  Cell_type=as.data.frame(Cell_type)
  rownames(Cell_type)=colnames(GSEA_cox_results_df)
  
  gaps_col_vector=as.numeric(table(Cell_type))
  
  for(i in 2:(length(gaps_col_vector))){
    gaps_col_vector[i]=gaps_col_vector[i]+gaps_col_vector[i-1]
  }
  
  ## heatmap
  p=pheatmap(as.data.frame(GSEA_cox_results_df),
             show_rownames = T,show_colnames = T,
             scale = "row",cluster_cols = F,cluster_rows = F,legend = T,
             gaps_col = c(gaps_col_vector[-length(gaps_col_vector)]),
             annotation_col = Cell_type,annotation_names_col = T,
             display_numbers=TRUE,cellwidth = 22, cellheight = 14,
             color = colorRampPalette(rev(brewer.pal(n=length(table(Cell_type[,1])),name = "RdYlBu")))(100))
  
  return(p)
  
}
