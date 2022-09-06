## part 2.0 functions
## load bulk fraction and clinical files
load_celltype_fraction=function(scRNA){
  return_list=list()
  path=paste0("./cibersortX_results/",scRNA,"/")
  
  fraction_files = list.files(path)
  names=sapply(fraction_files,function(x){strsplit(x,split = '_')[[1]][2]})
  
  fraction_list = list()
  for (i in 1:length(fraction_files)) {
    fraction_list[[names[i]]] = read.table(paste0(path,fraction_files[i]), sep = "\t", header = T,stringsAsFactors = F,row.names = 1)
  }
  
  clinic_list = load_bulk_Clinical(path="different")
  rm(i)
  return_list[["fraction_list"]]=fraction_list
  return_list[["clinic_list"]]=clinic_list
  
  return(return_list)
}

## cox
univariable_cox=function(fraction_list,clinic_list){
  cox_list = list()
  names=names(fraction_list)
  for (i in 1:length(fraction_list)) {
    fraction = fraction_list[[i]]
    clinic = clinic_list[[i]]
    cox_mat = matrix(data=NA, nrow = (ncol(fraction)-3), ncol = 3)
    for (j in 1:(ncol(fraction)-3)) {
      label = fraction[,j]
      clinic$label = label
      coxmodel <- coxph(Surv(clinic$dfs_time, dfs_event) ~ label, data=clinic)
      coxmodel <- summary(coxmodel)
      cox_mat[j,] = c(coxmodel$concordance[1],coxmodel$coefficients[2],coxmodel$logtest["pvalue"])
    }
    
    cox_mat = data.frame(cellType = colnames(fraction)[1:(ncol(fraction)-3)], 
                         Cindex = cox_mat[,1],
                         HR = cox_mat[,2],
                         p.value = cox_mat[,3],
                         stringsAsFactors = F)
    cox_list[[names[i]]] = cox_mat
  }
  return(cox_list)
}

## load cellsubtype and celltype vector
load_cellsubtype_vector=function(){
  B_cells <<- c("CD19+ CD20+ Plasma B cell","Cycling B cell","IgA+ IGLC2+ Plasma B cell",
            "IgA+ IGLL5+ Plasma B cell","IgA+ Plasma B cell","IgG+ Plasma B cell",
            "Naive B cell") ## 7
  
  T_cells <<- c("CD4+ IL7R+ T cells","CD4+ Tfh","CD4+ Th17","CD8+ Exhausted-CTL",
            "CD8+ GZMK+ CTL","CD8+ KLRD1+ CTL","Gamma-delta T cells","Naive T cells",
            "Regulatory T cells","NK cells") ## 10
  
  Myeloid_cells <<- c("cDC_CD1C","Macro_C1QC","Macro_INHBA",
                  "Macro_DNASE1L3","Macro_LYVE1","Macro_SPP1",
                  "Mono_FCN1","pDC_GZMB") ## 8
  
  Fibroblast_cells <<- c("eCAF","Fibro_ADAMDEC1","Fibro_BMP5","Fibro_GPM6B",
                     "Fibro_ID1","Fibro_LY6H","Fibro_OGN",
                     "Fibro_SGK1","iCAF","mCAF","myCAF_DES",
                     "myCAF_NOTCH3","Pericyte") ## 13
  
  Mast_cells <<- c("Mast cells") ## 1
  
  Endothelial_cells <<- c("EC_ACKR1","EC_ESM1","EC_GPIHBP1","EC_IGFBP3",
                      "EC_PLVAP_CD36","EC_PLVAP_PRCP",
                      "EC_TFF3","EC_Unspecific") ## 8
  
  Epithelial_cells <<- c("CMS1","CMS2","CMS3","CMS4") ## 4
}
