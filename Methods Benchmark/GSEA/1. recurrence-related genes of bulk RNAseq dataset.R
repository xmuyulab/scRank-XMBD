## part3.1 using GSEA to find the recurrence-related genes
## library
library(survival)
source("./0. functions.R")

## load exp data and clinical data (stage II and III)
exp_matrix_list <- load_exp_list()
clinic_matrix_list <- load_clinical_list()

## geoaccession
geo_accession <- names(clinic_matrix_list)

## cox
cox_result_list <- list()
for (i in 1:length(geo_accession)) {
  cox_result_list[[geo_accession[i]]] <- cox_model(exp_matrix_list[[i]], clinic_matrix_list[[i]])
}

## fileter significant
for (i in 1:length(geo_accession)) {
  cox_result_list[[i]] <- fileter_cox_results(cox_result_list[[i]], cutoff = 0)
}

## save cox result
saveRDS(cox_result_list, file = "CRC_bulk_RNAseq_cox_list.rds")
