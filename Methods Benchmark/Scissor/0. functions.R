## scissor-functions
## Load scRNA-seq dataset
load_scRNAseq_dataset <- function(GSEID) {
  path <- paste0(paste0(
    "./data/", GSEID, "(Tumor)_seurat_obj.rds"
  ))
  tem <- readRDS(path)

  return(tem)
}

