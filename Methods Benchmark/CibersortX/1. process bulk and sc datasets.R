## part2.1 process bulk and singlecell GEO dataasets
## library
library(Seurat)
source("../../dataPreprocess.r")

## bulk
## load bulk exp-matrix and clinical matrix (stage II and III)
clinic_list <- load_bulk_Clinical(path="different")
exp_matrix_list <- load_bulk_Exp(path="different")
names=names(exp_matrix_list)

## save the exp data
dir.create("./bulk_exp/")
for (bulkRNA in names) {
  exp_tem <- exp_matrix_list[[bulkRNA]]

  exp_tem <- as.data.frame(exp_tem)
  exp_tem[1:5, 1:5]

  exp_tem <- cbind(as.character(rownames(exp_tem)), exp_tem)
  exp_tem[1:5, 1:5]

  col_name <- colnames(exp_tem)
  exp_tem <- rbind(col_name, exp_tem)

  exp_tem[1:5, 1:5]

  exp_tem[1, 1] <- "Gene_Symbol"

  write.table(exp_tem,
    file = paste0("./bulk_exp/", bulkRNA, "_exp.txt"), row.names = F, col.names = F,
    quote = F, sep = "\t", fileEncoding = "utf-8"
  )
}

## sc
scRNA_seq_vector <- c("GSE144735", "GSE132465", "GSE132257")
for (scRNA in scRNA_seq_vector) {
  seurat_obj <- readRDS(paste0("../../result/", scRNA, "_SeuratObj_anno.rds"))
  seurat_obj <- subset(seurat_obj, (Class == "Tumor" | Class == "Border"))
  high_varia_genes <- seurat_obj@assays$RNA@var.features

  sc_exp <- seurat_obj@assays$RNA@counts[rownames(seurat_obj@assays$RNA@counts) %in% high_varia_genes, ]
  sc_exp <- exp(sc_exp) - 1
  sc_exp[1:5, 1:5]
  sc_exp <- as.data.frame(sc_exp)

  col_name <- c("GeneSymbol", seurat_obj@meta.data$Cell_subtype)
  row_name <- rownames(sc_exp)

  sc_exp <- cbind(row_name, sc_exp)
  sc_exp <- rbind(col_name, sc_exp)

  sc_exp[1:5, 1:5]

  write.table(sc_exp,
    file = paste0("./scRNA_mat/",scRNA, "_exp.txt"), row.names = F, col.names = F,
    sep = "\t", quote = F
  )
}
