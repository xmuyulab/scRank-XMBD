## using FindMarker as features to find the prognostic-related cell-subtypes
## library
library(survival)
library(Seurat)
library(RColorBrewer)
library(ggplot2)

source("./0. functions.r")
source("../../dataPreprocess.r")

## load scRNAseq data(tumor and border)
seurat_list <- list()
scRNA_name_vector <- c("GSE144735", "GSE132465", "GSE132257")

for (scRNA_name in scRNA_name_vector) {
  seurat_tem <- readRDS(paste0("../../result/", scRNA_name, "_SeuratObj_anno.rds"))

  seurat_tem <- subset(seurat_tem, (Class == "Tumor" | Class == "Border"))
  seurat_list[[scRNA_name]] <- seurat_tem
  rm(seurat_tem)
}

## load exp matrix
exp_matrix_list <- load_bulk_Exp(path="different")
##  load clincal data
clinic_list <- load_bulk_Clinical(path="different")

validate_id <- c(1, 3, 9)

list_ <- split_Datasets(exp_matrix_list, clinic_list, validate_id)

all_exp <- list_[[1]]
all_clinical <- list_[[2]]
rm(list_)

## find markers of each cell-subtype in all scRNAseq datasets
markers_list <- list()
for (scRNA_name in scRNA_name_vector) {
  markers_list[[scRNA_name]] <- find_markers(seurat_list[[scRNA_name]],
    logFC_cutoff = 1, padj_cutoff = 0.2
  )
}

saveRDS(markers_list, file = "./result/markers_list.rds")
markers_list <- readRDS("./result/markers_list.rds")

## the genes in bulk RNAseq datasets
bulk_genes <- rownames(all_exp)

## find the overlap genes in scRNA and bulk RNA
for (scRNA_name in scRNA_name_vector) {
  markers_list[[scRNA_name]] <- find_overlap_scRNA_bulkRNA(
    markers_list[[scRNA_name]],
    bulk_genes
  )
}

## univariable cox for all markers
cox_result_list <- list()
for (scRNA_name in scRNA_name_vector) {
  cox_result_list[[scRNA_name]] <- univariable_cox_for_marker(markers_list[[scRNA_name]],
    all_exp, all_clinical,
    padj_cutoff = 0
  )
}

saveRDS(cox_result_list, "./result/cox_result_list.rds")

## visualize
## load major celltype names
major_celltype_df <- load_major_celltype_name()
for (scRNA_name in scRNA_name_vector) {
  plot_df <- as.data.frame(cox_result_list[[scRNA_name]])

  ## some transform
  plot_df$lnHR <- log(as.numeric(plot_df$HR))
  plot_df$Cell_type <- major_celltype_df[match(plot_df$Cell_subtype, major_celltype_df$cell_subtype), 2]
  plot_df <- plot_df[order(plot_df$Cell_type), ]

  ## plot
  major_celltype <- names(table(plot_df$Cell_type))
  plot_list <- list()
  for (i in 1:length(major_celltype)) {
    plot_df_tem <- plot_df[which(plot_df$Cell_type == major_celltype[i]), ]
    p_tem <- ggplot(data = plot_df_tem, aes(x = Cell_subtype, y = lnHR, fill = Cell_subtype)) +
      geom_boxplot(alpha = 0.7) +
      scale_y_continuous(name = "ln(HR)") +
      scale_x_discrete(name = paste0("Cell subtype of ", major_celltype[i])) +
      ggtitle(paste0("Prognostic value boxplot of cell-subtype of ", major_celltype[i])) +
      theme_bw() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 9, angle = 0)
      )

    plot_list[[major_celltype[i]]] <- p_tem
  }
  if (length(plot_list) == 7) {
    p <- (plot_list[[1]]) /
      (plot_list[[2]]) /
      (plot_list[[3]] + plot_list[[5]]) /
      (plot_list[[4]]) /
      (plot_list[[6]]) /
      (plot_list[[7]])
  }
  if (length(plot_list) == 6) {
    p <- (plot_list[[1]]) /
      (plot_list[[2]]) /
      (plot_list[[3]]) /
      (plot_list[[5]]) /
      (plot_list[[4]]) /
      (plot_list[[6]])
  }
  pdf(
    file = paste0("./figure/", scRNA_name, "_prognostic_markers_boxplot", ".pdf"),
    width = 18, height = 22
  )
  print(p)
  dev.off()
}
