## part2.2 scissor-result visulization
library(Seurat)
library(Scissor)
library(clustree)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(ggpubr)
library(RColorBrewer)

source("../../dataPreprocess.r")

## there are three datasets were get the scissor results
scRNAseq_seurat_vector <- c("GSE132257", "GSE144735", "GSE132465") ## scRNA-seq datasets

exp_matrix_files <- list.files("../../data/bulkExp")
names <- sapply(strsplit(x = exp_matrix_files, split = "[.]"), function(x) x[1])
names <- sapply(strsplit(x = names, split = "_"), function(x) {
  x[1]
})
bulkRNAseq_seurat_vector <- names ## bulkRNA-seq datasets
rm(exp_matrix_files, names)

## load annotation seurat object
anno_seurat_list <- list()
for (scRNAseq_data in scRNAseq_seurat_vector) {
  anno_seurat_list[[scRNAseq_data]] <- readRDS(paste0("./data/", scRNAseq_data, "(Tumor)_seurat_obj.rds"))
}

## visulize each scRNAseq datasets' scissor results
for (scRNAseq_data in scRNAseq_seurat_vector) {
  bulkRNAseq_list <- list()

  ## attend scissor result to annotation results
  for (bulkRNAseq_data in bulkRNAseq_seurat_vector) {
    cat("scRNA-seq data: ", scRNAseq_data, "\n")
    cat("\t", "bulkRNA-seq data: ", bulkRNAseq_data, "\n")
    if (scRNAseq_data == "GSE132465") {
      all <- anno_seurat_list[["GSE132465"]]

      part_1 <- tryCatch(readRDS(paste0(
        "./result/",
        scRNAseq_data, "_1_", bulkRNAseq_data, "_scissor_seurat.rds"
      )),
      error = function(e) {
        return(NULL)
      },
      warning = function(w) {
        return(NULL)
      }
      )

      part_2 <- tryCatch(readRDS(paste0(
        "./result/",
        scRNAseq_data, "_2_", bulkRNAseq_data, "_scissor_seurat.rds"
      )),
      error = function(e) {
        return(NULL)
      },
      warning = function(w) {
        return(NULL)
      }
      )

      if (length(part_1) == 0 && length(part_2) == 0) {
        next
      }
      if (length(part_1) == 0 && length(part_2) != 0) {
        index2 <- match(rownames(part_2@meta.data), rownames(all@meta.data))

        all@meta.data$scissor <- 0
        all@meta.data$scissor[index2] <- part_2@meta.data$scissor

        bulkRNAseq_list[[bulkRNAseq_data]] <- all
      }
      if (length(part_1) != 0 && length(part_2) == 0) {
        index1 <- match(rownames(part_1@meta.data), rownames(all@meta.data))

        all@meta.data$scissor <- 0
        all@meta.data$scissor[index1] <- part_1@meta.data$scissor

        bulkRNAseq_list[[bulkRNAseq_data]] <- all
      }
      if (length(part_1) != 0 && length(part_2) != 0) {
        index1 <- match(rownames(part_1@meta.data), rownames(all@meta.data))
        index2 <- match(rownames(part_2@meta.data), rownames(all@meta.data))

        all@meta.data$scissor <- NA
        all@meta.data$scissor[index1] <- part_1@meta.data$scissor
        all@meta.data$scissor[index2] <- part_2@meta.data$scissor

        bulkRNAseq_list[[bulkRNAseq_data]] <- all
      }
      rm(part_1, part_2, index2, index1)
    } else {
      all <- anno_seurat_list[[scRNAseq_data]]

      tem <- tryCatch(readRDS(paste0(
        "./result/",
        scRNAseq_data, "_", bulkRNAseq_data, "_scissor_seurat.rds"
      )),
      error = function(e) {
        return(1)
      },
      warning = function(w) {
        return(1)
      }
      )
      if (length(names(tem)) == 0) {
        next
      }

      index <- match(rownames(tem@meta.data), rownames(all@meta.data))

      all@meta.data$scissor <- NA
      all@meta.data$scissor[index] <- tem@meta.data$scissor

      bulkRNAseq_list[[bulkRNAseq_data]] <- all
      rm(all, index, tem)
    }
  }

  ## extract dataframe to visualize
  scissor_result_matrix <- matrix(ncol = 6, nrow = 1)
  colnames(scissor_result_matrix) <- c("Cell_id", "Scissor_label", "Cluster_id", "Cell_type", "Cell_subtype", "BulkRNA_id")
  for (z in 1:length(bulkRNAseq_list)) {
    cat("Visualizations", "\n")
    cat("The dataset is ", names(bulkRNAseq_list)[z], "\n")

    ## class of scissor label
    cat(class(bulkRNAseq_list[[z]]@meta.data$scissor), "\n")

    ## table scissor label
    table(bulkRNAseq_list[[z]]@meta.data$scissor)

    ## transform to a matrix
    scissor_tem_matrix <- matrix(data = NA, ncol = 6, nrow = length(bulkRNAseq_list[[z]]@meta.data$scissor))
    colnames(scissor_tem_matrix) <- c("Cell_id", "Scissor_label", "Cluster_id", "Cell_type", "Cell_subtype", "BulkRNA_id")
    scissor_tem_matrix <- as.data.frame(scissor_tem_matrix)

    ## each column
    scissor_tem_matrix$Cell_id <- rownames(bulkRNAseq_list[[z]]@meta.data)
    scissor_tem_matrix$Scissor_label <- bulkRNAseq_list[[z]]@meta.data$scissor
    scissor_tem_matrix$Cluster_id <- bulkRNAseq_list[[z]]@meta.data$cell_typeid
    scissor_tem_matrix$Cell_type <- bulkRNAseq_list[[z]]@meta.data$Cell_type
    scissor_tem_matrix$Cell_subtype <- bulkRNAseq_list[[z]]@meta.data$Cell_subtype
    scissor_tem_matrix$BulkRNA_id <- rep(names(bulkRNAseq_list)[z], nrow(bulkRNAseq_list[[z]]@meta.data))

    scissor_result_matrix <- rbind(scissor_result_matrix, scissor_tem_matrix)
  }
  scissor_result_matrix <- scissor_result_matrix[-1, ]

  scissor_result_matrix1 <- scissor_result_matrix[which(scissor_result_matrix$Scissor_label == 1), ]
  scissor_result_matrix2 <- scissor_result_matrix[which(scissor_result_matrix$Scissor_label == 2), ]

  ## plot
  p1 <- ggplot(data = scissor_result_matrix1, mapping = aes(x = Cluster_id, fill = BulkRNA_id)) +
    geom_bar(stat = "count", width = 0.5, position = "dodge") + ## count表示条形图的高度是数量
    guides(fill = guide_legend(reverse = T)) +
    scale_fill_manual(values = brewer.pal(length(table(scissor_result_matrix1$BulkRNA_id)), "Set3")) +
    geom_text(stat = "count", aes(label = ..count..), color = "black", size = 3.5, position = position_dodge(0.5), vjust = -0.5) +
    labs(
      x = "Cluster ID", y = "Scissor label",
      title = paste0("Scissor+ cells of prognostic-risk(without adjuvant-chemotherapy) of ", scRNAseq_data, " Cluster id.pdf")
    ) +
    theme_classic() +
    coord_flip()


  p2 <- ggplot(data = scissor_result_matrix2, mapping = aes(x = Cluster_id, fill = BulkRNA_id)) +
    geom_bar(stat = "count", width = 0.5, position = "dodge") + ## count表示条形图的高度是数量
    guides(fill = guide_legend(reverse = T)) +
    scale_fill_manual(values = brewer.pal(length(table(scissor_result_matrix2$BulkRNA_id)), "Set3")) +
    geom_text(stat = "count", aes(label = ..count..), color = "black", size = 3.5, position = position_dodge(0.5), vjust = -0.5) +
    labs(x = "Cluster ID", y = "Scissor label", title = paste0("Scissor+ cells of prognostic-protected(without adjuvant-chemotherapy) of ", scRNAseq_data, " Cluter id.pdf")) +
    theme_classic() +
    coord_flip()

  p3 <- ggplot(data = scissor_result_matrix1, mapping = aes(x = Cell_type, fill = BulkRNA_id)) +
    geom_bar(stat = "count", width = 0.5, position = "dodge") + ## count表示条形图的高度是数量
    guides(fill = guide_legend(reverse = T)) +
    scale_fill_manual(values = brewer.pal(length(table(scissor_result_matrix1$BulkRNA_id)), "Set3")) +
    geom_text(stat = "count", aes(label = ..count..), color = "black", size = 3.5, position = position_dodge(0.5), vjust = -0.5) +
    labs(
      x = "Cell type", y = "Scissor label",
      title = paste0("Scissor+ cells of prognostic-risk(without adjuvant-chemotherapy) of ", scRNAseq_data, " Cell type.pdf")
    ) +
    theme_classic() +
    coord_flip()


  p4 <- ggplot(data = scissor_result_matrix2, mapping = aes(x = Cell_type, fill = BulkRNA_id)) +
    geom_bar(stat = "count", width = 0.5, position = "dodge") + ## count表示条形图的高度是数量
    guides(fill = guide_legend(reverse = T)) +
    scale_fill_manual(values = brewer.pal(length(table(scissor_result_matrix2$BulkRNA_id)), "Set3")) +
    geom_text(stat = "count", aes(label = ..count..), color = "black", size = 3.5, position = position_dodge(0.5), vjust = -0.5) +
    labs(x = "Cell type", y = "Scissor label", title = paste0("Scissor+ cells of prognostic-protected of ", scRNAseq_data, " Cell type.pdf")) +
    theme_classic() +
    coord_flip()

  p5 <- ggplot(data = scissor_result_matrix1, mapping = aes(x = Cell_subtype, fill = BulkRNA_id)) +
    geom_bar(stat = "count", width = 0.5, position = "dodge") + ## count表示条形图的高度是数量
    guides(fill = guide_legend(reverse = T)) +
    scale_fill_manual(values = brewer.pal(length(table(scissor_result_matrix1$BulkRNA_id)), "Set3")) +
    geom_text(stat = "count", aes(label = ..count..), color = "black", size = 3.5, position = position_dodge(0.5), vjust = -0.5) +
    labs(
      x = "Cell type", y = "Scissor label",
      title = paste0("Scissor+ cells of prognostic-risk of ", scRNAseq_data, " Cell subtype.pdf")
    ) +
    theme_classic() +
    coord_flip()

  p6 <- ggplot(data = scissor_result_matrix2, mapping = aes(x = Cell_subtype, fill = BulkRNA_id)) +
    geom_bar(stat = "count", width = 0.5, position = "dodge") + ## count表示条形图的高度是数量
    guides(fill = guide_legend(reverse = T)) +
    scale_fill_manual(values = brewer.pal(length(table(scissor_result_matrix2$BulkRNA_id)), "Set3")) +
    geom_text(stat = "count", aes(label = ..count..), color = "black", size = 3.5, position = position_dodge(0.5), vjust = -0.5) +
    labs(x = "Cell type", y = "Scissor label", title = paste0("Scissor+ cells of prognostic-protected of ", scRNAseq_data, " Cell subtype.pdf")) +
    theme_classic() +
    coord_flip()

  pdf(file = paste0("./figures/barplots/", scRNAseq_data, " Scissor results distribution Cluster ID.pdf"), width = 12, height = 9)
  print(p1 + p2)
  dev.off()

  pdf(file = paste0("./figures/barplots/", scRNAseq_data, " Scissor results distribution Cell type.pdf"), width = 12, height = 9)
  print(p3 + p4)
  dev.off()

  pdf(file = paste0("./figures/barplots/", scRNAseq_data, " Scissor results distribution Cell subtype.pdf"), width = 12, height = 9)
  print(p5 + p6)
  dev.off()
}
