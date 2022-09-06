## part2.2 using bubble plot to visualize the results of CibersortX
## library
library(survival)
library(ggplot2)
library(RColorBrewer)
source("./0. functions.R")

## load bulk abundence with each singlecell-RNAseq datasets and bulk RNAseq clinical data
GSEID_vector <- c("GSE144735", "GSE132465", "GSE132257")
for (scRNA in GSEID_vector) {
  fraction_list <- load_celltype_fraction(scRNA)[["fraction_list"]]
  clinic_list <- load_celltype_fraction(scRNA)[["clinic_list"]]

  ## univariable cox
  cox_list <- univariable_cox(fraction_list, clinic_list)

  ## xlab
  n_samples <- as.vector(sapply(clinic_list, function(x) nrow(x)))
  names_sub <- paste0(names(fraction_list), "_", n_samples)
  names_sub[1:5]

  cell_type_names <- colnames(fraction_list[[1]])[1:(ncol(fraction_list[[1]]) - 3)]
  for (i in 1:length(cell_type_names)) {
    cell_type_names[i] <- gsub(x = cell_type_names[i], pattern = "\\.\\.", replacement = "+ ")
    cell_type_names[i] <- gsub(x = cell_type_names[i], pattern = "\\.", replacement = " ")
  }

  bubble_plot <- data.frame(
    cellType = rep(cell_type_names, length(fraction_list)),
    dataSet = rep(names_sub, each = (ncol(fraction_list[[1]]) - 3)),
    HR = as.numeric(unlist(lapply(cox_list, function(x) x$HR))),
    p.value = as.numeric(unlist(lapply(cox_list, function(x) x$p.value))),
    stringsAsFactors = F
  )

  load_cellsubtype_vector()

  bubble_plot$cellType <- factor(bubble_plot$cellType,
    levels = rev(c(B_cells, T_cells, Myeloid_cells, Fibroblast_cells, Mast_cells, Endothelial_cells, Epithelial_cells))
  )
  bubble_plot$dataSet <- factor(bubble_plot$dataSet,
    levels = names_sub
  )

  ## bubble style
  bubble_plot$HR <- ifelse(bubble_plot$HR > 1, "Risk", "Protective")
  bubble_plot$p.value[bubble_plot$p.value < 0.001] <- 10
  bubble_plot$p.value[bubble_plot$p.value < 0.01] <- 8
  bubble_plot$p.value[bubble_plot$p.value < 0.05] <- 6
  bubble_plot$p.value[bubble_plot$p.value < 0.1] <- 4
  bubble_plot$p.value[bubble_plot$p.value >= 0.1 & bubble_plot$p.value <= 1] <- 2

  bubble_plot <- na.omit(bubble_plot)

  colors_y <- c(
    rep("#8DD3C7", 7), rep("#CCEBC5", 8), rep("#BEBADA", 7),
    rep("#80B1D3", 13), rep("#FDB462", 1), rep("#B3DE69", 8), rep("#BC80BD", 4)
  )
  colors_y <- rev(colors_y)

  p <- ggplot(data = bubble_plot, aes(x = dataSet, y = cellType, size = p.value, color = HR)) +
    scale_color_manual(values = c("#359830", "#bf0603")) +
    geom_point() +
    theme_bw(base_size = 25) +
    theme(
      axis.text.x = element_text(size = 15, face = "bold", vjust = 1, hjust = 1, angle = 45),
      axis.text.y = element_text(size = 15, color = colors_y, face = "bold"),
      axis.title = element_blank()
    )

  pdf(file = paste0("./figures/", scRNA, "_CibersortX_result.pdf"), height = 12, width = 10)
  print(p)
  dev.off()
}
