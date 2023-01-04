## plot cell fraction and markers in annotation Seurat object
plot_MarkersAndCellfraction <- function(dataset_list) {

  ## define colors
  colors <- c(
    brewer.pal(n = 99, name = "Set1"), brewer.pal(n = 99, name = "Set2"), brewer.pal(n = 99, name = "Set3"),
    brewer.pal(n = 99, name = "Pastel1"), brewer.pal(n = 99, name = "Pastel2"), brewer.pal(n = 99, name = "Paired"), brewer.pal(n = 99, name = "YlGn")
  )
  ## major colors
  celltypes <- names(table(dataset_list[[1]]@meta.data$Cell_type))
  for (i in 2:length(dataset_list)) {
    celltypes <- union(celltypes, dataset_list[[i]]@meta.data$Cell_type)
  }

  colors_df_major <- data.frame("colors" = colors[1:length(celltypes)], "Cell_type" = celltypes)

  ## subtype colors
  cellsubtypes <- names(table(dataset_list[[1]]@meta.data$Cell_subtype))
  for (i in 2:length(dataset_list)) {
    cellsubtypes <- union(cellsubtypes, dataset_list[[i]]@meta.data$Cell_subtype)
  }

  colors_df_subtype <- data.frame("colors" = colors[1:length(cellsubtypes)], "Cell_subtype" = cellsubtypes)

  ## annotation results
  for (i in datasets) {
    anno_seurat <- dataset_list[[i]]

    ## cell type makers and annotation
    Idents(anno_seurat) <- anno_seurat@meta.data$Cell_type

    ### cell type makers
    pdf(file = paste0("./figure/annotation/", i, " Major Celltype makers FeaturePlot.pdf"), width = 10, height = 7.5)
    p <- FeaturePlot(anno_seurat, features = c(
      "KIT", "CPA3", ## Mast cells
      "CD79A", "CD79B", ## B cells
      "CD3D", "CD3G", ## T cells
      "CD68", ## Myeloid cells
      "EPCAM", ## Epithelial cells
      "THY1", "COL3A1", ## Fibroblasts
      "PECAM1" ## Endothelial cells
    ))
    print(p)
    dev.off()

    ### cell type annotation
    pdf(file = paste0("./figure/annotation/", i, " Major Celltype Umap.pdf"), width = 15, height = 7.5)
    p1 <- DimPlot(anno_seurat, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 5, combine = T, group.by = "Cell_type")
    p2 <- DimPlot(anno_seurat, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 5, combine = T, group.by = "Patient")
    print(p1 + p2)
    dev.off()

    ### cell subtype annotation
    Idents(anno_seurat) <- anno_seurat@meta.data$Cell_subtype
    pdf(file = paste0("./figure/annotation/", i, " Cell subtype Umap.pdf"), width = 20, height = 7.5)
    p1 <- DimPlot(anno_seurat, reduction = "umap", label = F, pt.size = 0.5, label.size = 5, combine = T, group.by = "Cell_subtype")
    p2 <- DimPlot(anno_seurat, reduction = "umap", label = F, pt.size = 0.5, label.size = 5, combine = T, group.by = "Patient")
    print(p1 + p2)
    dev.off()

    ### cell type fraction
    p <- cell_subtype_percent_bar_ploting(anno_seurat, type = "major", colors_df_major)
    pdf(file = paste0("./figure/annotation/", i, " Major Celltype fraction.pdf"), width = 8, height = 6)
    print(p)
    dev.off()

    ### cell subtype fraction
    p <- cell_subtype_percent_bar_ploting(anno_seurat, type = "subtype", colors_df_subtype)
    pdf(file = paste0("./figure/annotation/", i, " Cell subtype fraction.pdf"), width = 12, height = 6)
    print(p)
    dev.off()
  }
  return(NULL)
}

## the cell subtype fractions in each tissue(barplot)
cell_subtype_percent_bar_ploting <- function(seurat_object, type, colors) {

  ### extract annotation information
  plot_df <- seurat_object@meta.data[, c("Class", "Cell_type", "Cell_subtype")]

  ### plot bar plot of major cell type
  if (type == "major") {
    p <- ggplot(data = plot_df, mapping = aes(x = Class, fill = Cell_type)) +
      geom_bar(stat = "count", width = 0.5, position = "fill") +
      guides(fill = guide_legend(reverse = T)) +
      scale_fill_manual(values = colors[match(names(table(plot_df$Cell_type)), colors[, 2]), 1]) +
      labs(x = "Tissues", y = "Celltype Proportion", title = "") +
      theme_classic() +
      coord_flip()

    return(p)
  }

  ### bar plot of cell subtype
  if (type == "subtype") {
    p <- ggplot(data = plot_df, mapping = aes(x = Class, fill = Cell_subtype)) +
      geom_bar(stat = "count", width = 0.5, position = "fill") +
      guides(fill = guide_legend(reverse = T)) +
      scale_fill_manual(values = colors[match(names(table(plot_df$Cell_subtype)), colors[, 2]), 1]) +
      labs(x = "Tissues", y = "Cell Subtype Proportion", title = "") +
      theme_classic() +
      coord_flip()

    return(p)
  }
}

## visualize the celltype-specific gene pairs
visualize_celltype_specific_genepairs <- function(average_exp_of_top_genepairs, specific_genepairs_list) {
  all_celltype_specific_genepairs <- c()
  for (i in 1:length(specific_genepairs_list)) {
    all_celltype_specific_genepairs <- c(
      all_celltype_specific_genepairs,
      names(specific_genepairs_list[[i]])
    )
  }
  all_celltype_specific_genepairs <- unique(all_celltype_specific_genepairs)

  print(table(all_celltype_specific_genepairs %in% rownames(average_exp_of_top_genepairs)))
  plot_df <- average_exp_of_top_genepairs[match(all_celltype_specific_genepairs, rownames(average_exp_of_top_genepairs)), ]

  colnames(plot_df) <- sapply(colnames(plot_df), function(x) {
    strsplit(x, split = "_Avg")[[1]][1]
  })

  color <- colorRampPalette(c("#436eee", "white", "#EE0000"))(100)
  p <- pheatmap(plot_df,
    color = color, scale = "row",
    cluster_rows = F, cluster_cols = F,
    legend_labels = c("ΔS high", "ΔS low"), legend = T,
    show_rownames = F, show_colnames = T
  )
  return(p)
}

## boxplot of prognostic gene-pairs in each cell subtype
boxplotForprognotic_CGPs <- function(prognostic_specific_genepairs_list, major_celltype_df, CellSubTypeList, scRNA_name, ncells, deltaS, clinical_cutoff) {

  ## boxplot of prognostic gene-pairs in each cell subtype
  HR_ <- c()
  cellsubtype_ <- c()
  for (i in 1:length(prognostic_specific_genepairs_list)) {
    HR_tem <- as.numeric(prognostic_specific_genepairs_list[[i]]$genepair_prognostic_info$HR)
    celltype_name_tem <- rep(as.character(names(prognostic_specific_genepairs_list)[[i]]), length(HR_tem))

    HR_ <- c(HR_, HR_tem)
    cellsubtype_ <- c(cellsubtype_, celltype_name_tem)
  }

  plot_df <- as.data.frame(matrix(data = NA, nrow = length(HR_), ncol = 2))
  colnames(plot_df) <- c("HR", "Cell_subtype")

  plot_df$HR <- log(HR_)
  for (i in 1:nrow(plot_df)) {
    if (plot_df[i, "HR"] > 1.25) {
      plot_df[i, "HR"] <- 1.25
    }
    if (plot_df[i, "HR"] < (-3)) {
      plot_df[i, "HR"] <- (-3)
    }
  }

  plot_df$Cell_subtype <- cellsubtype_
  plot_df$Cell_type <- major_celltype_df[match(plot_df$Cell_subtype, major_celltype_df$cell_subtype), 2]

  plot_df <- plot_df[order(plot_df$Cell_type), ]

  ## plot
  major_celltype <- names(table(plot_df$Cell_type))

  plot_df_ <- plot_df[plot_df$Cell_subtype %in% CellSubTypeList, ]

  p_tem <- ggplot(data = plot_df_, aes(x = Cell_subtype, y = HR, fill = Cell_subtype)) +
    geom_boxplot(alpha = 0.7) +
    scale_y_continuous(name = "ln(HR)") +
    scale_x_discrete(name = paste0("Cell subtypes")) +
    ggtitle(paste0("Prognostic value boxplot of cell subpopulations")) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      text = element_text(size = 12),
      axis.title = element_text(face = "bold"),
      axis.text.x = element_text(size = 9, angle = 10)
    )
  pdf(paste0("./figure/The prognostic value of cell-subtype-specific gene pairs(", scRNA_name, ",ncells=", ncells, ",deltaS=", deltaS, ",cutoff=", clinical_cutoff, ").pdf"), width = 15, height = 6)
  print(p_tem)
  dev.off()
  return(NULL)
}

## visualize some C-GPs for certain cell subpopulation
plot_n_CGPs <- function(specific_gene_pairs_list, scRNA_seurat_list, consistance_list, celltype_names, n_pairs = 5, n_cells = 30) {
  system("mkdir ./figure/CGPs_REO")
  for (celltype in celltype_names) {
    setwd("./figure/CGPs_REO")
    dir.create(paste0("./", celltype))
    setwd(paste0("./", celltype))
    visualize_specific_pairs_rank(specific_gene_pairs_list, scRNA_seurat_list, celltype, consistance_list,
      n_pairs = n_pairs, n_cells = n_cells
    )
  setwd("../../../")
  }
  return(NULL)
}

## visualize some classic pairs of cell-subtype
visualize_specific_pairs_rank <- function(specific_gene_pairs_list, scRNA_seurat_list, celltype, consistance_list,
                                          n_pairs = 15, n_cells = 30) {
  set.seed(619)
  consistance_df <- consistance_list[[celltype]]

  if (nrow(consistance_df) == 0) {
    cat(paste0("Cell-subtype: ", celltype, " do not have stable gene pairs.", "\n"))
    return(NULL)
  }

  ## omit the specific gene pair list which dataset without celltype
  for (i in length(specific_gene_pairs_list):1) {
    if (!celltype %in% names(table(names(specific_gene_pairs_list[[i]])))) {
      specific_gene_pairs_list[[i]] <- NULL
    }
  }
  scRNA_name_vector <- names(specific_gene_pairs_list)

  ## gene pairs to be ploted
  if (n_pairs > nrow(consistance_df)) {
    n_pairs <- nrow(consistance_df)
  }
  top_n_pairs <- rownames(consistance_df[1:n_pairs, ])
  top_n_pairs_df <- matrix(NA, nrow = n_pairs, ncol = 2)

  colnames(top_n_pairs_df) <- c("gene1", "gene2")
  for (i in 1:nrow(top_n_pairs_df)) {
    top_n_pairs_df[i, ] <- c(strsplit(top_n_pairs[i], "_")[[1]][1], strsplit(top_n_pairs[i], "_")[[1]][2])
  }

  plot_list <- list()
  for (i in 1:nrow(top_n_pairs_df)) {
    plot_list_subset <- list()
    for (scRNA_name in scRNA_name_vector) {
      ## exp matrix of each scRNAseq dataset
      exp <- subset(scRNA_seurat_list[[scRNA_name]], Cell_subtype == celltype)@assays$RNA@counts
      if (n_cells > ncol(exp)) {
        n_cells <- ncol(exp)
      }
      ## sample cells to plot
      sample_cell_index <- sample(c(1:ncol(exp)), replace = F, n_cells)
      exp <- exp[union(top_n_pairs_df[, 1], top_n_pairs_df[, 2]), sample_cell_index]

      gene1 <- top_n_pairs_df[i, 1]
      gene2 <- top_n_pairs_df[i, 2]
      plot_df <- exp[union(gene1, gene2), ]

      plot_df2 <- matrix(NA, nrow = n_cells * 2, ncol = 5)
      colnames(plot_df2) <- c("cellid", "expression", "gene", "dataset", "pair")

      plot_df2[, 1] <- c(rep(1:n_cells, 2))
      plot_df2[, 2] <- c(plot_df[1, ], plot_df[2, ])
      plot_df2[, 3] <- c(rep(rownames(plot_df)[1], n_cells), rep(rownames(plot_df)[2], n_cells))
      plot_df2[, 4] <- c(rep(scRNA_name, nrow(plot_df2)))
      plot_df2[, 5] <- c(rep(top_n_pairs[i], nrow(plot_df2)))

      ## plot
      plot_df2 <- as.data.frame(plot_df2)
      plot_df2$expression <- as.numeric(plot_df2$expression)

      p <- ggplot(plot_df2, aes(x = cellid, y = expression, group = gene)) +
        geom_line(aes(linetype = gene, color = gene, size = gene)) +
        geom_point() +
        scale_linetype_manual(values = c("solid", "solid")) +
        scale_color_manual(values = c("#999999", "#E69F00")) +
        scale_size_manual(values = c(1, 1.5)) +
        theme_classic() +
        labs(title = paste0(scRNA_name))
      plot_list_subset[[scRNA_name]] <- p
    }
    p <- plot_list_subset[[1]]
    if (length(plot_list_subset) > 1) {
      for (p_tem in 2:length(plot_list_subset)) {
        p <- p | plot_list_subset[[p_tem]]
      }
    }
    pdf(file = paste0("Cell-subtype: ", celltype, " ", top_n_pairs[i], ".pdf"), width = 20, height = 6)
    print(p)
    dev.off()
    plot_list[[top_n_pairs[i]]] <- p
  }
}

gene_pair_exp_distribution=function(genes,cellsubtypes,seurat){
  plot_df=matrix(data = NA,nrow = 0,ncol = 3)
  for(celltype in cellsubtypes){
    scrna_tem=subset(seurat,subset = (Cell_subtype==celltype))
    cat("The cell number of ",celltype," is: ",ncol(scrna_tem),"\n")
    exp=as.data.frame(scrna_tem@assays$RNA@counts[genes,])
    mat_tem=matrix(data = NA,nrow = nrow(exp)*ncol(exp),ncol = 3)
    
    tem=c()
    for(a in 1:ncol(exp)){
      tem=c(tem,exp[,a])
    }
    
    mat_tem[,1]=as.character(tem)
    mat_tem[,2]=rep(rownames(exp),each=ncol(exp))
    mat_tem[,3]=rep(celltype,nrow(mat_tem))
    
    plot_df=rbind(plot_df,mat_tem)
  }
  plot_df=as.data.frame(plot_df)
  plot_df[,1]=as.numeric(plot_df[,1])
  plot_df[,2]=as.character(plot_df[,2])
  plot_df[,3]=as.factor(plot_df[,3])
  
  colnames(plot_df)=c("Expression","Gene","Celltype")
  
  p <- ggboxplot(plot_df, x = "Celltype", y = "Expression",
                 color = "Gene", palette = "jco")
  p2 = p + stat_compare_means(aes(group = Gene),paired = TRUE)
  
  pdf(file = paste0(genes[1],"-",genes[2]," expression level.pdf"))
  print(p2)
  dev.off()
}


## visualize the certain gene pair rank in different cell subtype
visualize_certain_gene_pair <- function(gene1, gene2, cellsubtypes, scRNA_seurat_list, ori_cells = 15) {
  ## gene1, gene2 : chacracter, represents the gene pair to plot
  ## cell subtypes : a vector contain all the cell subtype in scRNA_seurat_list
  ## scRNA_seurat_list : seurat list of the scRNA-seq data
  ## cells represenet the cell number to visualize the gene pair rank
  dir.create(paste0(gene1, "_", gene2))
  setwd(paste0(gene1, "_", gene2))

  set.seed(619)

  for (i in 1:length(scRNA_seurat_list)) {
    dir.create(names(scRNA_seurat_list)[i])
    setwd(names(scRNA_seurat_list)[i])
    exp <- scRNA_seurat_list[[i]]@assays$RNA@counts
    exp <- exp[c(gene1, gene2), ]

    for (j in cellsubtypes) {
      dir.create(j)
      setwd(j)

      exp_tem <- exp[, which(scRNA_seurat_list[[i]]@meta.data$Cell_subtype == j)]

      if (ori_cells > ncol(exp_tem)) {
        cells <- ncol(exp_tem)
      } else {
        cells <- ori_cells
      }

      exp_tem <- exp_tem[, sample(c(1:ncol(exp_tem)), replace = F, cells)]

      plot_df <- matrix(NA, nrow = cells * 2, ncol = 5)
      colnames(plot_df) <- c("cellid", "expression", "gene", "dataset", "pair")

      plot_df[, 1] <- c(rep(1:cells, 2))
      plot_df[, 2] <- c(exp_tem[1, ], exp_tem[2, ])
      plot_df[, 3] <- c(rep(gene1, cells), rep(gene2, cells))
      plot_df[, 4] <- c(rep(names(scRNA_seurat_list)[i], nrow(plot_df)))
      plot_df[, 5] <- c(rep(paste0(gene1, "-", gene2), nrow(plot_df)))

      plot_df <- as.data.frame(plot_df)
      plot_df$expression <- as.numeric(plot_df$expression)

      p <- ggplot(plot_df, aes(x = cellid, y = expression, group = gene)) +
        geom_line(aes(linetype = gene, color = gene, size = gene)) +
        geom_point() +
        scale_linetype_manual(values = c("solid", "solid")) +
        scale_color_manual(values = c("#999999", "#E69F00")) +
        scale_size_manual(values = c(1, 1.5)) +
        theme_classic() +
        labs(title = paste0(names(scRNA_seurat_list)[i]))

      pdf(file = paste0("Cell-subtype: ", j, " ", paste0(gene1, "-", gene2), ".pdf"), width = 4, height = 3)
      print(p)
      dev.off()
      setwd("../")
    }
    setwd("../")
  }
  gene_pair_exp_distribution(c(gene1, gene2), cellsubtypes, scRNA_seurat_list[[1]])
  setwd("../")
  return(NULL)
}

## Lollipop chart
Lollipop_chart <- function(iteration_result) {
  data_df <- shape_to_dataframe(iteration_result)
  data_df <- as.data.frame(data_df)
  data_df <- cbind("Cell_subtype" = rownames(data_df), data_df)

  major_celltype <- load_major_celltype_name()
  data_df$Cell_type <- major_celltype[match(data_df$Cell_subtype, major_celltype$cell_subtype), 2]

  ## plot
  p <- ggdotchart(data_df,
    x = "Cell_subtype", y = "Test_validate",
    color = met.brewer("Tara", n = 1), # Color by groups
    palette = c("#00AFBB", "#E7B800", "#FC4E07"), # Custom color palette
    sorting = "descending", # Sort value in descending order
    add = "segments", # Add segments from y = 0 to dots
    rotate = TRUE, # Rotate vertically
    # group = "Cell_type",                          # Order by groups
    dot.size = 12, # Large dot size
    label = round(data_df$Test_validate, digits = 4), # Add mpg values as dot labels
    font.label = list(
      color = "white", size = 10,
      vjust = 0.5
    ), # Adjust label parameters
    ggtheme = theme_pubr() # ggplot2 theme
  ) + ylab("Testing success itration number / Validation success itration number") + xlab("Cell subtypes")

  return(p)
}

## use best c-idnex in each model to visualize the Lollipop chart
Lollipop_chart_cindex <- function() {
  path <- "./model/"
  cell_types <- list.files(path)

  plot_df <- matrix(data = NA, nrow = 0, ncol = 2)
  colnames(plot_df) <- c("Cell_subtype", "C_index")

  for (i in cell_types) {
    iteratinos <- list.files(paste0(path, i, "/"))
    if (length(iteratinos) == 0) {
      next
    }
    max_cindex <- 0
    id <- -1
    for (j in iteratinos) {
      model <- readRDS(paste0(path, i, "/", j, "/", "Traning set(size=367) prognostic-signature of ", i, "(ncells=10)_model.rds"))
      c_index_tem <- model$cindex
      if (max_cindex < c_index_tem) {
        max_cindex <- c_index_tem
        id <- j
      }
    }
    plot_df <- rbind(plot_df, c(i, max_cindex))
  }

  plot_df <- as.data.frame(plot_df)
  plot_df[, 2] <- as.numeric(plot_df[, 2])

  p <- ggdotchart(plot_df,
    x = "Cell_subtype", y = "C_index",
    color = met.brewer("Tara", n = 1), # Color by groups
    palette = c("#00AFBB", "#E7B800", "#FC4E07"), # Custom color palette
    sorting = "descending", # Sort value in descending order
    add = "segments", # Add segments from y = 0 to dots
    rotate = TRUE, # Rotate vertically
    # group = "Cell_type",                          # Order by groups
    dot.size = 12, # Large dot size
    label = round(plot_df$C_index, digits = 2), # Add mpg values as dot labels
    font.label = list(
      color = "white", size = 10,
      vjust = 0.5
    ), # Adjust label parameters
    ggtheme = theme_pubr() # ggplot2 theme
  ) + ylab("C-index") + xlab("Cell subtype")

  pdf(file = "Lipport chart of C-index.pdf", width = 8, height = 12)
  print(p)
  dev.off()
}

## multiple variable cox analysis
MultiCox <- function(vali_clinical,iteration) {
   label1 <- 0
   label2 <- 0
   label3 <- 0

  ## mul-cox GSE14333
  if (1) {
    mul_cox <- coxph(Surv(dfs_time, dfs_event) ~
      Location + Duke_stage + Age + Gender + label_tem, data = vali_clinical[[1]])
    mul_cox1_result <- summary(mul_cox)
    colnames(mul_cox1_result$conf.int)

    ## forest
    multi1 <- as.data.frame(round(mul_cox1_result$conf.int[, c(1, 3, 4)], 2))
    multi2 <- ShowRegTable(mul_cox,
      exp = TRUE,
      digits = 2, pDigits = 3,
      printToggle = TRUE, quote = FALSE, ciFun = confint
    )
    GSE14333_result <- cbind(multi1, multi2)
    GSE14333_result <- cbind(Features = rownames(GSE14333_result), GSE14333_result)

    head(GSE14333_result)
    GSE14333_result[, 1] <- as.character(GSE14333_result[, 1])
    GSE14333_result[, 5] <- as.character(GSE14333_result[, 5])
    GSE14333_result[, 6] <- as.character(GSE14333_result[, 6])

    rm(multi1, multi2, mul_cox, mul_cox1_result)
    if (as.numeric(GSE14333_result[nrow(GSE14333_result), ncol(GSE14333_result)]) <= 0.05 | GSE14333_result[nrow(GSE14333_result), ncol(GSE14333_result)] == "<0.001") {
      label1 <- 1
    }
  }

  ## mul-cox GSE17536
  if (1) {
    print(colnames(vali_clinical[["GSE17536"]]))
    mul_cox <- coxph(Surv(dfs_time, dfs_event) ~
      Age + AJCC_stage + Gender + Histologic_grade + label_tem, data = vali_clinical[[2]])
    mul_cox1_result <- summary(mul_cox)
    colnames(mul_cox1_result$conf.int)

    ## forest
    multi1 <- as.data.frame(round(mul_cox1_result$conf.int[, c(1, 3, 4)], 2))
    multi2 <- ShowRegTable(mul_cox,
      exp = TRUE,
      digits = 2, pDigits = 3,
      printToggle = TRUE, quote = FALSE, ciFun = confint
    )
    GSE17536_result <- cbind(multi1, multi2)
    GSE17536_result <- cbind(Features = rownames(GSE17536_result), GSE17536_result)

    head(GSE17536_result)
    GSE17536_result[, 1] <- as.character(GSE17536_result[, 1])
    GSE17536_result[, 5] <- as.character(GSE17536_result[, 5])
    GSE17536_result[, 6] <- as.character(GSE17536_result[, 6])

    if (as.numeric(GSE17536_result[nrow(GSE17536_result), ncol(GSE17536_result)]) <= 0.05 | GSE17536_result[nrow(GSE17536_result), ncol(GSE17536_result)] == "<0.001") {
      label2 <- 1
    }

    rm(multi1, multi2, mul_cox, mul_cox1_result)
  }

  ## mul-cox GSE39582
  if (1) {
    print(colnames(vali_clinical[["GSE39582"]]))
    mul_cox <- coxph(Surv(dfs_time, dfs_event) ~
      Age + Gender + TNM_stage + BRAF_mutation + KRAS_mutation + TP53_mutation + MMR_status + label_tem, data = vali_clinical[[3]])
    mul_cox1_result <- summary(mul_cox)
    colnames(mul_cox1_result$conf.int)

    ## forest
    multi1 <- as.data.frame(round(mul_cox1_result$conf.int[, c(1, 3, 4)], 2))
    multi2 <- ShowRegTable(mul_cox,
      exp = TRUE,
      digits = 2, pDigits = 3,
      printToggle = TRUE, quote = FALSE, ciFun = confint
    )
    GSE39582_result <- cbind(multi1, multi2)
    GSE39582_result <- cbind(Features = rownames(GSE39582_result), GSE39582_result)

    head(GSE39582_result)
    GSE39582_result[, 1] <- as.character(GSE39582_result[, 1])
    GSE39582_result[, 5] <- as.character(GSE39582_result[, 5])
    GSE39582_result[, 6] <- as.character(GSE39582_result[, 6])

    if (as.numeric(GSE39582_result[nrow(GSE39582_result), ncol(GSE39582_result)]) <= 0.05 | GSE39582_result[nrow(GSE39582_result), ncol(GSE39582_result)] == "<0.001") {
      label3 <- 1
    }

    rm(multi1, multi2, mul_cox, mul_cox1_result)
  }

  ## Definate space
  ins <- function(x) {
    c(as.character(x), rep(NA, ncol(GSE14333_result) - 1))
  }

  ## result matrix
  if (1) {
    result_df <- rbind(
      c("Features", NA, NA, NA, "HR(95%CI)", "p-value"),
      ins("GSE14333"),
      ins("Location(Left/Rectum/Right)"),
      GSE14333_result[c(1, 2), ],
      ins("Bulk Stage"),
      GSE14333_result[3, ],
      ins("Age"),
      GSE14333_result[4, ],
      ins("Gender"),
      GSE14333_result[5, ],
      ins("Signature Label"),
      GSE14333_result[6, ],
      ins("GSE17536"),
      ins("AJCC Stage"),
      GSE17536_result[2, ],
      ins("Age"),
      GSE17536_result[1, ],
      ins("Gender"),
      GSE17536_result[3, ],
      ins("Histologic Grade"),
      GSE17536_result[4:5, ],
      ins("Signature Label"),
      GSE17536_result[6, ],
      ins("GSE39582"),
      ins("TNM Stage"),
      GSE39582_result[3, ],
      ins("Age"),
      GSE39582_result[1, ],
      ins("Gender"),
      GSE39582_result[2, ],
      ins("BRAF Mutation"),
      GSE39582_result[5, ],
      ins("KRAS Mutation"),
      GSE39582_result[7, ],
      ins("TP53 Mutation"),
      GSE39582_result[9, ],
      ins("MMR status"),
      GSE39582_result[11, ],
      ins("Signature Label"),
      GSE39582_result[12, ],
      c(NA, NA, NA, NA, NA, NA)
    )
  }

  ## hight-light rows
  is_summary_vector <- c()
  for (i in 1:nrow(result_df)) {
    if (is.na(result_df[i, 2])) {
      is_summary_vector <- c(is_summary_vector, TRUE)
    } else {
      is_summary_vector <- c(is_summary_vector, FALSE)
    }
  }
  rm(i)

  ## rename
  if (1) {
    rownames(result_df) <- 1:nrow(result_df)
    result_df[, 1] <- as.character(result_df[, 1])
    ## GSE14333
    result_df[4, 1] <- "Rectum versus Left"
    result_df[5, 1] <- "Right versus Left"
    result_df[7, 1] <- "StageIII versus StageII"
    result_df[9, 1] <- "Age>=70 versus Age<70"
    result_df[11, 1] <- "Male versus Female"
    result_df[13, 1] <- "Low-risk versus High-risk"
    ## GSE17536
    result_df[16, 1] <- "StageIII versus StageII"
    result_df[18, 1] <- "Age>=70 versus Age<70"
    result_df[20, 1] <- "Male versus Female"
    result_df[22, 1] <- "2-MD(Moderately differentiated) versus 1-WD(Well differentiated)"
    result_df[23, 1] <- "3-PD(Poorly differentiated) versus 1-WD(Well differentiated)"
    result_df[25, 1] <- "Low-risk versus High-risk"
    ## GSE39582
    result_df[28, 1] <- "StageIII versus StageII"
    result_df[30, 1] <- "Age>=70 versus Age<70"
    result_df[32, 1] <- "Male versus Female"
    result_df[34, 1] <- "BRAF Mutation versus Wild-type"
    result_df[36, 1] <- "KRAS Mutation versus Wild-type"
    result_df[38, 1] <- "TP53 versus Wild-type"
    result_df[40, 1] <- "pMMR versus Wild-dMMR"
    result_df[42, 1] <- "Low-risk versus High-risk"
  }
  ## forest plot
  p <- forestplot(result_df[, c(1, 5, 6)],
    mean = as.numeric(result_df[, 2]),
    lower = as.numeric(result_df[, 3]),
    upper = as.numeric(result_df[, 4]),
    zero = 1,
    boxsize = 0.6,
    graph.pos = "right",
    hrzl_lines = list(
      "1" = gpar(lty = 1, lwd = 2),
      "2" = gpar(lty = 2),
      "43" = gpar(lwd = 2, lty = 1, columns = c(1:4))
    ),
    graphwidth = unit(.25, "npc"),
    xlab = "HR(exp(coef))",
    xticks = c(0.4, 1, 3, 5, 7, 10),
    is.summary = is_summary_vector,
    txt_gp = fpTxtGp(
      label = gpar(cex = 1),
      ticks = gpar(cex = 1),
      xlab = gpar(cex = 1.5),
      title = gpar(cex = 2)
    ),
    lwd.zero = 1,
    lwd.ci = 1.5,
    lwd.xaxis = 2,
    lty.ci = 1.5,
    ci.vertices = T,
    ci.vertices.height = 0.2,
    clip = c(0.1, 8),
    ineheight = unit(8, "mm"),
    line.margin = unit(8, "mm"),
    colgap = unit(6, "mm"),
    fn.ci_norm = "fpDrawDiamondCI",
    title = "Multi-variables forest plot",
    col = fpColors(
      box = "#021eaa",
      lines = "#021eaa",
      zero = "black"
    )
  )

  pdf(
    file = paste0("Multi-variable forest of recurrence predictor of ", celltype, ".pdf"),
    width = 12, height = 12
  )
  print(p)
  dev.off()
  if (label1 & label2 & label3) {
    setwd("../")
    file.rename(paste0("iteration", iteration, "/"), paste0("iteration", iteration, "_1/"))
    setwd(paste0("iteration", iteration, "_1"))
  }
  return(NULL)
}
## plot cell fraction and markers in annotation Seurat object
plot_MarkersAndCellfraction <- function(dataset_list) {

    ## define colors
    colors <- c(
        brewer.pal(n = 99, name = "Set1"), brewer.pal(n = 99, name = "Set2"), brewer.pal(n = 99, name = "Set3"),
        brewer.pal(n = 99, name = "Pastel1"), brewer.pal(n = 99, name = "Pastel2"), brewer.pal(n = 99, name = "Paired"), brewer.pal(n = 99, name = "YlGn")
    )

    cellsubtypes <- names(table(dataset_list[[1]]@meta.data$Cell_subtype))
    for (i in 2:length(dataset_list)) {
        cellsubtypes <- union(cellsubtypes, dataset_list[[i]]@meta.data$Cell_subtype)
    }

    colors_df <- data.frame("colors" = colors[1:length(cellsubtypes)], "Cell_subtype" = cellsubtypes)

    ## annotation results
    for (i in datasets) {
        anno_seurat <- dataset_list[[i]]

        ## cell type makers and annotation
        Idents(anno_seurat) <- anno_seurat@meta.data$Cell_type

        ### cell type makers
        pdf(file = paste0("./figure/annotation/", i, " Major Celltype makers FeaturePlot.pdf"), width = 10, height = 7.5)
        p <- FeaturePlot(anno_seurat, features = c(
            "KIT", "CPA3", ## Mast cells
            "CD79A", "CD79B", ## B cells
            "CD3D", "CD3G", ## T cells
            "CD68", ## Myeloid cells
            "EPCAM", ## Epithelial cells
            "THY1", "COL3A1", ## Fibroblasts
            "PECAM1" ## Endothelial cells
        ))
        print(p)
        dev.off()

        ### cell type annotation
        pdf(file = paste0("./figure/annotation/", i, " Major Celltype Umap.pdf"), width = 15, height = 7.5)
        p1 <- DimPlot(anno_seurat, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 5, combine = T, group.by = "Cell_type")
        p2 <- DimPlot(anno_seurat, reduction = "umap", label = TRUE, pt.size = 0.5, label.size = 5, combine = T, group.by = "Patient")
        print(p1 + p2)
        dev.off()

        ### cell subtype annotation
        Idents(anno_seurat) <- anno_seurat@meta.data$Cell_subtype
        pdf(file = paste0("./figure/annotation/", i, " Cell subtype Umap.pdf"), width = 20, height = 7.5)
        p1 <- DimPlot(anno_seurat, reduction = "umap", label = F, pt.size = 0.5, label.size = 5, combine = T, group.by = "Cell_subtype")
        p2 <- DimPlot(anno_seurat, reduction = "umap", label = F, pt.size = 0.5, label.size = 5, combine = T, group.by = "Patient")
        print(p1 + p2)
        dev.off()

        ### cell type fraction
        p <- cell_subtype_percent_bar_ploting(anno_seurat, type = "major", colors_df)
        pdf(file = paste0("./figure/annotation/", i, " Major Celltype fraction.pdf"), width = 8, height = 6)
        print(p)
        dev.off()

        ### cell subtype fraction
        p <- cell_subtype_percent_bar_ploting(anno_seurat, type = "subtype", colors_df)
        pdf(file = paste0("./figure/annotation/", i, " Cell subtype fraction.pdf"), width = 12, height = 6)
        print(p)
        dev.off()
    }
    return(NULL)
}

## the cell subtype fractions in each tissue(barplot)
cell_subtype_percent_bar_ploting <- function(seurat_object, type, colors) {

    ### extract annotation information
    plot_df <- seurat_object@meta.data[, c("Class", "Cell_type", "Cell_subtype")]

    ### plot bar plot of major cell type
    if (type == "major") {
        p <- ggplot(data = plot_df, mapping = aes(x = Class, fill = Cell_type)) +
            geom_bar(stat = "count", width = 0.5, position = "fill") +
            guides(fill = guide_legend(reverse = T)) +
            scale_fill_manual(values = colors[1:length(table(plot_df$Cell_type))]) +
            labs(x = "Tissues", y = "Celltype Proportion", title = "") +
            theme_classic() +
            coord_flip()

        return(p)
    }

    ### bar plot of cell subtype
    if (type == "subtype") {
        p <- ggplot(data = plot_df, mapping = aes(x = Class, fill = Cell_subtype)) +
            geom_bar(stat = "count", width = 0.5, position = "fill") +
            guides(fill = guide_legend(reverse = T)) +
            scale_fill_manual(values = colors_df[match(names(table(plot_df$Cell_subtype)), colors_df[, 2]), 1]) +
            labs(x = "Tissues", y = "Cell Subtype Proportion", title = "") +
            theme_classic() +
            coord_flip()

        return(p)
    }
}
