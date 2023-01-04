## form reference matrix in SciBet annotation
SciBetRef_transfrom <- function(seuratObj) {
    ## training dataset transform
    exp_training <- seuratObj@assays$RNA@counts
    print(dim(exp_training))
    exp_training[1:5, 1:5] ## row as gene features, column as patients

    ## view the celltype(major)
    celltype_label <- seuratObj@meta.data$Cell_type
    print(table(celltype_label))

    ## traning and query exp matrix
    exp_training <- t(as.matrix(exp_training))

    ## GSE144735 reference
    GSE144735_ref <- as.data.frame(exp_training)
    GSE144735_ref <- cbind(GSE144735_ref, label = as.character(celltype_label))

    ## training
    etest_gene <- SelectGene(GSE144735_ref, k = 50)

    ## visulize markers

    pdf(file = "./figure/annotation/MajorCelltype ref genes.pdf", width = 12, height = 5)
    p <- Marker_heatmap(GSE144735_ref, etest_gene)
    print(p)
    dev.off()

    return(GSE144735_ref)
}

## load and predict query datasets of major celltype via scibet
load_and_predict_query <- function(ref, dataset) {
    ## load
    tem_seurat <- readRDS(paste0("./result/", dataset, "_SeuratObj.rds"))
    query_exp <- tem_seurat@assays$RNA@counts
    print(dim(query_exp))

    query_exp <- t(as.matrix(query_exp))

    ## predict
    prd_query <- SciBet(ref, query_exp)

    tem_seurat@meta.data$Cell_type <- prd_query
    saveRDS(tem_seurat, file = paste0("./result/", dataset, "_SeuratObj_anno.rds"))

    return(NULL)
}

predict_query_subtype <- function(GSE144735_seurat, query_list, celltype_label_list) {

    ## splite the celltype
    for (celltype in celltype_label_list) {
        if (celltype == "Epithelial cells") {
            ## normal
            if (1) {
                GSE144735_seurat <- subset(GSE144735_seurat, Class == "Normal")
                training_exp <- split_training_celltype(GSE144735_seurat, celltype)


                ## markers
                etest_gene <- SelectGene(training_exp, k = 50)

                ## visulize markers
                pdf(file = paste0("./figure/annotation/", celltype, "(Normal) ref genes.pdf"), width = 12, height = 5)
                p <- Marker_heatmap(training_exp, etest_gene)

                ## makers
                etest_gene <- SelectGene(training_exp, k = 50)

                ## visulize makers
                pdf(file = "./figure/annotation/", celltype, "(Normal) ref genes.pdf", width = 12, height = 5)
                p <- Marker_heatmap(GSE144735_ref, etest_gene)

                print(p)
                dev.off()

                for (i in 1:length(query_list)) {
                    if (!(celltype %in% query_list[[i]]@meta.data$Cell_type)) {
                        next
                    }
                    query_list_tem <- list()
                    query_list_tem[[i]] <- subset(query_list[[i]], Class == "Normal")
                    query_exp <- split_query_celltype(query_list_tem[[i]], celltype)

                    query_label <- SciBet(training_exp, query_exp)

                    query_list[[i]]@meta.data[rownames(query_list[[i]]@meta.data) %in% rownames(query_exp), "Cell_subtype"] <- query_label
                }
            }

            ## tumor border
            if (1) {
                for (i in 1:length(query_list)) {
                    if (!(celltype %in% query_list[[i]]@meta.data$Cell_type)) {
                        next
                    }
                    query_list_tem <- list()
                    query_list_tem[[i]] <- subset(query_list[[i]], (Class == "Tumor" | Class == "Border"))
                    query_exp <- split_query_celltype(query_list_tem[[i]], celltype)

                    ## extract exp matrix
                    query_exp <- t(query_exp)

                    ## ID transform
                    row_name <- rownames(query_exp)
                    row_name <- bitr(row_name, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
                    query_exp <- query_exp[rownames(query_exp) %in% row_name$SYMBOL, ]
                    rownames(query_exp) <- row_name$ENTREZID[match(rownames(query_exp), row_name$SYMBOL)]

                    ## scale
                    query_exp <- scale(query_exp)

                    query_label <- CMSclassifier::classifyCMS(query_exp, method = "SSP")[[3]]$SSP.nearestCMS

                    query_list[[i]]@meta.data[rownames(query_list[[i]]@meta.data) %in% colnames(query_exp), "Cell_subtype"] <- query_label
                }
            }
            next
        } else {
            training_exp <- split_training_celltype(GSE144735_seurat, celltype)

            ## markers0
            etest_gene <- SelectGene(training_exp, k = 50)

            ## visulize markers
            pdf(file = paste0("./figure/annotation/", celltype, " ref genes.pdf"), width = 12, height = 5)
            p <- Marker_heatmap(training_exp, etest_gene)

            ## makers0
            etest_gene <- SelectGene(training_exp, k = 50)

            ## visulize makers
            pdf(file = "./figure/annotation/", celltype, " ref genes.pdf", width = 12, height = 5)
            p <- Marker_heatmap(GSE144735_ref, etest_gene)

            print(p)
            dev.off()

            for (i in 1:length(query_list)) {
                if (!(celltype %in% query_list[[i]]@meta.data$Cell_type)) {
                    next
                }
                query_exp <- split_query_celltype(query_list[[i]], celltype)

                query_label <- SciBet(training_exp, query_exp)

                query_list[[i]]@meta.data$Cell_subtype[query_list[[i]]@meta.data$Cell_type == celltype] <- query_label
            }
        }
    }

    for (celltype in c("Mast cells")) {
        for (i in 1:length(query_list)) {
            query_list[[i]]@meta.data$Cell_subtype[query_list[[i]]@meta.data$Cell_type == celltype] <- query_list[[i]]@meta.data$Cell_type[query_list[[i]]@meta.data$Cell_type == celltype]
        }
    }


    saveRDS(query_list[[1]], file = paste0("./result/", strsplit(names(query_list), split = "_")[[1]][1], "_SeuratObj_anno.rds"))
    saveRDS(query_list[[2]], file = paste0("./result/", strsplit(names(query_list), split = "_")[[2]][1], "_SeuratObj_anno.rds"))

    saveRDS(query_list[[1]], file = "./result/", strsplit(names(query_list[[1]]), split = "_")[[1]][1], "_SeuratObj_anno.rds")
    saveRDS(query_list[[2]], file = "./result/", strsplit(names(query_list[[2]]), split = "_")[[1]][1], "_SeuratObj_anno.rds")

    return(NULL)
}

## split the major celltype and predict cell subtype by scibet
split_training_celltype <- function(seurat, celltype) {
    cat("The celltype of GSE144735 is ", celltype, "\n")
    seurat_tem <- subset(seurat, Cell_type == celltype)

    tem_exp <- seurat_tem@assays$RNA@counts
    tem_exp <- as.data.frame(t(as.matrix(tem_exp)))

    tem_exp <- cbind(tem_exp, label = seurat_tem@meta.data$Cell_subtype)
    return(tem_exp)
}

split_query_celltype <- function(seurat, celltype) {
    cat("The celltype of query dataset is ", celltype, "\n")
    seurat_tem <- subset(seurat, Cell_type == celltype)

    tem_exp <- seurat_tem@assays$RNA@counts
    tem_exp <- as.data.frame(t(as.matrix(tem_exp)))

    return(tem_exp)
}

## remove the cell-subtype which cell number less than cutoff n
remove_celltype_less_n <- function(exp_matrix, anno_matrix, n) {
    return_list <- list()

    tem <- as.data.frame(table(anno_matrix$Cell_subtype) < n)
    tem <- cbind(rownames(tem), tem)
    tem <- rownames(tem[which(tem[, 2] == TRUE), ])

    index_remove <- anno_matrix$Cell_subtype %in% tem
    anno_matrix_2 <- anno_matrix[!index_remove, ]

    exp_matrix <- exp_matrix[, !index_remove]

    ## remove na
    na_index <- !(is.na(anno_matrix_2$Cell_subtype))
    anno_matrix_2 <- anno_matrix_2[na_index, ]
    exp_matrix <- exp_matrix[, na_index]

    ## remove empyty
    empyty_index <- anno_matrix_2$Cell_subtype %in% ""

    anno_matrix_2 <- anno_matrix_2[!empyty_index, ]
    exp_matrix <- exp_matrix[, !empyty_index]

    rownames(anno_matrix_2) <- colnames(exp_matrix)

    return_list[[1]] <- exp_matrix
    return_list[[2]] <- anno_matrix_2

    return(return_list)
}

## get C-GPs
get_C_Gps <- function(average_exp_of_top_genepairs, scRNA_name, deltaS, ncells) {
    specific_genepairs_list <- list()
    for (i in 1:ncol(average_exp_of_top_genepairs)) {
        celltype_tem <- colnames(average_exp_of_top_genepairs)[i]
        specific_genepairs_list[[i]] <- select_celltype_specific_genepairs(
            average_exp_of_top_genepairs, celltype_tem,
            cutoff = deltaS
        )
    }
    names(specific_genepairs_list) <- sapply(
        colnames(average_exp_of_top_genepairs), function(x) {
            strsplit(x, split = "_Avg")[[1]][1]
        }
    )

    ## save the celltype specific gene pairs
    saveRDS(specific_genepairs_list, file = paste0(
        "./result/Cell-subtype-specific gene pairs(", scRNA_name,
        ",ncells=", ncells, ",deltaS=", deltaS, ").rds"
    ))



    ## remove the cell-subtype which celltype specific gene pairs less than n pairs
    specific_genepairs_list <- filter_cellsubtype_less_n(specific_genepairs_list, n = 10)

    return(specific_genepairs_list)
}

## select celltype-specific classcifier gene pairs
select_celltype_specific_genepairs <- function(gpTab, celltype, cutoff) {
    celltype_gene_pairs_exp <- gpTab[, celltype]
    if (length(celltype_gene_pairs_exp) == 0) {
        return(NULL)
    }

    cat(class(celltype_gene_pairs_exp), "\n")
    cat(length(celltype_gene_pairs_exp), "\n")

    ## calculate ΔS
    mean_exp_in_all_celltype <- apply(gpTab, MARGIN = 1, FUN = mean) ## mean of a gene-pair expression in all celltype
    cat("Length of Gene-pairs", length(mean_exp_in_all_celltype), "\n")

    ## ΔS
    delta_S <- celltype_gene_pairs_exp - mean_exp_in_all_celltype
    names(delta_S) <- rownames(gpTab)
    cat("ΔS of ", celltype, " :", "\n", delta_S[1:5], "\n")

    ## sort by ΔS
    delta_S <- sort(delta_S, decreasing = T)

    return(delta_S[which(delta_S > cutoff)])
}

## get stable C-GPs
get_stable_C_Gps <- function(scRNA_name_vector, scRNA_seurat_list, ncells, deltaS, stable_cutoff = 0.6) {
    ## load exp and specific_genepairs_list
    exp_list <- list()
    anno_list <- list()
    specific_gene_pairs_list <- list()
    for (scRNA_name in scRNA_name_vector) {
        exp_list[[scRNA_name]] <- scRNA_seurat_list[[scRNA_name]]@assays$RNA@counts
        anno_list[[scRNA_name]] <- scRNA_seurat_list[[scRNA_name]]@meta.data
        specific_gene_pairs_list[[scRNA_name]] <- readRDS(
        file = paste0(
        "./result/Cell-subtype-specific gene pairs(", scRNA_name,
        ",ncells=", ncells, ",deltaS=", deltaS, ").rds"))
    }

    ## overlap of specific gene pairs in all scRNA-seq datasets
    celltype_names <- c()
    for (i in 1:length(specific_gene_pairs_list)) {
        celltype_names <- union(celltype_names, names(specific_gene_pairs_list[[i]]))
    }

    ## calculate the consistance ratio of gene pairs from different scRNA-seq
    consistance_list <- list()
    for (celltype in celltype_names) {
        consistance_list[[celltype]] <- calculate_consistance_ratio(exp_list, anno_list, specific_gene_pairs_list, celltype)
    }

    ## filter the celltype-specific gene pairs consistance ratio higher than cutoff
    specific_gene_pairs_list <- stable_specific_gene_pairs(consistance_list, specific_gene_pairs_list, stable_cutoff = stable_cutoff)
    for (i in 1:length(consistance_list)) {
        consistance_list[[i]] <- consistance_list[[i]][which(consistance_list[[i]][, "mean"] > stable_cutoff), ]
    }
    saveRDS(specific_gene_pairs_list, file = paste0(
        "./result/Stable-cell-subtype-specific gene pairs(stable_cutoff=", stable_cutoff,
        ",ncells=", ncells, ",deltaS=", deltaS, ").rds"
    ))

    return(list(specific_gene_pairs_list, consistance_list, celltype_names))
}

## remove the cell-subtype which celltype specific gene pairs less than n pairs
filter_cellsubtype_less_n <- function(specific_genepairs_list, n) {
    remain_specific_genepairs_list <- list()
    remain_index <- c()
    j <- 1
    for (i in 1:length(specific_genepairs_list)) {
        if (length(specific_genepairs_list[[i]]) >= n) {
            remain_specific_genepairs_list[[j]] <- specific_genepairs_list[[i]]
            remain_index <- c(remain_index, i)
            j <- j + 1
        }
    }
    names(remain_specific_genepairs_list) <- names(specific_genepairs_list)[remain_index]
    return(remain_specific_genepairs_list)
}

filter_cellsubtype_less_n_v2 <- function(stable_pairs_list, n) {
    for (i in length(stable_pairs_list):1) {
        if (length(stable_pairs_list[[i]]) < n) {
            cat("The cell-subtype ", names(stable_pairs_list)[i], " only has ", length(stable_pairs_list[[i]]), " gene pairs, less than ", n, ".", "\n")
            stable_pairs_list[[i]] <- NULL
        }
    }
    return(stable_pairs_list)
}


## identifying the C-GPs with prognositic value
get_C_Gps_with_prognosis <- function(training_exp,training_clinical,specific_genepairs_list, scRNA_name, ncells, deltaS, clinical_cutoff) {
    prognostic_specific_genepairs_list <- list()
    if (clinical_cutoff == 0) {
        p_adj <- FALSE
    } else {
        p_adj <- TRUE
    }
    for (i in 1:length(specific_genepairs_list)) {
        prognostic_specific_genepairs_list[[i]] <- select_each_celltype_prognostic_specific_genepairs(
            celltype_specific_genepair = specific_genepairs_list[[i]],
            exp = training_exp,
            clinical = training_clinical,
            p_adj = p_adj,
            cutoff = clinical_cutoff
        )
    }
    names(prognostic_specific_genepairs_list) <- names(specific_genepairs_list)

    ## save the prognostic gene pairs
    saveRDS(prognostic_specific_genepairs_list, file = paste0(
        "./result/Cellsubtype prognostic gene pairs(", scRNA_name,
        ",ncells=", ncells, ",padj=", clinical_cutoff, ",deltaS=", deltaS, ").rds"
    ))

    return(prognostic_specific_genepairs_list)
}

split_genepairs=function(celltype_specific_genepair){
  genepairs=names(celltype_specific_genepair)
  gene_pair_df=matrix(data = NA,nrow = length(genepairs),ncol = 2)
  colnames(gene_pair_df)=c("Gene1","Gene2")
  gene1=sapply(genepairs,function(x){strsplit(x,split = '_')[[1]][1]})
  gene2=sapply(genepairs,function(x){strsplit(x,split = '_')[[1]][2]})
  gene_pair_df[,1]=gene1
  gene_pair_df[,2]=gene2
  
  return(gene_pair_df)
}

select_prognostic_genepair=function(exp,clinical,p_adj,cutoff){
  ## survival time equal to zero should be remove
  
  zero_time=clinical[,1]%in%0
  
  ## x matrix(features)
  exp=exp[,!zero_time]
  
  ## y matrix(survival time and status)
  y=as.matrix(Surv(time = as.double(clinical[,1]) ,event = as.double(clinical[,2])))
  y=y[!zero_time,]
  
  ## cox celltype-x specific gene pairs
  cox_for_each_features=matrix(data = NA,nrow = nrow(exp),ncol = 2)
  rownames(cox_for_each_features)=rownames(exp)
  colnames(cox_for_each_features)=c("pvalue","HR")
  for(i in 1:nrow(exp)){
    feature=rownames(exp)[i]
    data=cbind(y,'feature'=as.numeric(exp[i,]))
    univeriable_cox_formula=as.formula(Surv(y[,1],y[,2])~feature)
    
    univeriable_cox_models=coxph(formula = univeriable_cox_formula,data = as.data.frame(data))
    
    info=summary(univeriable_cox_models)
    
    pvalue=signif(as.matrix(info$coefficients)[,5],2)
    HR=signif(as.matrix(info$coefficients)[,2],2)
    
    cox_for_each_features[i,1]=pvalue
    cox_for_each_features[i,2]=HR
    
  }
  if(p_adj==TRUE){
    ## p value adjust
    cox_for_each_features = na.omit(cox_for_each_features)
    cox_for_each_features=as.data.frame(cox_for_each_features)
    cox_for_each_features$p_adj=p.adjust(cox_for_each_features$pvalue,method = 'BH')
    
    ## selcet p-adj minor than cutoff
    cox_for_each_features=cox_for_each_features[which(cox_for_each_features$p_adj<=cutoff),]
    
  }
  if(p_adj==FALSE){
    cox_for_each_features = na.omit(cox_for_each_features)
    cox_for_each_features=as.data.frame(cox_for_each_features)
    
    ## selcet p-value minor than cutoff
    cox_for_each_features=cox_for_each_features[which(cox_for_each_features$pvalue<=cutoff),]
    
  }
  ## features
  x=as.matrix(t(exp[rownames(exp)%in%rownames(cox_for_each_features),]))
  
  return_list=list()
  return_list[["genepairs"]]=x
  return_list[["genepair_prognostic_info"]]=cox_for_each_features
  return(return_list)
}

select_each_celltype_prognostic_specific_genepairs=function(celltype_specific_genepair,exp,clinical,p_adj,cutoff){
  ## splite gene-pairs into two genes
  celltype_specific_genepair=split_genepairs(celltype_specific_genepair)
  
  ## calculate the gene-pairs rank in bulk exp data
  bulk_exp_genepairs_rank=calculate_genepair_rank(exp,celltype_specific_genepair)
  
  ## removing some features which all equal 1 or 0
  bulk_exp_genepairs_rank=remove_all_1or0_features(bulk_exp_genepairs_rank)
  
  ## select prognostic gene pairs
  Celltype_prognostic_genepair=select_prognostic_genepair(bulk_exp_genepairs_rank,clinical,p_adj,cutoff)
  
  return(Celltype_prognostic_genepair)
}

## load major celltype names
load_major_celltype_name <- function() {
    B_cells <- c(
        "CD19+ CD20+ Plasma B cell", "Cycling B cell", "IgA+ IGLC2+ Plasma B cell",
        "IgA+ IGLL5+ Plasma B cell", "IgA+ Plasma B cell", "IgG+ Plasma B cell",
        "Naive B cell"
    ) ## 7

    T_cells <- c(
        "CD4+ IL7R+ T cells", "CD4+ Tfh", "CD4+ Th17", "CD8+ Exhausted-CTL",
        "CD8+ GZMK+ CTL", "CD8+ KLRD1+ CTL", "Gamma-delta T cells", "Naive T cells",
        "Regulatory T cells", "NK cells"
    ) ## 10

    Myeloid_cells <- c(
        "cDC_CD1C", "Macro_C1QC", "Macro_INHBA",
        "Macro_DNASE1L3", "Macro_LYVE1", "Macro_SPP1",
        "Mono_FCN1", "pDC_GZMB"
    ) ## 8

    Fibroblast_cells <- c(
        "eCAF", "Fibro_ADAMDEC1", "Fibro_BMP5", "Fibro_GPM6B",
        "Fibro_ID1", "Fibro_LY6H", "Fibro_OGN",
        "Fibro_SGK1", "iCAF", "mCAF", "myCAF_DES",
        "myCAF_NOTCH3", "Pericyte"
    ) ## 13

    Mast_cells <- c("Mast cells") ## 1

    Endothelial_cells <- c(
        "EC_ACKR1", "EC_ESM1", "EC_GPIHBP1", "EC_IGFBP3",
        "EC_PLVAP_CD36", "EC_PLVAP_PRCP",
        "EC_TFF3", "EC_Unspecific"
    ) ## 8

    Epithelial_cells <- c("CMS1", "CMS2", "CMS3", "CMS4") ## 4

    major_celltype_df <- matrix(nrow = 51, ncol = 2)
    major_celltype_df <- as.data.frame(major_celltype_df)
    major_celltype_df[, 1] <- c(
        B_cells, T_cells, Mast_cells, Fibroblast_cells, Endothelial_cells,
        Epithelial_cells, Myeloid_cells
    )
    major_celltype_df[, 2] <- c(
        rep("B cells", length(B_cells)),
        rep("T cells", length(T_cells)),
        rep("Mast cells", length(Mast_cells)),
        rep("Fibroblasts", length(Fibroblast_cells)),
        rep("Endothelial cells", length(Endothelial_cells)),
        rep("Epithelial cells", length(Epithelial_cells)),
        rep("Myeloid cells", length(Myeloid_cells))
    )
    colnames(major_celltype_df) <- c("cell_subtype", "celltype")
    return(major_celltype_df)
}

## calculate the consistance ratio of gene pairs from different scRNA-seq
calculate_consistance_ratio <- function(exp_list, anno_list, specific_gene_pairs_list, celltype) {
    ## remove the datasets which do not have the cell-subtype specific gene pairs
    for (i in length(specific_gene_pairs_list):1) {
        if (!celltype %in% names(specific_gene_pairs_list[[i]])) {
            cat(celltype, " do not be contained in dataset ", names(specific_gene_pairs_list)[i], ".", "\n")
            specific_gene_pairs_list[[i]] <- NULL
            exp_list[[i]] <- NULL
            anno_list[[i]] <- NULL
        }
    }

    onelabel <- F
    if (length(specific_gene_pairs_list) < 2) {
        cat("There is only one dataset contain ", celltype, ".", "\n")
        onelabel <- T
    }

    ## extrace celltype specific exp matrix
    celltype_exp_list <- list()
    for (z in 1:length(exp_list)) {
        exp_tem <- exp_list[[z]]
        anno_tem <- anno_list[[z]]

        index_tem <- which(anno_tem$Cell_subtype == celltype)
        celltype_exp_list[[z]] <- exp_tem[, index_tem]
    }
    names(celltype_exp_list) <- names(exp_list)
    rm(exp_tem, anno_tem, index_tem, z, exp_list, anno_list)

    ## extract cell-subtype specific gene pairs in each datasets
    celltype_specific_genepairs_vector <- c()
    for (i in 1:length(specific_gene_pairs_list)) {
        celltype_specific_genepairs_vector <- c(
            celltype_specific_genepairs_vector,
            specific_gene_pairs_list[[i]][[celltype]]
        )
    }
    rm(specific_gene_pairs_list)

    ## split gene-pairs name
    gene_pairs <- names(celltype_specific_genepairs_vector)
    genePair_df <- matrix(NA, nrow = length(celltype_specific_genepairs_vector), ncol = 2)
    colnames(genePair_df) <- c("gene1", "gene2")
    for (i in 1:nrow(genePair_df)) {
        genePair_df[i, 1] <- strsplit(gene_pairs[i], split = "_")[[1]][1]
        genePair_df[i, 2] <- strsplit(gene_pairs[i], split = "_")[[1]][2]
    }

    ## omit the gene pairs do not be contained in datasets exp
    remain_index <- rep(T, nrow(genePair_df))
    for (j in 1:length(celltype_exp_list)) {
        remain_index <- remain_index & genePair_df[, 1] %in% rownames(celltype_exp_list[[j]])
        remain_index <- remain_index & genePair_df[, 2] %in% rownames(celltype_exp_list[[j]])
    }
    genePair_df <- genePair_df[remain_index, ]
    gene_pairs <- gene_pairs[remain_index]
    rm(j, remain_index)

    ## gene pairs rank list
    genepair_rank_list <- list()
    for (j in 1:length(celltype_exp_list)) {
        exp_tem <- celltype_exp_list[[j]]
        exp_tem <- as.matrix(exp_tem)

        gene_pair_rank_df_tem <- exp_tem[genePair_df[, 1], ] - exp_tem[genePair_df[, 2], ]
        gene_pair_rank_df_tem <- ifelse(gene_pair_rank_df_tem > 0, 1, 0)

        rownames(gene_pair_rank_df_tem) <- gene_pairs
        colnames(gene_pair_rank_df_tem) <- colnames(exp_tem)
        genepair_rank_list[[j]] <- gene_pair_rank_df_tem
    }
    names(genepair_rank_list) <- names(celltype_exp_list)
    rm(gene_pair_rank_df_tem, exp_tem, j)

    ## calculate consistance ration
    rate_df <- matrix(NA, nrow = nrow(genepair_rank_list[[1]]), ncol = (length(genepair_rank_list) + 1))
    colnames(rate_df) <- c(names(genepair_rank_list), "mean")
    rownames(rate_df) <- gene_pairs

    for (i in 1:length(genepair_rank_list)) {
        rank_df_tem <- genepair_rank_list[[i]]
        rate_df[, i] <- apply(rank_df_tem, 1, "sum")
        rate_df[, i] <- rate_df[, i] / ncol(rank_df_tem)
    }

    if (onelabel == T) {
        rate_df[, "mean"] <- rate_df[, 1]
    }
    if (onelabel == F) {
        rate_df[, "mean"] <- apply(rate_df[, 1:length(genepair_rank_list)], 1, "mean")
    }
    rate_df <- rate_df[order(rate_df[, ncol(rate_df)], decreasing = T), ]
    rm(rank_df_tem)

    return(rate_df)
}

## filter the celltype-specific gene pairs consistance ratio higher than cutoff
stable_specific_gene_pairs <- function(consistance_list, specific_gene_pairs_list, stable_cutoff = 0.75) {
    for (i in 1:length(consistance_list)) {
        celltype <- names(consistance_list)[i]
        tem_df <- consistance_list[[celltype]]
        tem_df <- tem_df[which(tem_df[, "mean"] > stable_cutoff), ]
        genepairs <- rownames(tem_df)
        rm(tem_df)
        for (j in 1:length(specific_gene_pairs_list)) {
            specific_gene_pairs_list[[j]][[celltype]] <- specific_gene_pairs_list[[j]][[celltype]][names(specific_gene_pairs_list[[j]][[celltype]]) %in% genepairs]
        }
        rm(genepairs, j)
    }
    return(specific_gene_pairs_list)
}

## combine stable C-GPs from different scRNA-seq dataset
combine_StableCGPs <- function(scRNA_name_vector, ncells, deltaS, numPairs_cutoff = 50) {
    stable_celltype_specific_genepairs <- list()
    for (scRNA in scRNA_name_vector) {
        stable_celltype_specific_genepairs[[scRNA]] <- readRDS(paste0(
            "./result/Cell-subtype-specific gene pairs(", scRNA, ",ncells=", ncells, ",deltaS=", deltaS, ").rds"
        ))
    }

    celltype_list <- c()
    for (each_ in stable_celltype_specific_genepairs) {
        celltype_list <- union(celltype_list, names(each_))
    }

    stable_pairs_list <- list()
    for (celltype in celltype_list) {
        for (i in 1:length(stable_celltype_specific_genepairs)) {
            if (celltype %in% names(stable_celltype_specific_genepairs[[i]])) {
                stable_pairs_list[[celltype]] <- union(stable_pairs_list[[celltype]], names(stable_celltype_specific_genepairs[[i]][[celltype]]))
            }
        }
    }

    stable_pairs_list <- filter_cellsubtype_less_n_v2(stable_pairs_list, numPairs_cutoff)
    return((stable_pairs_list))
}

## K-M plot
KM_plot <- function(clinical, label) {
    clinical[, 1] <- as.numeric(clinical[, 1])
    clinical[, 2] <- as.numeric(clinical[, 2])

    data <- cbind(clinical, score = label)

    fit <- surv_fit(Surv(clinical[, 1], clinical[, 2]) ~ label, data = data)
    logrank_pvalue <- survdiff(Surv(clinical[, 1], clinical[, 2]) ~ label, data = data)
    logrank_pvalue <- 1 - pchisq(logrank_pvalue$chisq, length(logrank_pvalue$n) - 1)

    fit2 <- coxph(Surv(clinical[, 1], clinical[, 2]) ~ label, data = data)
    fit2_sum <- summary(fit2)
    c_index <- fit2_sum$concordance[1]
    c_index <- round(digits = 3, c_index)

    p <- ggsurvplot(fit,
        data = data,
        linetype = c("solid", "solid"),
        surv.median.line = "hv", surv.scale = "percent",
        pval = T, risk.table = T,
        conf.int = T, conf.int.alpha = 0.1, conf.int.style = "ribbon",
        risk.table.y.text = T,
        palette = c("#CC3300", "#3300CC"),
        xlab = "Disease-Free survival time(month)"
    )

    p$plot <- p$plot + labs(caption = paste0("C-idnex: ", c_index))

    return_list <- list()
    return_list[["p"]] <- p
    return_list[["logrank_p"]] <- logrank_pvalue
    return_list[["c-index"]] <- c_index

    return(return_list)
}

## analysis the KMs result
analysis_KMs <- function(ncells, celltype_list, path) {
    interation_result <- list()

    for (celltype in celltype_list) {
        setwd(path)
        setwd(paste0("./", celltype, "/"))

        iteration_files_path <- getwd()

        interation_files <- list.files(iteration_files_path)
        all_interations_num <- length(interation_files)
        training_success <- 0
        testing_success <- 0
        all_success <- 0

        for (i in interation_files) {
            files_num_tem <- length(list.files(paste0(iteration_files_path, "/", i, "/")))
            if (files_num_tem < 5) {
                unlink(i, recursive = TRUE)
            }

            if ((files_num_tem > 5 & files_num_tem < 6) | files_num_tem == 5) {
                training_success <- training_success + 1
                unlink(i, recursive = TRUE)
            }
            if (files_num_tem == 6 | (files_num_tem > 6 & files_num_tem < 9)) {
                training_success <- training_success + 1
                testing_success <- testing_success + 1
                unlink(i, recursive = TRUE)
            }
            if (files_num_tem >= 9) {
                training_success <- training_success + 1
                testing_success <- testing_success + 1
                all_success <- all_success + 1
            }
        }

        interation_result[[celltype]]$all_interations_num <- all_interations_num
        interation_result[[celltype]]$training_success <- training_success
        interation_result[[celltype]]$testing_success <- testing_success
        interation_result[[celltype]]$all_success <- all_success
    
    setwd("../../")
    }
    return(interation_result)
}

## write the signature information into a table
save_signature <- function(celltype, iteration) {
    RDS <- readRDS(paste0(
        "./model/",
        celltype, "/iteration", iteration, "_1/Traning set(size=367) prognostic-signature of ", celltype, "(ncells=10)_model.rds"
    ))
    model <- RDS[["lasso-cox model"]]

    beta <- coef(model)
    features <- rownames(beta)[which(beta[, 1] != 0)]

    final_table <- matrix(data = NA, nrow = length(features), ncol = 5)
    colnames(final_table) <- c("Signature", "REOs (Ra>Rb)", "Coefficient", "Lambda", "C-Index")

    final_table[, 1] <- paste0("Gene pair ", 1:length(features))

    final_table[, 2] <- sapply(features, function(x) {
        gene1 <- strsplit(x, split = "__")[[1]][1]
        gene2 <- strsplit(x, split = "__")[[1]][2]

        return(paste0(gene1, ">", gene2))
    })

    final_table[, 3] <- round(beta[which(beta[, 1] != 0), ], digits = 7)

    final_table[, 4] <- rep(model$lambda.1se, length(features))

    final_table[, 5] <- rep(RDS[["cindex"]], length(features))

    write.table(final_table, file = paste0("./result/",celltype, " signatures.txt"), sep = "\t", col.names = T, row.names = F, fileEncoding = "UTF-8")

    
    return(NULL)

}


## reshape to dataframe
shape_to_dataframe=function(iteration_result){
  celltype_names=names(iteration_result)
  
  ### form data frame
  result_df=matrix(data = NA,ncol = length(iteration_result[[1]]),nrow = length(iteration_result))
  colnames(result_df)=names(iteration_result[[1]])
  rownames(result_df)=celltype_names
  
  for(i in 1:length(celltype_names)){
    result_df[i,]=c(iteration_result[[celltype_names[i]]]$all_interations_num,
                    iteration_result[[celltype_names[i]]]$training_success,
                    iteration_result[[celltype_names[i]]]$testing_success,
                    iteration_result[[celltype_names[i]]]$all_success)
  }
  
  result_df=cbind(result_df,"Train_Test"=result_df[,3]/result_df[,2])
  result_df=cbind(result_df,"Test_validate"=result_df[,4]/result_df[,3])
  
  return(result_df)
}
