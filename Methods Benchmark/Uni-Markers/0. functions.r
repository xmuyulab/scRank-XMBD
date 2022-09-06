## find makers of each cell-subtype in all scRNA-seq datasets
find_markers <- function(scRNA_seurat, logFC_cutoff = 1, padj_cutoff = 0.2) {
    celltype <- union(scRNA_seurat@meta.data$Cell_subtype, NULL)
    Idents(scRNA_seurat) <- scRNA_seurat@meta.data$Cell_subtype
    return_list <- list()

    for (i in celltype) {
        if (length(which(scRNA_seurat@meta.data$Cell_subtype == i)) < 5) {
            next
        }
        makers_df <- FindMarkers(scRNA_seurat, ident.1 = i, min.pct = 0.25)
        if (padj_cutoff == 0) {
            markers_df <- makers_df[which(makers_df$p_val <= 0.05), ]
        } else {
            markers_df <- makers_df[which(makers_df$p_val_adj < padj_cutoff), ]
        }
        markers_df <- markers_df[which(markers_df$avg_log2FC > logFC_cutoff), ]
        return_list[[i]] <- markers_df
    }
    return(return_list)
}

## find the overlap genes in scRNA and bulk RNA
find_overlap_scRNA_bulkRNA <- function(sc_markers_list, bulk_genes) {
    for (i in 1:length(sc_markers_list)) {
        sc_markers_list[[i]] <- sc_markers_list[[i]][rownames(sc_markers_list[[i]])
        %in% bulk_genes, ]
    }
    return(sc_markers_list)
}

## univariable cox for all markers
univariable_cox_for_marker <- function(sc_markers_list,
                                       exp, clinical, padj_cutoff) {
    cox_df <- matrix(data = NA, ncol = 4, nrow = 0)
    colnames(cox_df) <- c("Cell_subtype", "genes", "HR", "pvalue")

    for (i in 1:length(sc_markers_list)) {
        for (j in 1:nrow(sc_markers_list[[i]])) {
            gene_tem <- rownames(sc_markers_list[[i]])[j]

            clinical_df <- cbind(clinical, label = exp[gene_tem, ])

            univeriable_cox_models <- coxph(
                formula = Surv(dfs_time, dfs_event) ~ label,
                data = as.data.frame(clinical_df)
            )

            info <- summary(univeriable_cox_models)

            pvalue <- signif(as.matrix(info$coefficients)[, 5], 2)
            HR <- signif(as.matrix(info$coefficients)[, 2], 2)

            cox_df <- rbind(cox_df, c(names(sc_markers_list)[i], gene_tem, HR, pvalue))
        }
    }
    cox_df <- as.data.frame(cox_df)
    cox_df <- cbind(cox_df, "padj" = p.adjust(as.numeric(cox_df[, 4])))

    if (padj_cutoff == 0) {
        cox_df <- cox_df[which(cox_df$pvalue <= 0.05), ]
    } else {
        cox_df <- cox_df[which(cox_df$padj <= padj_cutoff), ]
    }
    return(cox_df)
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