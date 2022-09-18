# functions for dataPreprocess


## pre-process of seurat pipeline
pre_processing_scRNAseq <- function(exp, clinical) {

    ## creat seurat object
    seurat_object <- CreateSeuratObject(counts = exp, project = "seurat_object", min.cells = 3)

    ## add meta data
    seurat_object <- AddMetaData(object = seurat_object, metadata = clinical$Class, col.name = "Class")
    seurat_object <- AddMetaData(object = seurat_object, metadata = clinical$Patient, col.name = "Patient")

    ## qc
    seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT")
    seurat_object <- subset(seurat_object, subset = nFeature_RNA < 6000 & nFeature_RNA > 200 & percent.mt < 20)

    ## normalize
    ## since import log(TPM+1) matrix, did not need to normalize
    ## seurat_object = NormalizeData(seurat_object,normalization.method = 'LogNormalize', scale.factor = 10000)

    ## find highly variable genes
    seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)

    ## sacle
    seurat_object <- ScaleData(object = seurat_object)

    ## PCA
    seurat_object <- RunPCA(object = seurat_object) ## 对于差异表达基因进行PCA降维

    return(seurat_object)
}

## apply CMSclassifier to annotate cell-subpopulation
CMSclassifier_annoFn <- function(seuratObj) {
    ## exp matrix
    exp <- seuratObj@assays$RNA@counts
    ## ID transfer
    row_name <- rownames(exp)
    row_name <- bitr(row_name, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
    exp <- exp[rownames(exp) %in% row_name$SYMBOL, ]
    rownames(exp) <- row_name$ENTREZID[match(rownames(exp), row_name$SYMBOL)]

    ## normalizing and scaling
    exp <- scale(exp)

    ## predicting CMS label
    res <- CMSclassifier::classifyCMS(exp, method = "SSP")[[3]]
    cms_label <- res$SSP.nearestCMS

    ## gene sets from CMScaller::geneSets.CMS
    gene_sets <- list()
    for (i in 1:length(geneSets.CMS)) {
        gene_sets[[i]] <- bitr(geneSets.CMS[[i]], fromType = "ENTREZID", toType = c("SYMBOL"), OrgDb = org.Hs.eg.db)
        gene_sets[[i]] <- gene_sets[[i]]$SYMBOL
    }
    names(gene_sets) <- names(geneSets.CMS)

    ## GSVA to identify the pathway
    exp <- seuratObj@assays$RNA@counts ## all genes
    exp <- scale(exp)

    GSVA_df <- matrix(NA, nrow = nrow(exp), ncol = length(table(cms_label)))
    for (i in 1:length(table(cms_label))) {
        label_tem <- names(table(cms_label))[i]
        exp_tem <- exp[, which(cms_label == label_tem)]
        GSVA_df[, i] <- apply(exp_tem, MARGIN = 1, FUN = "mean")
    }
    rownames(GSVA_df) <- rownames(exp)
    colnames(GSVA_df) <- sapply(c(1:4), function(x) {
        paste0("CMS", x)
    })

    GSVA_result_df <- gsva(as.matrix(GSVA_df), gene_sets, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE)

    p <- pheatmap(GSVA_result_df,
        scale = "row", show_colnames = T, show_rownames = T,
        cluster_cols = T, cluster_rows = F,
        clustering_method = "complete"
    )
    pdf(file = "./figure/annotation/GSE144735 Epithelial(Tumor) CMS pathway GSVA.pdf", width = 10, height = 7.5)
    print(p)
    dev.off()

    return(cms_label)
}

## this function uses the normal Seurat pipeline to annotate scRNA-seq dataset GSE144735
GSE144735_annoFn <- function() {
    ## major cell type annotation
    ## load data

    exp <- read.table(file = "./data/GSE144735_exp.txt", header = T, row.names = 1, sep = "\t")
    clinical <- read.table(file = "./data/GSE144735_anno.txt", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)

    ## pre-process of seurat pipeline
    seuratobj <- pre_processing_scRNAseq(exp, clinical)
    seuratobj <- FindNeighbors(seuratobj, dims = 1:11)
    seuratobj <- FindClusters(seuratobj, resolution = 0.9)
    seuratobj <- RunUMAP(seuratobj, dims = 1:11)

    ## plot makers
    pdf("./figure/annotation/GSE144735 Celltype makers FeaturePlot.pdf", width = 10, height = 7.5)
    p <- FeaturePlot(seuratobj, features = c(
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

    ## add cluster id into metadata
    seuratobj@meta.data$cell_typeid <- as.character(seuratobj@active.ident)
    seuratobj@meta.data$Cell_subtype <- NA

    ## rename clusters
    new.cluster.ids <- c(
        "T cells", "B cells", "T cells", "B cells", "Myeloid cells", "Fibroblasts",
        "Endothelial cells", "Epithelial cells", "Epithelial cells", "Fibroblasts", "Epithelial cells", "Epithelial cells",
        "Fibroblasts", "T cells", "Myeloid cells", "Fibroblasts", "Epithelial cells", "Epithelial cells",
        "Fibroblasts", "Fibroblasts", "B cells", "Fibroblasts", "Mast cells"
    )
    names(new.cluster.ids) <- levels(seuratobj)
    seuratobj <- RenameIdents(seuratobj, new.cluster.ids)

    ## add cluster label into metadata
    Cell_type <- as.character(seuratobj@active.ident)
    seuratobj@meta.data$Cell_type <- Cell_type

    p <- DimPlot(seuratobj, reduction = "umap", label = T, repel = F, pt.size = 0.5, label.size = 5, combine = T) +
        theme(plot.title = element_text(hjust = 0.5))
    ggsave(p, file = "./figure/annotation/GSE144735 Celltype annotation.pdf", width = 12, height = 10)

    saveRDS(seuratobj, file = "./result/GSE144735_SeuratObj.rds")

    ## the follow part is annotating the cell-subpopulation, the pipelien is similar to above
    ## B cell
    if (T) {
        B_seurat <- subset(seuratobj, subset = Cell_type == "B cells")
        B_seurat <- RunPCA(object = B_seurat)
        B_seurat <- FindNeighbors(B_seurat, dims = 1:10)
        B_seurat <- FindClusters(B_seurat, resolution = 1.1)
        B_seurat <- RunUMAP(B_seurat, dims = 1:10)

        ## rename clusters
        new.cluster.ids <- c(
            "CD19+ CD20+ Plasma B cell", "IgA+ Plasma B cell", "IgA+ IGLC2+ Plasma B cell", "IgA+ Plasma B cell", "Naive B cell", "IgA+ Plasma B cell",
            "IgG+ Plasma B cell", "Cycling B cell", "CD19+ CD20+ Plasma B cell", "IgA+ IGLC2+ Plasma B cell", "IgG+ Plasma B cell", "IgA+ Plasma B cell",
            "IgA+ IGLL5+ Plasma B cell", "CD19+ CD20+ Plasma B cell"
        )

        names(new.cluster.ids) <- levels(B_seurat)
        B_seurat <- RenameIdents(B_seurat, new.cluster.ids)
        Cell_subtype <- as.character(B_seurat@active.ident)
        B_seurat@meta.data$Cell_subtype <- Cell_subtype

        ## save the cell subtype information into whole seurat object
        seuratobj@meta.data$Cell_subtype[match(rownames(B_seurat@meta.data), rownames(seuratobj@meta.data))] <- Cell_subtype

        ## umap plot
        pdf(file = "./figure/annotation/GSE144735 Bcells annotation.pdf", width = 12, height = 7.5)
        p1 <- DimPlot(B_seurat, reduction = "umap", label = F, pt.size = 0.5, label.size = 5, combine = T) +
            theme(plot.title = element_text(hjust = 0.5))
        p2 <- DimPlot(B_seurat, reduction = "umap", label = F, pt.size = 0.5, label.size = 5, combine = T, split.by = "Class") +
            theme(plot.title = element_text(hjust = 0.5))
        p3 <- DimPlot(B_seurat, reduction = "umap", label = F, pt.size = 0.5, label.size = 5, combine = T, group.by = "Patient") +
            theme(plot.title = element_text(hjust = 0.5))
        print((p1 + p3) / p2)
        dev.off()

        ### makers dotplot
        if (1) {
            makers_list <- list()
            makers_list[["CD19+"]] <- c("MS4A1", "CD19")
            makers_list[["IgA+"]] <- c("IGHA1", "IGLC2", "IGLL5")
            makers_list[["Naive"]] <- c("IGHM", "IGHD")
            makers_list[["IgG+"]] <- c("IGHG1")
            makers_list[["Cycling"]] <- c("HMGN2")
        }

        pdf(file = "./figure/annotation/GSE144735 Bcells maker dotplot.pdf", width = 10, height = 7.5)
        p <- DotPlot(B_seurat, features = makers_list) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        print(p)
        dev.off()

        ## save seurat object
        # saveRDS(B_seurat, file = "./result/GSE144735_Bcell_SeuratObj.rds")
    }

    ## T cells
    if (T) {
        ### subset T cells
        T_seurat <- subset(seuratobj, subset = Cell_type == "T cells")
        T_seurat <- RunPCA(object = T_seurat)
        T_seurat <- FindNeighbors(T_seurat, dims = 1:10)
        T_seurat <- FindClusters(T_seurat, resolution = 1.1)
        T_seurat <- RunUMAP(T_seurat, dims = 1:10)

        ## rename clusters
        new.cluster.ids <- c(
            "Naive T cells", "CD8+ GZMK+ CTL", "CD4+ IL7R+ T cells", "Regulatory T cells", "Naive T cells", "CD8+ KLRD1+ CTL",
            "CD4+ Th17", "NK cells", "CD4+ Tfh", "CD8+ Exhausted-CTL", "Gamma-delta T cells"
        )
        names(new.cluster.ids) <- levels(T_seurat)
        T_seurat <- RenameIdents(T_seurat, new.cluster.ids)
        Cell_subtype <- as.character(T_seurat@active.ident)
        T_seurat@meta.data$Cell_subtype <- Cell_subtype

        ## save the cell subtype information into whole seurat object
        seuratobj@meta.data$Cell_subtype[match(rownames(T_seurat@meta.data), rownames(seuratobj@meta.data))] <- Cell_subtype

        ### umap plot
        pdf("./figure/annotation/GSE144735 Tcells annotation.pdf", width = 12, height = 7.5)
        p1 <- DimPlot(T_seurat, reduction = "umap", label = F, pt.size = 0.5, label.size = 5, combine = T) +
            theme(plot.title = element_text(hjust = 0.5))
        p2 <- DimPlot(T_seurat, reduction = "umap", label = F, pt.size = 0.5, label.size = 5, combine = T, split.by = "Class") +
            theme(plot.title = element_text(hjust = 0.5))
        p3 <- DimPlot(T_seurat, reduction = "umap", label = F, pt.size = 0.5, label.size = 5, combine = T, group.by = "Patient") +
            theme(plot.title = element_text(hjust = 0.5))
        print((p1 + p3) / p2)
        dev.off()

        ### makers dotplot
        if (1) {
            makers_list <- list()
            makers_list[["CD4+"]] <- c("CD4", "ICOS")
            makers_list[["CD8+"]] <- c("CD8A", "CD8B")
            makers_list[["Gamma-Delta"]] <- c("TRGC1", "TRGC2", "TRDC")
            makers_list[["Exhausted"]] <- c("LAG3", "TIGIT", "HAVCR2", "CTLA4")
            makers_list[["IL7R"]] <- c("IL7R")
            makers_list[["Tfh"]] <- c("MAF", "CXCL13", "PDCD1")
            makers_list[["NK"]] <- c("XCL1", "FCGR3A", "KLRF1")
            makers_list[["CTL"]] <- c("GZMA", "GNLY", "PRF1", "GZMB", "GZMK", "NKG7", "KLRD1", "KLRG1", "ITGAE")
            makers_list[["Regulatory T"]] <- c("IL2RA", "FOXP3", "TGFB1", "TGFB3", "TGFBI", "TGFBR1")
            makers_list[["Th17"]] <- c("IL17A", "IL17F", "IL22", "IL21")
            makers_list[["Naive"]] <- c("TCF7", "SELL", "LEF1", "CCR7")
        }
        pdf(file = "./figure/annotation/GSE144735 Tcells maker dotplot.pdf", width = 20, height = 7.5)
        p <- DotPlot(T_seurat, features = makers_list) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        print(p)
        dev.off()

        ## save seurat object
        # saveRDS(T_seurat, file = "GSE144735(All_tissue)_T cell_seurat_obj.rds")
    }

    ## Myeloid cells
    if (T) {
        Myeloid_seurat <- subset(seuratobj, subset = Cell_type == "Myeloid cells")

        Myeloid_seurat <- RunPCA(object = Myeloid_seurat)

        Myeloid_seurat <- FindNeighbors(Myeloid_seurat, dims = 1:10)
        Myeloid_seurat <- FindClusters(Myeloid_seurat, resolution = 1.1)
        Myeloid_seurat <- RunUMAP(Myeloid_seurat, dims = 1:10)

        ## rename clusters
        new.cluster.ids <- c(
            "Macro_C1QC", "Macro_C1QC", "Macro_SPP1", "Macro_LYVE1", "Macro_DNASE1L3", "Mono_FCN1",
            "cDC_CD1C", "Mono_FCN1", "Macro_SPP1", "Macro_SPP1", "Macro_LYVE1", "pDC_GZMB", "Macro_INHBA"
        )
        names(new.cluster.ids) <- levels(Myeloid_seurat)
        Myeloid_seurat <- RenameIdents(Myeloid_seurat, new.cluster.ids)
        Cell_subtype <- as.character(Myeloid_seurat@active.ident)
        Myeloid_seurat@meta.data$Cell_subtype <- Cell_subtype

        ### save the cell subtype information into whole seurat object
        seuratobj@meta.data$Cell_subtype[match(rownames(Myeloid_seurat@meta.data), rownames(seuratobj@meta.data))] <- Cell_subtype

        ### umap plot
        pdf(file = "./figure/annotation/GSE144735 Myeloid annotation.pdf", width = 12, height = 7.5)
        p1 <- DimPlot(Myeloid_seurat, reduction = "umap", label = F, pt.size = 0.5, label.size = 5, combine = T) +
            theme(plot.title = element_text(hjust = 0.5))
        p2 <- DimPlot(Myeloid_seurat, reduction = "umap", label = F, pt.size = 0.5, label.size = 5, combine = T, split.by = "Class") +
            theme(plot.title = element_text(hjust = 0.5))
        p3 <- DimPlot(Myeloid_seurat, reduction = "umap", label = F, pt.size = 0.5, label.size = 5, combine = T, group.by = "Patient") +
            theme(plot.title = element_text(hjust = 0.5))
        print((p1 + p3) / p2)
        dev.off()

        ## makers
        pdf(file = "./figure/annotation/GSE144735 Myeloid maker dotplot.pdf", width = 10, height = 7.5)
        if (1) {
            makers_list <- list()
            makers_list[["Macrophages"]] <- c("LYVE1", "PLTP", "SPP1", "APOE", "C1QC", "C1QA", "DNASE1L3", "INHBA", "IL1RN", "CCL4")
            makers_list[["Monocytes"]] <- c("FCN1", "S100A9", "S100A8")
            makers_list[["DCs"]] <- c("CD1C", "FCER1A", "HLA-DQA1", "GZMB", "LILRA4", "IL3RA")
        }
        p <- DotPlot(Myeloid_seurat, features = makers_list) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        print(p)
        dev.off()

        ## save seurat object
        # saveRDS(Myeloid_seurat,file = "GSE144735(All_tissue)_Myeloid cell_seurat_obj.rds")
    }

    ## Fibroblasts
    if (T) {
        Fibroblasts_seurat <- subset(seuratobj, subset = Cell_type == "Fibroblasts")

        Fibroblasts_seurat <- RunPCA(object = Fibroblasts_seurat)
        Fibroblasts_seurat <- FindNeighbors(Fibroblasts_seurat, dims = 1:15)
        Fibroblasts_seurat <- FindClusters(Fibroblasts_seurat, resolution = 0.9)
        Fibroblasts_seurat <- RunUMAP(Fibroblasts_seurat, dims = 1:15)

        ## rename clusters
        new.cluster.ids <- c(
            "Fibro_OGN", "eCAF", "Fibro_OGN", "mCAF", "Fibro_ADAMDEC1", "Fibro_GPM6B",
            "Pericyte", "myCAF_NOTCH3", "Fibro_OGN", "Fibro_LY6H", "Fibro_SGK1", "Fibro_OGN",
            "Fibro_BMP5", "iCAF", "myCAF_DES", "eCAF", "Fibro_ID1"
        )
        names(new.cluster.ids) <- levels(Fibroblasts_seurat)
        Fibroblasts_seurat <- RenameIdents(Fibroblasts_seurat, new.cluster.ids)
        Cell_subtype <- as.character(Fibroblasts_seurat@active.ident)
        Fibroblasts_seurat@meta.data$Cell_subtype <- Cell_subtype

        ## save the cell subtype information into whole seurat object
        seuratobj@meta.data$Cell_subtype[match(rownames(Fibroblasts_seurat@meta.data), rownames(seuratobj@meta.data))] <- Cell_subtype

        ## umap plot
        pdf(file = "./figure/annotation/GSE144735 Fibroblast annotation.pdf", width = 10, height = 7.5)
        p1 <- DimPlot(Fibroblasts_seurat, reduction = "umap", label = F, pt.size = 0.5, label.size = 5, combine = T) +
            theme(plot.title = element_text(hjust = 0.5))
        p2 <- DimPlot(Fibroblasts_seurat, reduction = "umap", label = F, pt.size = 0.5, label.size = 5, combine = T, split.by = "Class") +
            theme(plot.title = element_text(hjust = 0.5))
        p3 <- DimPlot(Fibroblasts_seurat, reduction = "umap", label = F, pt.size = 0.5, label.size = 5, combine = T, group.by = "Patient") +
            theme(plot.title = element_text(hjust = 0.5))
        print((p1 + p3) / p2)
        dev.off()

        ## makers dotplot
        pdf(file = "./figure/annotation/GSE144735 Fibroblast maker dotplot.pdf", width = 20, height = 7.5)
        if (1) {
            makers_list <- list()
            makers_list[["myCAF"]] <- c("MYL9", "TAGLN", "ACTA2", "MYH11", "MYLK", "ACTG2", "DES", "NOTCH3")
            makers_list[["iCAF"]] <- c("IL24", "MMP1", "MMP3", "CXCL6", "WNT5A")
            makers_list[["mCAF"]] <- c("CST1", "MMP11", "GRP")
            makers_list[["eCAF"]] <- c("MMP14", "LOXL2", "POSTN")
            makers_list[["Pericyte"]] <- c("RGS5", "CSPG4", "ABCC9", "KCNJ8")
            makers_list[["Fibroblast"]] <- c(
                "CCDC80", "OGN", "EFEMP1", "APOE", "CCL8", "ADAMDEC1",
                "CTSC", "BMP2", "BMP5", "FRZB", "GPM6B", "LY6H", "SGK1", "ID1"
            )
        }
        p <- DotPlot(Fibroblasts_seurat, features = makers_list) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        print(p)
        dev.off()
        ## save seurat object
        # saveRDS(Fibroblasts_seurat,file = "GSE144735(All_tissue)_Fibroblasts_seurat_obj.rds")
    }

    ## Endothelial cells
    if (T) {
        Endothelial_seurat <- subset(seuratobj, subset = Cell_type == "Endothelial cells")
        Endothelial_seurat <- RunPCA(object = Endothelial_seurat)
        Endothelial_seurat <- FindNeighbors(Endothelial_seurat, dims = 1:10)
        Endothelial_seurat <- FindClusters(Endothelial_seurat, resolution = 0.7)
        Endothelial_seurat <- RunUMAP(Endothelial_seurat, dims = 1:10)

        ## rename clusters
        new.cluster.ids <- c(
            "EC_PLVAP_PRCP", "EC_ACKR1", "EC_ACKR1", "EC_IGFBP3", "EC_GPIHBP1",
            "EC_ESM1", "EC_TFF3", "EC_PLVAP_CD36", "EC_Unspecific"
        )
        names(new.cluster.ids) <- levels(Endothelial_seurat)
        Endothelial_seurat <- RenameIdents(Endothelial_seurat, new.cluster.ids)
        Cell_subtype <- as.character(Endothelial_seurat@active.ident)
        Endothelial_seurat@meta.data$Cell_subtype <- Cell_subtype

        ## save the cell subtype information into whole seurat object
        seuratobj@meta.data$Cell_subtype[match(rownames(Endothelial_seurat@meta.data), rownames(seuratobj@meta.data))] <- Cell_subtype

        ## umap plot
        pdf(file = "./figure/annotation/GSE144735 Endothelial annotation.pdf", width = 12, height = 7.5)
        p1 <- DimPlot(Endothelial_seurat, reduction = "umap", label = F, pt.size = 0.5, label.size = 5, combine = T) +
            theme(plot.title = element_text(hjust = 0.5))
        p2 <- DimPlot(Endothelial_seurat, reduction = "umap", label = F, pt.size = 0.5, label.size = 5, combine = T, split.by = "Class") +
            theme(plot.title = element_text(hjust = 0.5))
        p3 <- DimPlot(Endothelial_seurat, reduction = "umap", label = F, pt.size = 0.5, label.size = 5, combine = T, group.by = "Patient") +
            theme(plot.title = element_text(hjust = 0.5))
        print((p1 + p3) / p2)
        dev.off()

        ## makers DotPlot
        pdf(file = "./figure/annotation/GSE144735 Endothelial maker dotplot.pdf", width = 15, height = 7.5)
        if (1) {
            makers_list <- list()
            makers_list[["EC_ACKR1"]] <- c("ACKR1", "SELP")
            makers_list[["EC_PLVAP"]] <- c("PLVAP", "PRCP", "CD36")
            makers_list[["EC_ESM1"]] <- c("ESM1", "APLN")
            makers_list[["EC_GPIHBP1"]] <- c("GPIHBP1", "FABP4")
            makers_list[["EC_TFF3"]] <- c("TFF3", "CCL21", "PROX1", "LYVE1")
            makers_list[["EC_IGFBP3"]] <- c("IGFBP3", "FBLN5")
        }
        p <- DotPlot(Endothelial_seurat, features = makers_list) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        print(p)
        dev.off()

        ## save seurat object
        # saveRDS(Endothelial_seurat,file = "GSE144735(All_tissue)_Endothelial_seurat_obj.rds")
    }

    ## Epithelial cells
    if (T) {
        ### subset Epithelial cells
        Epithelial_seurat <- subset(seuratobj, subset = Cell_type == "Epithelial cells")

        ## normal tissue annotation
        if (1) {
            Epithelial_seurat <- subset(Epithelial_seurat, subset = Class == "Normal")
            Epithelial_seurat <- RunPCA(object = Epithelial_seurat)
            Epithelial_seurat <- FindNeighbors(Epithelial_seurat, dims = 1:15)
            Epithelial_seurat <- FindClusters(Epithelial_seurat, resolution = 0.5)
            Epithelial_seurat <- RunUMAP(Epithelial_seurat, dims = 1:15)

            ## rename clusters
            new.cluster.ids <- c(
                "Intermediate", "Intermediate", "Stem-like(TA) cells", "Mature colonocytes", "BEST4+ colonocytes", "Goblet cells",
                "CDKN2B-AS1+ cells", "Goblet cells", "Tuft cells", "CHGA+ cells"
            )
            names(new.cluster.ids) <- levels(Epithelial_seurat)
            Epithelial_seurat <- RenameIdents(Epithelial_seurat, new.cluster.ids)
            Cell_subtype <- as.character(Epithelial_seurat@active.ident)

            ## save the cell subtype information into whole seurat object
            seuratobj@meta.data$Cell_subtype[match(rownames(Epithelial_seurat@meta.data), rownames(seuratobj@meta.data))] <- Cell_subtype

            ## umap plot
            pdf(file = "./figure/annotation/GSE144735 Epithelial(Normal) annotation.pdf", width = 10, height = 5)
            p1 <- DimPlot(Epithelial_seurat, reduction = "umap", label = F, pt.size = 0.5, label.size = 5, combine = T) +
                theme(plot.title = element_text(hjust = 0.5))
            p3 <- DimPlot(Epithelial_seurat, reduction = "umap", label = F, pt.size = 0.5, label.size = 5, combine = T, group.by = "Patient") +
                theme(plot.title = element_text(hjust = 0.5))
            print((p1 + p3))
            dev.off()

            ## makers DotPlot
            pdf(file = "./figure/annotation/GSE144735 Epithelial(Normal) maker dotplot.pdf", width = 15, height = 7.5)
            if (1) {
                makers_list <- list()
                makers_list[["Stem-like"]] <- c("LGR5", "OLFM4", "MKI67", "PCNA")
                makers_list[["Goblet"]] <- c("TFF3", "SPINK1", "REG4", "AGR2")
                makers_list[["Colonocytes"]] <- c("GUCA2B", "SLC26A3", "BEST4")
                makers_list[["Tuft"]] <- c("SOX9", "HPGDS", "PTGS1")
                makers_list[["Others"]] <- c("CDKN2B-AS1", "CHGA")
            }
            p <- DotPlot(Epithelial_seurat, features = makers_list) +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
            print(p)
            dev.off()
        }

        ## tumor and border tissue annotation
        if (1) {

            Epithelial_seurat <- subset(seuratobj, subset = Cell_type == "Epithelial cells")

            Epithelial_seurat <- subset(GSE144735_seurat, subset = Cell_type == "Epithelial cells")

            Epithelial_seurat <- subset(Epithelial_seurat, subset = (Class == "Tumor" | Class == "Border"))
            ## using CMSclassifier to annotate tumor and border cells via CMS
            cms_label <- CMSclassifier_annoFn(Epithelial_seurat)
            ## rename clusters
            Idents(Epithelial_seurat) <- cms_label
            Cell_subtype <- as.character(Epithelial_seurat@active.ident)
            ## save the cell subtype information into whole seurat object
            seuratobj@meta.data$Cell_subtype[match(rownames(Epithelial_seurat@meta.data), rownames(seuratobj@meta.data))] <- Cell_subtype
        }

        ## Overall
        Epithelial_seurat <- subset(seuratobj, subset = Cell_type == "Epithelial cells")
        Epithelial_seurat <- RunPCA(object = Epithelial_seurat)
        Epithelial_seurat <- FindNeighbors(Epithelial_seurat, dims = 1:15)
        Epithelial_seurat <- FindClusters(Epithelial_seurat, resolution = 0.7)
        Epithelial_seurat <- RunUMAP(Epithelial_seurat, dims = 1:15)

        Idents(Epithelial_seurat) <- Epithelial_seurat@meta.data$Cell_subtype

        ### umap plot
        pdf(file = "./figure/annotation/GSE144735 Epithelial(All) annotation.pdf", width = 12, height = 7.5)
        p1 <- DimPlot(Epithelial_seurat, reduction = "umap", label = F, pt.size = 0.5, label.size = 5, combine = T) +
            theme(plot.title = element_text(hjust = 0.5))
        p2 <- DimPlot(Epithelial_seurat, reduction = "umap", label = F, pt.size = 0.5, label.size = 5, combine = T, split.by = "Class") +
            theme(plot.title = element_text(hjust = 0.5))
        p3 <- DimPlot(Epithelial_seurat, reduction = "umap", label = F, pt.size = 0.5, label.size = 5, combine = T, group.by = "Patient") +
            theme(plot.title = element_text(hjust = 0.5))
        print((p1 + p3) / p2)
        dev.off()
    }

    ## save the annotation results
    seuratobj@meta.data$Cell_subtype[is.na(seuratobj@meta.data$Cell_subtype)] <- seuratobj@meta.data[is.na(seuratobj@meta.data$Cell_subtype), "Cell_type"]
    saveRDS(seuratobj, file = "./result/GSE144735_SeuratObj_anno.rds")
    write.table(seuratobj@meta.data,
        file = "./result/GSE144735_SeuratObj_annodata.txt",
        fileEncoding = "utf-8", sep = "\t"
    )

    return(NULL)
}

GSE132465_preProcessFn <- function() {
    exp <- read.table(file = "./data/GSE132465_exp.txt", header = T, row.names = 1, sep = "\t")
    clinical <- read.table(file = "./data/GSE132465_anno.txt", header = T, row.names = 1, sep = "\t", stringsAsFactors = F)

    GSE132465_seurat <- pre_processing_scRNAseq(exp, clinical)
    GSE132465_seurat <- FindNeighbors(GSE132465_seurat, dims = 1:12)
    GSE132465_seurat <- FindClusters(GSE132465_seurat, resolution = 0.3)

    GSE132465_seurat <- RunUMAP(GSE132465_seurat, dims = 1:12)

    GSE132465_seurat@meta.data$Cell_type <- NA
    GSE132465_seurat@meta.data$Cell_subtype <- NA

    ## save the un-annotation seurat object
    saveRDS(GSE132465_seurat, file = "./result/GSE132465_SeuratObj.rds")
    return(NULL)
}

GSE132257_preProcessFn <- function() {
    exp <- read.table(file = "./data/GSE132257_exp.csv", header = T, sep = ",")
    exp <- as.data.frame(t(exp))

    ### load clinical
    clinical <- read.table(file = "./data/GSE132257_anno.csv", header = T, sep = ",")

    GSE132257_seurat <- pre_processing_scRNAseq(exp, clinical)
    GSE132257_seurat <- FindNeighbors(GSE132257_seurat, dims = 1:10)
    seurat_test <- FindClusters(object = GSE132257_seurat, resolution = c(seq(0.1, 1.6, 0.2)))
    GSE132257_seurat <- FindClusters(GSE132257_seurat, resolution = 0.9)
    GSE132257_seurat <- RunUMAP(GSE132257_seurat, dims = 1:10)

    ## add cluster id into metadata
    GSE132257_seurat@meta.data$Cell_type <- NA
    GSE132257_seurat@meta.data$Cell_subtype <- NA

    ### save the un-annotation seurat object
    saveRDS(GSE132257_seurat, file = "./result/GSE132257_SeuratObj.rds")
    return(NULL)

}

## load bulk RNA-seq expression profiles and clinical dataframe
load_bulk_Exp <- function(path="same") {
    if (path=="same"){
        exp_matrix_files <- list.files("./data/bulkExp")
        names <- sapply(strsplit(x = exp_matrix_files, split = "[.]"), function(x) x[1])
        names <- sapply(strsplit(x = names, split = "_"), function(x) {
            x[1]
        })
        exp_matrix_list <- list()
        for (i in 1:length(exp_matrix_files)) {
            exp_matrix_list[[names[i]]] <- readRDS(paste0(
                "./data/bulkExp/",
                exp_matrix_files[i]
            ))
        }
        return(exp_matrix_list)
    }
    if (path=="different"){
        exp_matrix_files <- list.files("../../data/bulkExp")
        names <- sapply(strsplit(x = exp_matrix_files, split = "[.]"), function(x) x[1])
        names <- sapply(strsplit(x = names, split = "_"), function(x) {
            x[1]
        })
        exp_matrix_list <- list()
        for (i in 1:length(exp_matrix_files)) {
            exp_matrix_list[[names[i]]] <- readRDS(paste0(
                "../../data/bulkExp/",
                exp_matrix_files[i]
            ))
        }
        return(exp_matrix_list)
    }        

}



load_bulk_Clinical <- function(path="same") {
    if (path=="same"){
        clinical_files <- list.files("./data/bulkClinical")
        names <- sapply(strsplit(x = clinical_files, split = "[.]"), function(x) x[1])
        names <- sapply(strsplit(x = names, split = "_"), function(x) {
            x[2]
        })
        clinic_list <- list()

        for (i in 1:length(clinical_files)) {
            clinic_list[[names[i]]] <- read.table(paste0(
                "./data/bulkClinical/",
                clinical_files[i]
            ), sep = "\t", header = T, stringsAsFactors = F, row.names = 1)

            if (dim(clinic_list[[names[i]]])[2] == 3) {
                clinic_list[[names[i]]] <- clinic_list[[names[i]]][, 2:3]
            }
        }
        return(clinic_list)
    }

    if (path=="different"){
        clinical_files <- list.files("../../data/bulkClinical")
        names <- sapply(strsplit(x = clinical_files, split = "[.]"), function(x) x[1])
        names <- sapply(strsplit(x = names, split = "_"), function(x) {
            x[2]
        })
        clinic_list <- list()

        for (i in 1:length(clinical_files)) {
            clinic_list[[names[i]]] <- read.table(paste0(
                "../../data/bulkClinical/",
                clinical_files[i]
            ), sep = "\t", header = T, stringsAsFactors = F, row.names = 1)

            if (dim(clinic_list[[names[i]]])[2] == 3) {
                clinic_list[[names[i]]] <- clinic_list[[names[i]]][, 2:3]
            }
        }
        return(clinic_list)
    }
}

## load raw clinical data
load_bulk_RawClinical <- function() {
    raw_clinical <- list()
    raw_clinical_files <- list.files("./data/bulkRawClinical/")
    names <- sapply(strsplit(x = raw_clinical_files, split = "[.]"), function(x) x[1])
    names <- sapply(strsplit(x = names, split = "_"), function(x) {
        x[1]
    })
    for (i in 1:length(raw_clinical_files)) {
        raw_clinical[[names[i]]] <- read.table(paste0(
            "./data/bulkRawClinical/",
            raw_clinical_files[i]
        ), sep = "\t", header = T, stringsAsFactors = F, row.names = 1)
    }
    return(raw_clinical)
}

## seperate all datasets into training, testing and validating sets
combine_Datasets <- function(exp_matrix_list, clinic_list) {
    all_exp <- matrix(nrow = nrow(exp_matrix_list[[1]]), ncol = 0)
    for (i in 1:length(exp_matrix_list)) {
        all_exp <- cbind(all_exp, exp_matrix_list[[i]])
    }

    all_clinical <- matrix(ncol = ncol(clinic_list[[1]]), nrow = 0)
    for (i in 1:length(clinic_list)) {
        all_clinical <- rbind(all_clinical, clinic_list[[i]])
    }

    return(list(all_exp, all_clinical))
}

split_Datasets <- function(exp_matrix_list, clinic_list, validate_id) {
    ## split training dataset and test-set and validate-set
    all_exp <- matrix(nrow = nrow(exp_matrix_list[[1]]), ncol = 0)
    for (i in 1:length(exp_matrix_list)) {
        if (i %in% validate_id) {
            next
        }
        all_exp <- cbind(all_exp, exp_matrix_list[[i]])
    }

    all_clinical <- matrix(ncol = ncol(clinic_list[[1]]), nrow = 0)
    for (i in 1:length(clinic_list)) {
        if (i %in% validate_id) {
            next
        }
        all_clinical <- rbind(all_clinical, clinic_list[[i]])
    }

    training_exp <- all_exp
    training_clinical <- all_clinical

    ## training the models
    ## traning set size : testing set size= 8:2
    training_size <- round(digits = 0, x = nrow(all_clinical) / 10 * 8)
    test_size <- nrow(all_clinical) - training_size

    return(list(training_exp, training_clinical, training_size, test_size))
}

## processing clinical datasets and applied multi-cox
MultiCox_Data_Transform <- function(path_, validation_sets, exp_matrix_list, clinic_list, raw_clinical) {
    interations <- list.files(path_)
    if (length(interations) == 0) {
        return("next")
    } else {
        interations <- sapply(interations, function(x) {
            strsplit(x, "tion")[[1]][2]
        })
        interations <- as.numeric(interations)
        names(interations) <- NULL
    }
    for (iteration in interations) {
        label1 <- 0
        label2 <- 0
        label3 <- 0
        setwd(paste0(path_, "/", "iteration", iteration, "/"))
        model_path <- paste0(
            "./Traning set(size=367) prognostic-signature of ",
            celltype, "(ncells=10)_model.rds"
        )

        train_model <- readRDS(model_path)
        cutoff <- train_model$cutoff

        vali_exp <- list()
        vali_clinical <- list()
        raw_clinical_list <- list()

        for (i in validation_sets) {
            vali_exp[[i]] <- exp_matrix_list[[i]]
            vali_clinical[[i]] <- clinic_list[[i]]
            raw_clinical_list[[i]] <- raw_clinical[[i]]
            raw_clinical_list[[i]] <- raw_clinical_list[[i]][match(rownames(vali_clinical[[i]]), raw_clinical[[i]]$Patient_ID), ]
        }

        ## calculate gene-pairs rank
        for (i in validation_sets) {
            label_tem <- grouping(vali_exp[[i]], vali_clinical[[i]], cutoff)
            vali_clinical[[i]] <- cbind(vali_clinical[[i]], label_tem)
        }

        ## pre-processing the clinical files
        if (1) {
            ## Age,Gender,Location,Mutations,Stage to analysis
            ## GSE144735
            if (1) {
                head(raw_clinical_list[[1]], 5)
                raw_clinical_list[[1]]$Age <- ifelse(raw_clinical_list[[1]]$Age >= 70, ">=70", "<70")
                raw_clinical_list[[1]]$Gender <- ifelse(raw_clinical_list[[1]]$Gender == " F", "Female", "Male")
                raw_clinical_list[[1]] <- raw_clinical_list[[1]][, c(1:5)]

                vali_clinical[[1]]$Patient_ID <- rownames(vali_clinical[[1]])

                vali_clinical[[1]] <- full_join(vali_clinical[[1]], raw_clinical_list[[1]], by = "Patient_ID")
                rownames(vali_clinical[[1]]) <- vali_clinical[[1]]$Patient_ID
                vali_clinical[[1]] <- vali_clinical[[1]][, -4]
            }
            ## GSE17536
            if (1) {
                head(raw_clinical_list[[2]], 5)
                raw_clinical_list[[2]]$Age <- ifelse(raw_clinical_list[[2]]$Age >= 70, ">=70", "<70")
                raw_clinical_list[[2]]$Gender <- ifelse(raw_clinical_list[[2]]$Gender == "female", "Female", "Male")
                raw_clinical_list[[2]] <- raw_clinical_list[[2]][, c(1:3, 8:10)]

                vali_clinical[[2]]$Patient_ID <- rownames(vali_clinical[[2]])

                vali_clinical[[2]] <- full_join(vali_clinical[[2]], raw_clinical_list[[2]], by = "Patient_ID")
                rownames(vali_clinical[[2]]) <- vali_clinical[[2]]$Patient_ID
                vali_clinical[[2]] <- vali_clinical[[2]][, -4]
            }
            ## GSE39582
            if (1) {
                head(raw_clinical_list[[3]], 5)
                raw_clinical_list[[3]]$Age <- ifelse(raw_clinical_list[[3]]$Age >= 70, ">=70", "<70")

                raw_clinical_list[[3]] <- raw_clinical_list[[3]][, c(1:3, 8:9, 14, 17, 19:20)]

                vali_clinical[[3]]$Patient_ID <- rownames(vali_clinical[[3]])

                vali_clinical[[3]] <- full_join(vali_clinical[[3]], raw_clinical_list[[3]], by = "Patient_ID")
                rownames(vali_clinical[[3]]) <- vali_clinical[[3]]$Patient_ID
                vali_clinical[[3]] <- vali_clinical[[3]][, -4]
            }
            rm(vali_exp, raw_clinical_list)
        }

        ## multi-cox
        for (i in 1:length(vali_clinical)) {
            vali_clinical[[i]] <- remove_less1_level(vali_clinical[[i]])
        }
        setwd("../../../")
        return(vali_clinical)
    }
}

## remove factor which less than 1
remove_less1_level <- function(vali_df) {
    cat(colnames(vali_df), "\n")
    for (i in 3:ncol(vali_df)) {
        x <- table(vali_df[, i])
        print(table(vali_df[, i]))
        x <- as.data.frame(x)
        delete <- c(as.character(x$Var1)[which(x$Freq < 2)])
        if (length(delete)) {
            cat(as.character(x$Var1)[which(x$Freq < 2)], "is less than 2, delete it.", "\n")
            index <- !vali_df[, i] %in% delete
            vali_df <- vali_df[index, ]
        }
        vali_df[, i] <- as.factor(vali_df[, i])
    }
    return(vali_df)
}
