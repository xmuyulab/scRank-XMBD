## This part is to evaluate the REO of some gene-pairs in bulk RNA-seq dataset with ground truth # nolint

getwd()
workdir <- "./" # nolint
setwd(workdir)

library(ggplot2)
library(ggpubr)
library(ggthemes)

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
} ## The cell sub-population and major population

## CRC bulk ground truth
bulk_ID <- "GSE39396"
bulk_exp <- readRDS(paste0("./data/", bulk_ID, "_exp.rds"))
bulk_clinical <- readRDS(paste0("./data/", bulk_ID, "_clinical.rds"))

bulk_exp[1:5, 1:5]
head(bulk_clinical)
print(table(bulk_clinical$cell_population))

## gene pairs
genePairs <- readRDS("./data/Stable-cell-subtype-specific gene pairs(stable_cutoff=0.7,ncells=10,deltaS=0.5).rds") # nolint
class(genePairs)
class(genePairs[["GSE144735"]])
class(genePairs[["GSE144735"]][["CMS3"]])

names(genePairs[["GSE144735"]])

major_df <- load_major_celltype_name()

mergePair <- function(list1, list2) {
    new_list <- list()

    subtypeID <- union(names(list1), names(list2))

    for (subtype in subtypeID) {
        Gene_pair <- union(names(list1[[subtype]]), names(list2[[subtype]]))
        new_list[[subtype]] <- Gene_pair
        print(paste0("The number of gene pairs in ", subtype, " is ", length(Gene_pair))) # nolint
    }
    return(new_list)
}

allPair <- mergePair(genePairs[[1]], genePairs[[2]])

ToMajorPair <- function(list_, major_df) {
    subtypeID <- names(list_)
    new_list <- list()

    for (subtype in subtypeID) {
        major_tem <- major_df[match(subtype, major_df[, 1]), 2]
        new_list[[major_tem]] <- union(new_list[[major_tem]], list_[[subtype]])

        print(paste0("The number of gene pairs in ", major_tem, " is ", length(new_list[[major_tem]]))) # nolint
    }
    return(new_list)
}

allPair <- ToMajorPair(allPair, major_df)

splitPair <- function(list_) {
    new_list <- list()

    for (i in names(list_)) {
        tem <- list_[[i]]
        df_tem <- matrix(data = NA, ncol = 2, nrow = length(tem))
        gene1 <- sapply(tem, function(x) {
            return(strsplit(x, split = "_")[[1]][1])
        })
        gene2 <- sapply(tem, function(x) {
            return(strsplit(x, split = "_")[[1]][2])
        })

        df_tem[, 1] <- gene1
        df_tem[, 2] <- gene2
        colnames(df_tem) <- c("gene1", "gene2")
        new_list[[i]] <- df_tem
    }
    return(new_list)
}

allPair <- splitPair(allPair)

## epithelial cell
ReoInBulk <- function(exp_, clinical_, genePair, celltype) {

    ## get the gene pair of certain celltyep
    genePair <- genePair[[celltype]]

    allgenes <- c()

    allgenes <- as.character(genePair)
    allgenes <- unique(allgenes)

    ## get the exp
    tem <- exp_[match(allgenes, rownames(exp_)), ]
    print(dim(tem))

    plot_df <- matrix(data = NA, nrow = nrow(genePair) * ncol(tem) * 2, ncol = 5)
    colnames(plot_df) <- c("expression", "gene_name", "gene_group", "patient", "celltype") # nolint
    plot_df <- as.data.frame(plot_df)
    print(dim(plot_df))

    ## stack by each column
    plot_df[, 2] <- c(rep(genePair[, 1], times = nrow(clinical_)), rep(genePair[, 2], times = nrow(clinical_))) ## gene_name # nolint
    plot_df[, 3] <- rep(c("gene1", "gene2"), each = nrow(genePair) * ncol(tem)) ## gene group # nolint
    plot_df[, 4] <- rep(rep(rownames(clinical_), each = nrow(genePair)), times = 2) ## patient # nolint
    plot_df[, 5] <- clinical_[match(plot_df[, 4], clinical_[, "geo_accession"]), "cell_population"] ## celltype # nolint

    for (i in 1:nrow(plot_df)) {
        plot_df[i, 1] <- tem[match(plot_df[i, 2], rownames(tem)), match(plot_df[i, 4], colnames(tem))] ## expression # nolint
    }
    print(dim(plot_df))
    plot_df <- na.omit(plot_df)
    print(dim(plot_df))

    p <- ggplot(plot_df, aes(celltype, expression)) +
        geom_boxplot(aes(color = gene_group)) +
        scale_fill_discrete(labels = c("Gene1", "Gene2")) +
        scale_color_manual(values = c("#b31919", "#1d1da0")) +
        theme_wsj(color = "white") +
        stat_compare_means(aes(group = gene_group), label = "p.format")

    ggsave(paste0("./figure/",celltype, " REOs in GSE39396.pdf"), plot = p, width = 10, height = 8, dpi = 150, units = "in", device = "pdf") # nolint

    return(NULL)
}

print(names(allPair))

for (i in names(allPair)) {
    ReoInBulk(bulk_exp, bulk_clinical, allPair, i)
}

