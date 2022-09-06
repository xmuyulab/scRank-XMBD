## part2.0 extract tumor object
library(Seurat)
library(clustree)
GSE132257_seurat <- readRDS("../../result/GSE132257_SeuratObj_anno.rds")
GSE144735_seurat <- readRDS("../../result/GSE144735_SeuratObj_anno.rds")
GSE132465_seurat <- readRDS("../../result/GSE132465_SeuratObj_anno.rds")

set.seed(619)
## GSE132257
GSE132257_seurat <- subset(GSE132257_seurat, subset = (Class == "Tumor" | Class == "Border"))
GSE132257_seurat <- FindVariableFeatures(GSE132257_seurat, selection.method = "vst", nfeatures = 2000)
GSE132257_seurat <- ScaleData(object = GSE132257_seurat) ## 中心化
GSE132257_seurat <- RunPCA(object = GSE132257_seurat) ## 对于差异表达基因进行PCA降维
seurat_obj2 <- FindNeighbors(GSE132257_seurat, dims = 1:10)
seurat_obj2 <- FindClusters(seurat_obj2, resolution = 0.7)
seurat_obj_tsne <- RunTSNE(seurat_obj2, dims = 1:10)
seurat_obj_tsne@meta.data$cell_typeid <- as.character(seurat_obj_tsne@active.ident)
saveRDS(seurat_obj_tsne, file = "./data/GSE132257(Tumor)_seurat_obj.rds")


## GSE144735
GSE144735_seurat <- subset(GSE144735_seurat, subset = (Class == "Tumor" | Class == "Border"))
GSE144735_seurat <- FindVariableFeatures(GSE144735_seurat, selection.method = "vst", nfeatures = 2000)
GSE144735_seurat <- ScaleData(object = GSE144735_seurat) ## 中心化
GSE144735_seurat <- RunPCA(object = GSE144735_seurat)
seurat_obj2 <- FindNeighbors(GSE144735_seurat, dims = 1:10)
seurat_obj2 <- FindClusters(seurat_obj2, resolution = 0.9)
seurat_obj_tsne <- RunTSNE(seurat_obj2, dims = 1:10)
seurat_obj_tsne@meta.data$cell_typeid <- as.character(seurat_obj_tsne@active.ident)
saveRDS(seurat_obj_tsne, file = "./data/GSE144735(Tumor)_seurat_obj.rds")

## GSE132465
GSE132465_seurat <- subset(GSE132465_seurat, subset = (Class == "Tumor" | Class == "Border"))
GSE132465_seurat <- FindVariableFeatures(GSE132465_seurat, selection.method = "vst", nfeatures = 2000)
GSE132465_seurat <- ScaleData(object = GSE132465_seurat)
GSE132465_seurat <- RunPCA(object = GSE132465_seurat)
seurat_obj2 <- FindNeighbors(GSE132465_seurat, dims = 1:10)
seurat_obj2 <- FindClusters(seurat_obj2, resolution = 0.5)
seurat_obj_tsne <- RunTSNE(seurat_obj2, dims = 1:10)
seurat_obj_tsne@meta.data$cell_typeid <- as.character(seurat_obj_tsne@active.ident)
saveRDS(seurat_obj_tsne, file = "./data/GSE132465(Tumor)_seurat_obj.rds")

## for GSE132465 split into to part
seurat_obj_tsne <- readRDS("./data/GSE132465(Tumor)_seurat_obj.rds")

seurat_obj_tsne@meta.data$orig.ident <- as.character(seurat_obj_tsne@meta.data$orig.ident)

seurat_obj_tsne_1 <- subset(
  seurat_obj_tsne,
  orig.ident == "SMC01.T" |
    orig.ident == "SMC02.T" |
    orig.ident == "SMC03.T" |
    orig.ident == "SMC04.T" |
    orig.ident == "SMC05.T" |
    orig.ident == "SMC06.T" |
    orig.ident == "SMC07.T" |
    orig.ident == "SMC08.T" |
    orig.ident == "SMC09.T" |
    orig.ident == "SMC10.T" |
    orig.ident == "SMC11.T"
)

seurat_obj_tsne_2 <- subset(
  seurat_obj_tsne,
  orig.ident == "SMC14.T" |
    orig.ident == "SMC15.T" |
    orig.ident == "SMC16.T" |
    orig.ident == "SMC17.T" |
    orig.ident == "SMC18.T" |
    orig.ident == "SMC19.T" |
    orig.ident == "SMC20.T" |
    orig.ident == "SMC21.T" |
    orig.ident == "SMC22.T" |
    orig.ident == "SMC23.T" |
    orig.ident == "SMC24.T" |
    orig.ident == "SMC25.T"
)

seurat_obj_tsne_1 <- FindNeighbors(seurat_obj_tsne_1, dims = 1:10)
seurat_obj_tsne_2 <- FindNeighbors(seurat_obj_tsne_2, dims = 1:10)

saveRDS(seurat_obj_tsne_1, "./data/GSE132465_1(Tumor)_seurat_obj.rds")
saveRDS(seurat_obj_tsne_2, "./data/GSE132465_2(Tumor)_seurat_obj.rds")
