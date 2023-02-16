# integration analysis


#library(fishpond)
#library(data.table)
library(Seurat)
#library(miQC)
library(SeuratWrappers)
library(DoubletFinder)
library(glmGamPoi)
library(flexmix)
library(SingleCellExperiment)
library(Matrix)
library(stringr)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(tidytext)
library(here)

#library(plyranges)

library(mikelaffr)

# OUTPUT FILES #########################################################################################################
# integrated dataset
seurat.integrated.rds <- here("results/seurat/20230214_PGP1_10percent_of_cells_integrated_seurat_object.rds")

# INPUT FILES ##########################################################################################################
# merged seurat filtered file
seurat.merged.filtered.rds <- here("results/seurat/20230201_PGP1_filtered_seurat_object.rds")

# gencode gtf file
#gencode.gtf <- here("data/refgenome/gencode/gencode.v40.annotation.gtf.gz")
# gencode gene names, types, and ensgid
df.gencode.csv <- here("data/refgenome/gencode/gencode.v40.genes.csv.gz")

# GLOBALS ##############################################################################################################
#



# Load Seurat Object ########
printMessage("Loading Filtered Seurat Object...")
seur.filtered <- readRDS(seurat.merged.filtered.rds)
dim(seur.filtered)
seur.filtered@meta.data$Sample <- paste(seur.filtered@meta.data$Site, seur.filtered@meta.data$Day, seur.filtered@meta.data$Rep, sep = "_")

# filter for at least 100 cells per sample

df.cellData <- as_tibble(seur.filtered@meta.data)

df.cellData %>%
    group_by(Sample, WTK_ID, Site, Day, Rep) %>%
    dplyr::summarise(cell_count = dplyr::n(),
                     num_wells = n_distinct(well)) -> df.samples

df.samples %>%
    filter(cell_count > 100) %>%
    pull(Sample) -> samples.to.keep


seur.filtered <- subset(seur.filtered, subset = Sample %in% samples.to.keep)
dim(seur.filtered)

# split object into list by sample
seur.list <- SplitObject(seur.filtered, split.by = "Sample")

# sample to x cells per samples
# numCells <- 400
# printMessage(paste("Sampling to", numCells, "cells per sample"))
# seur.list <- lapply(X = seur.list, FUN = function(x) {
#     if (length(colnames(x)) < numCells) {
#         x <- x
#     } else {
#         x <- x[, sample(colnames(x), size = numCells, replace = FALSE)]
#     }
# })

# sample to 10% cells per samples (unless thats under 100)
percentCells <- 0.1
printMessage(paste("Sampling to", percentCells*100, "percent cells per sample"))
seur.list <- lapply(X = seur.list, FUN = function(x) {
    numCells <- length(colnames(x))
    targetCells <- numCells * percentCells
    if (targetCells < 100) {
        x <- x[, sample(colnames(x), size = 100, replace = FALSE)]
    } else {
        x <- x[, sample(colnames(x), size = targetCells, replace = FALSE)]
    }
})

# SCTransform each object in list
printMessage("SCTransform")
seur.list <- lapply(X = seur.list, FUN = function(x) {
    x <- SCTransform(x,
                     vars.to.regress = c("percent.mt"),
                     return.only.var.genes = FALSE,
                     vst.flavor = "v2",
                     method = 'glmGamPoi',
                     verbose = FALSE)
})

# select features that are variable across datasets
printMessage("Selecting integration features")
features <- SelectIntegrationFeatures(object.list = seur.list)

# prep integration
printMessage("Prepping integration")
seur.list <- PrepSCTIntegration(object.list = seur.list, anchor.features = features)

# find anchors
printMessage("Finding integration anchors")
pgp1.anchors <- FindIntegrationAnchors(object.list = seur.list, normalization.method = "SCT", anchor.features = features)

# integrate data
printMessage("Integrating data")
seur.integrated <- IntegrateData(anchorset = pgp1.anchors, normalization.method = "SCT")

# save integrated data
printMessage("Saving integrated data")
saveRDS(seur.integrated, seurat.integrated.rds)
printMessage("Done saving data")



# DefaultAssay(seur.integrated) <- "integrated"
#
# pdf("integrated_ALLcells.pdf", height = 7, width = 10)
#
#
# seur.integrated <- RunPCA(seur.integrated, verbose = TRUE)
# ElbowPlot(seur.integrated)
#
# DimPlot(seur.integrated, reduction = "pca", group.by = "Day")
# DimPlot(seur.integrated, reduction = "pca", group.by = "Site")
# DimPlot(seur.integrated, reduction = "pca", group.by = "WTK_ID")
# FeaturePlot(seur.integrated, reduction = "pca", features = c("percent.mt"))
#
# seur.integrated <- RunUMAP(seur.integrated, dims = 1:10, verbose = TRUE)
#
# seur.integrated <- FindNeighbors(seur.integrated, dims = 1:10, verbose = TRUE)
# seur.integrated <- FindClusters(seur.integrated, verbose = TRUE, resolution = 0.8)
#
# DimPlot(seur.integrated,
#         label = FALSE)
#
# DimPlot(seur.integrated,
#         group.by = "Day",
#         label = FALSE)
# DimPlot(seur.integrated,
#         group.by = "Site",
#         label = FALSE)
# DimPlot(seur.integrated,
#         group.by = "WTK_ID",
#         label = FALSE)
# DimPlot(seur.integrated,
#         group.by = "Sample",
#         label = FALSE)
#
# ptsize <- NULL
#
# DefaultAssay(seur.integrated) <- "SCT"
# FeaturePlot(seur.integrated, reduction = "umap", features = c("percent.mt"), pt.size = ptsize)
#
# FeaturePlot(seur.integrated, reduction = "umap", features = c("ENSG00000181449.4"), pt.size = ptsize) + ggtitle("SOX2")
# FeaturePlot(seur.integrated, reduction = "umap", features = c("ENSG00000148773.14"), pt.size = ptsize) + ggtitle("MKI67")
# FeaturePlot(seur.integrated, reduction = "umap", features = c("ENSG00000007372.25"), pt.size = ptsize) + ggtitle("PAX6")
#
# FeaturePlot(seur.integrated, reduction = "umap", features = c("ENSG00000136535.15"), pt.size = ptsize) + ggtitle("TBR1")
# FeaturePlot(seur.integrated, reduction = "umap", features = c("ENSG00000163508.13"), pt.size = ptsize) + ggtitle("EOMES")
#
# FeaturePlot(seur.integrated, reduction = "umap", features = c("ENSG00000171476.22"), pt.size = ptsize) + ggtitle("HOPX")
# FeaturePlot(seur.integrated, reduction = "umap", features = c("ENSG00000134853.12"), pt.size = ptsize) + ggtitle("PDGFRA")
#
# FeaturePlot(seur.integrated, reduction = "umap", features = c("ENSG00000077279.20"), pt.size = ptsize) + ggtitle("DCX")
# FeaturePlot(seur.integrated, reduction = "umap", features = c("ENSG00000258947.8"), pt.size = ptsize) + ggtitle("TUBB3")
# FeaturePlot(seur.integrated, reduction = "umap", features = c("ENSG00000164600.7"), pt.size = ptsize) + ggtitle("NEUROD6")
#
# FeaturePlot(seur.integrated, reduction = "umap", features = c("ENSG00000127152.18"), pt.size = ptsize) + ggtitle("BCL11B")
# FeaturePlot(seur.integrated, reduction = "umap", features = c("ENSG00000119042.17"), pt.size = ptsize) + ggtitle("SATB2")
#
# FeaturePlot(seur.integrated, reduction = "umap", features = c("ENSG00000204531.21"), pt.size = ptsize) + ggtitle("OCT4")
# FeaturePlot(seur.integrated, reduction = "umap", features = c("ENSG00000144355.15"), pt.size = ptsize) + ggtitle("DLX1")
# FeaturePlot(seur.integrated, reduction = "umap", features = c("ENSG00000115844.11"), pt.size = ptsize) + ggtitle("DLX2")
# FeaturePlot(seur.integrated, reduction = "umap", features = c("ENSG00000231764.11"), pt.size = ptsize) + ggtitle("DLX6-AS1")
#
# FeaturePlot(seur.integrated, reduction = "umap", features = c("ENSG00000118271.12"), pt.size = ptsize) + ggtitle("TTR")
# FeaturePlot(seur.integrated, reduction = "umap", features = c("ENSG00000176165.12"), pt.size = ptsize) + ggtitle("FOXG1")
# FeaturePlot(seur.integrated, reduction = "umap", features = c("ENSG00000180613.11"), pt.size = ptsize) + ggtitle("GSX2")
# FeaturePlot(seur.integrated, reduction = "umap", features = c("ENSG00000106852.16"), pt.size = ptsize) + ggtitle("LHX6")
#
#
# dev.off()















