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
seurat.integrated.rds <- here("results/seurat/20230310_PGP1_ALL_cells_UNCR3_and_Bhaduri_reference_integrated_seurat_object.rds")

# INPUT FILES ##########################################################################################################
# merged seurat filtered file
seurat.merged.filtered.rds <- here("results/seurat/20230227_PGP1_QC_filtered_seurat_object.rds")

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
# df.samples %>%
#     filter(cell_count > 200) %>%
#     pull(Sample) -> samples.to.keep

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

# sample to 10% cells per samples (unless thats under 200)
# percentCells <- 0.1
# printMessage(paste("Sampling to", percentCells*100, "percent cells per sample"))
# seur.list <- lapply(X = seur.list, FUN = function(x) {
#     numCells <- length(colnames(x))
#     targetCells <- numCells * percentCells
#     if (targetCells < 200) {
#         x <- x[, sample(colnames(x), size = 200, replace = FALSE)]
#     } else {
#         x <- x[, sample(colnames(x), size = targetCells, replace = FALSE)]
#     }
# })

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

# run PCA
printMessage("Running PCA")
seur.list <- lapply(X = seur.list, FUN = function(x) {
    x <- RunPCA(x, features = features, verbose = FALSE)
})

# find anchors
printMessage("Finding integration anchors")
pgp1.anchors <- FindIntegrationAnchors(object.list = seur.list,
                                       reference = c(2,7),
                                       reduction = "cca",
                                       dims = 1:50,
                                       normalization.method = "SCT",
                                       anchor.features = features)

# integrate data
printMessage("Integrating data")
seur.integrated <- IntegrateData(anchorset = pgp1.anchors,
                                 normalization.method = "SCT",
                                 dims = 1:50)

# save integrated data
printMessage("Saving integrated data")
saveRDS(seur.integrated, seurat.integrated.rds)
printMessage("Done saving data")



















