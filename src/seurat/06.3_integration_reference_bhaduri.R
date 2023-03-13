# integration analysis

library(Seurat)
library(glmGamPoi)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(here)

library(mikelaffr)

# OUTPUT FILES #########################################################################################################
# integrated dataset
seurat.integrated.bhaduriREF.rds <- here("results/seurat/20230312_PGP1_ALL_cells_Bhaduri_reference_integrated_seurat_object.rds")

# INPUT FILES ##########################################################################################################
# merged seurat filtered file (QC filtered and doublet filtered)
seurat.pgp1.rds <- here("results/seurat/20230227_PGP1_QC_DF_filtered_seurat_object.rds")

# Bhaduri (kreigstein) seurat object
seur.bhaduri.rds <- "/proj/steinlab/projects/IVIV_scRNA/youngsook_pine/extData_annotation/kreigstein_primary_seuratObject_SCTv2/Kreigstein_primary_integrated.rds"

# gencode gtf file
#gencode.gtf <- here("data/refgenome/gencode/gencode.v40.annotation.gtf.gz")
# gencode gene names, types, and ensgid
df.gencode.csv <- here("data/refgenome/gencode/gencode.v40.genes.csv.gz")

# GLOBALS ##############################################################################################################
#



# Load Seurat Objects ########
printMessage("Loading PGP1 Seurat Object...")
seur.pgp1 <- readRDS(seurat.pgp1.rds)
dim(seur.pgp1)

printMessage("Loading Bhaduri Seurat Object...")
seur.bhaduri <- readRDS(seur.bhaduri.rds)
dim(seur.bhaduri)




# split object into list by sample
seur.list <- SplitObject(seur.filtered, split.by = "Sample")

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



















