# plot UMAP



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

library(plyranges)

library(mikelaffr)

# OUTPUT FILES #########################################################################################################
dir.pdf <- here("doc/seurat/pdf/")
dir.create(dir.pdf, showWarnings = FALSE, recursive = TRUE)

# INPUT FILES ##########################################################################################################
# merged and filtered, sctransformed with sample regressed
seurat.transformed.rds <- here("results/seurat/20230205_PGP1_QCfiltered_SCTransform_sample_regressed_seurat_object.rds")

# integrated dataset 10% of cells
seurat.transformed.rds <- here("results/seurat/20230214_PGP1_10percent_of_cells_integrated_seurat_object.rds")

# gencode gene names, types, and ensgid
df.gencode.csv <- here("data/refgenome/gencode/gencode.v40.genes.csv.gz")

# GLOBALS ##############################################################################################################
#

# load gene names
df.gencode <- read_csv(df.gencode.csv)

printMessage("Loading Seurat Object...")
seur.filtered <- readRDS(seurat.transformed.rds)
printMessage("Finished Loading.")

# Subsample Cells
# seur.sample <- seur.filtered[, sample(colnames(seur.filtered), size = 10000, replace = FALSE)]
# rm(seur.filtered)
# seur.filtered <- seur.sample
# rm(seur.sample)

DefaultAssay(seur.filtered) <- "integrated"
seur.filtered <- RunPCA(seur.filtered, verbose = TRUE)
seur.filtered <- RunUMAP(seur.filtered, dims = 1:10, verbose = TRUE)

seur.filtered <- FindNeighbors(seur.filtered, dims = 1:10, verbose = TRUE)
seur.filtered <- FindClusters(seur.filtered, verbose = TRUE, resolution = 0.8)

pdf(paste0(dir.pdf, "20230220_10percent_cells_full_integrated_0.8_clustered.pdf"), height = 5, width = 7)

ElbowPlot(seur.filtered)

printMessage("Plotting PCA")
DimPlot(seur.filtered, reduction = "pca", group.by = "Day")
DimPlot(seur.filtered, reduction = "pca", group.by = "Site")
DimPlot(seur.filtered, reduction = "pca", group.by = "WTK_ID")
FeaturePlot(seur.filtered, reduction = "pca", features = c("percent.mt"))


DefaultAssay(seur.filtered) <- "SCT"
printMessage("Plotting UMAPs")
#pdf("scaled_transformed_resolution_0.8.pdf", height = 7, width = 10)
DimPlot(seur.filtered,
        label = FALSE)

DimPlot(seur.filtered,
        group.by = "Day",
        label = FALSE)
DimPlot(seur.filtered,
        group.by = "Site",
        label = FALSE)
DimPlot(seur.filtered,
        group.by = "WTK_ID",
        label = FALSE)
#ggsave("iddrc_umap_by_site.pdf", height = 5, width = 7)
DimPlot(seur.filtered,
        group.by = "Sample",
        label = FALSE)

#pdf("scaled_transformed_resolution_10kCells_featurePlots.pdf", height = 7, width = 10)

ptsize <- 0.3

FeaturePlot(seur.filtered, reduction = "umap", features = c("percent.mt"), pt.size = ptsize)

FeaturePlot(seur.filtered, reduction = "umap", features = c("ENSG00000181449.4"), pt.size = ptsize) + ggtitle("SOX2")
FeaturePlot(seur.filtered, reduction = "umap", features = c("ENSG00000148773.14"), pt.size = ptsize) + ggtitle("MKI67")
FeaturePlot(seur.filtered, reduction = "umap", features = c("ENSG00000007372.25"), pt.size = ptsize) + ggtitle("PAX6")

FeaturePlot(seur.filtered, reduction = "umap", features = c("ENSG00000136535.15"), pt.size = ptsize) + ggtitle("TBR1")
FeaturePlot(seur.filtered, reduction = "umap", features = c("ENSG00000163508.13"), pt.size = ptsize) + ggtitle("EOMES")

FeaturePlot(seur.filtered, reduction = "umap", features = c("ENSG00000171476.22"), pt.size = ptsize) + ggtitle("HOPX")
FeaturePlot(seur.filtered, reduction = "umap", features = c("ENSG00000134853.12"), pt.size = ptsize) + ggtitle("PDGFRA")

FeaturePlot(seur.filtered, reduction = "umap", features = c("ENSG00000077279.20"), pt.size = ptsize) + ggtitle("DCX")
FeaturePlot(seur.filtered, reduction = "umap", features = c("ENSG00000258947.8"), pt.size = ptsize) + ggtitle("TUBB3")
FeaturePlot(seur.filtered, reduction = "umap", features = c("ENSG00000164600.7"), pt.size = ptsize) + ggtitle("NEUROD6")

FeaturePlot(seur.filtered, reduction = "umap", features = c("ENSG00000127152.18"), pt.size = ptsize) + ggtitle("BCL11B")
FeaturePlot(seur.filtered, reduction = "umap", features = c("ENSG00000119042.17"), pt.size = ptsize) + ggtitle("SATB2")

FeaturePlot(seur.filtered, reduction = "umap", features = c("ENSG00000204531.21"), pt.size = ptsize) + ggtitle("OCT4")
FeaturePlot(seur.filtered, reduction = "umap", features = c("ENSG00000144355.15"), pt.size = ptsize) + ggtitle("DLX1")
FeaturePlot(seur.filtered, reduction = "umap", features = c("ENSG00000115844.11"), pt.size = ptsize) + ggtitle("DLX2")
FeaturePlot(seur.filtered, reduction = "umap", features = c("ENSG00000231764.11"), pt.size = ptsize) + ggtitle("DLX6-AS1")

FeaturePlot(seur.filtered, reduction = "umap", features = c("ENSG00000118271.12"), pt.size = ptsize) + ggtitle("TTR")
FeaturePlot(seur.filtered, reduction = "umap", features = c("ENSG00000176165.12"), pt.size = ptsize) + ggtitle("FOXG1")
FeaturePlot(seur.filtered, reduction = "umap", features = c("ENSG00000180613.11"), pt.size = ptsize) + ggtitle("GSX2")
FeaturePlot(seur.filtered, reduction = "umap", features = c("ENSG00000106852.16"), pt.size = ptsize) + ggtitle("LHX6")

#HTR1A, HTR1B, HTR1D, HTR2A, HTR2B, HTR2C, HTR3A, HTR3B, HTR5A, HTR6, HTR7, ADRB1, ADRB2, DRD1, DRD2, DRD3, DRD4, DRD5, HRH3, HRH4

# gene.names <- c("DCX", "BCL11B", "SATB2", "HTR1A", "HTR1B", "HTR1D", "HTR2A", "HTR2B", "HTR2C", "HTR3A", "HTR3B", "HTR5A", "HTR6", "HTR7",
#                 "ADRB1", "ADRB2", "DRD1", "DRD2", "DRD3", "DRD4", "DRD5", "HRH3", "HRH4")
#
# gene.ids <- df.gencode$gene_id[match(gene.names, df.gencode$gene_name)]
#
# DefaultAssay(seur.filtered) <- "SCT"
#
# plots <- FeaturePlot(seur.filtered, features = gene.ids, combine = FALSE, raster = TRUE)
# plots <- lapply(1:length(gene.ids), function(x) {plots[[x]] + labs(title = gene.names[x])})
#
# Reduce( `+`, plots ) +
#     patchwork::plot_layout( ncol = 4 )
#
# ggsave("UMAP.for.clozapine.PGP1.pdf", height = 20, width = 20)
#
# CombinePlots(plots, ncol = 4)
#
# print(plots[[1]])
#
# gene.to.plot.by.name <- c("SOX2", "MKI67", "PAX6", "ZBED1")



dev.off()


printMessage("Finished")

