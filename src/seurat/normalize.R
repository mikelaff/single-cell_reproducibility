# Normalize and doublet filter PGP1 samples


library(fishpond)
#library(data.table)
library(Seurat)
library(miQC)
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



# INPUT FILES ##########################################################################################################
# merged seurat filtered file
seurat.merged.filtered.rds <- here("results/seurat/20230126_PGP1_filtered_seurat_object.rds")

# gencode gtf file
#gencode.gtf <- here("data/refgenome/gencode/gencode.v40.annotation.gtf.gz")
# gencode gene names, types, and ensgid
df.gencode.csv <- here("data/refgenome/gencode/gencode.v40.genes.csv.gz")

# GLOBALS ##############################################################################################################
#

# Load GTF File #######
#gencode <- read_gff(gencode.gtf)

#df.gencode <- as_tibble(gencode)

#rm(gencode)

#df.gencode %<>%
#    filter(type == "gene") %>%
#    select(gene_id, gene_type, gene_name, chrom = seqnames)

df.gencode <- read_csv(df.gencode.csv)

# Load Seurat Object ########
seur.filtered <- readRDS(seurat.merged.filtered.rds)
dim(seur.filtered)
seur.filtered@meta.data$Sample <- paste(seur.filtered@meta.data$Site, seur.filtered@meta.data$Day, seur.filtered@meta.data$Rep, sep = "_")

# Remove Zelda
seur.filtered <- subset(seur.filtered, subset = Site %in% c("UNC", "CN", "CHOP"))

#str(seur.filtered)

# Subsample Cells
seur.sample <- seur.filtered[, sample(colnames(seur.filtered), size = 20000, replace = FALSE)]
dim(seur.sample)

#seur.sample@meta.data$Sample <- paste(seur.sample@meta.data$Site, seur.sample@meta.data$Day, seur.sample@meta.data$Rep, sep = "_")
rm(seur.filtered)
seur.filtered <- seur.sample
rm(seur.sample)
# Normalize
#seur.filtered <- NormalizeData(seur.filtered)
#str(seur.filtered)

# Find Variable Features
#seur.filtered <- FindVariableFeatures(seur.filtered)
#str(seur.filtered)

# Scale Data
#all.features <- rownames(seur.filtered@assays$RNA@counts)
#seur.filtered <- ScaleData(seur.filtered, features = all.features)

df.gencode %>%
    filter(chrom == "chrM") %>%
    pull(gene_id) -> ensg.mt

seur.filtered <- PercentageFeatureSet(object = seur.filtered, features = ensg.mt, col.name = "percent.mt")

# Run SCTransform
seur.filtered <- SCTransform(seur.filtered, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = TRUE)

#pdf("scaled_transformed_resolution_0.8_full.pdf", height = 7, width = 10)
seur.filtered <- RunPCA(seur.filtered, verbose = TRUE)
ElbowPlot(seur.filtered)

DimPlot(seur.filtered, reduction = "pca", group.by = "Day")
DimPlot(seur.filtered, reduction = "pca", group.by = "Site")
DimPlot(seur.filtered, reduction = "pca", group.by = "WTK_ID")
FeaturePlot(seur.filtered, reduction = "pca", features = c("percent.mt"))

seur.filtered <- RunUMAP(seur.filtered, dims = 1:10, verbose = TRUE)

seur.filtered <- FindNeighbors(seur.filtered, dims = 1:10, verbose = TRUE)
seur.filtered <- FindClusters(seur.filtered, verbose = TRUE, resolution = 0.8)

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
ggsave("iddrc_umap_by_site.pdf", height = 5, width = 7)
DimPlot(seur.filtered,
        group.by = "Sample",
        label = FALSE)
#FeaturePlot(seur.filtered, reduction = "umap", features = c("percent.mt"))

dev.off()

