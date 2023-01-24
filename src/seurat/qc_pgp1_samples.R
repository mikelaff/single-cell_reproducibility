# QC pgp1 samples


library(fishpond)
#library(data.table)
library(Seurat)
library(miQC)
library(SeuratWrappers)
library(flexmix)
library(SingleCellExperiment)
library(Matrix)
library(stringr)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(here)

library(plyranges)

library(mikelaffr)

# OUTPUT FILES #########################################################################################################

# merged seurat filtered file
seurat.merged.filtered.rds <- here("results/seurat/20230123_PGP1_filtered_seurat_object.rds")

# INPUT FILES ##########################################################################################################
# merged seurat file, unfiltered
seurat.merged.rds <- here("results/seurat/PGP1_seurat_object.rds")

# gencode gtf file
gencode.gtf <- here("data/refgenome/gencode/gencode.v40.annotation.gtf.gz")

# GLOBALS ##############################################################################################################
#



# Load GTF File #######
gencode <- read_gff(gencode.gtf)

df.gencode <- as_tibble(gencode)

rm(gencode)

df.gencode %<>%
    filter(type == "gene") %>%
    select(gene_id, gene_type, gene_name, chrom = seqnames)


# Load Seurat Object ########
seurat.merged <- readRDS(seurat.merged.rds)


# percent mt reads
df.gencode %>%
    filter(chrom == "chrM") %>%
    pull(gene_id) -> ensg.mt

seurat.merged <- PercentageFeatureSet(object = seurat.merged, features = ensg.mt, col.name = "percent.mt")

df.cellData <- as_tibble(seurat.merged@meta.data)

df.cellData %>%
    ggplot(aes(x = nCount_RNA, y = nFeature_RNA)) +
    geom_point()

df.cellData %>%
    mutate(SampleName = paste(Site, Day, Rep, sep = "_")) %>%
    group_by(SampleName, WTK_ID) %>%
    dplyr::summarise(cell_count = dplyr::n(),
                     num_wells = n_distinct(well)) %>%
    ggplot(aes(x = SampleName, y = cell_count)) +
    geom_col()

df.cellData %>%
    mutate(SampleName = paste(Site, Day, Rep, sep = "_")) %>%
    group_by(SampleName, WTK_ID) %>%
    dplyr::summarise(cell_count = dplyr::n(),
                     num_wells = n_distinct(well)) %>%
    ggplot(aes(x = num_wells, y = cell_count)) +
    geom_point()

df.cellData %>%
    ggplot(aes(y = percent.mt, x = orig.ident)) +
    geom_violin() +
    geom_jitter()

df.cellData %>%
    ggplot(aes(x = nCount_RNA, y = percent.mt)) +
    geom_hex()

df.cellData %>%
    ggplot(aes(x = nCount_RNA)) +
    geom_density()

df.cellData %>%
    ggplot(aes(x = nFeature_RNA)) +
    geom_density()

sum(df.cellData$nCount_RNA < 100)
sum(df.cellData$nFeature_RNA < 100)

sum(df.cellData$nFeature_RNA < 100 | df.cellData$nCount_RNA < 100)

seurat.merged.filtered <- subset(seurat.merged, subset = nCount_RNA > 100 & nFeature_RNA > 100)

df.cellData.filtered <- as_tibble(seurat.merged.filtered@meta.data)

df.cellData.filtered %>%
    mutate(SampleName = paste(Site, Day, Rep, sep = "_")) %>%
    group_by(SampleName, WTK_ID) %>%
    dplyr::summarise(cell_count = dplyr::n(),
                     num_wells = n_distinct(well)) %>%
    ggplot(aes(x = SampleName, y = cell_count)) +
    geom_col()

seurat.merged.filtered <- RunMiQC(seurat.merged.filtered, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA",
                                  posterior.cutoff = 0.75, model.slot = "flexmix_model")

seurat.merged.filtered.again <- subset(seurat.merged.filtered, miQC.keep == "keep")
dim(seurat.merged.filtered.again)

#VlnPlot(seurat.merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)



