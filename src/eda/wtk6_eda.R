# wtk6 analysis in seurat

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
library(magrittr)
library(ggplot2)
library(reshape2)

# OUTPUT FILES #########################################################################################################
#

# INPUT FILES ##########################################################################################################
# alevin-fry quantification root directory
al.fry.root.dir <- "/proj/steinlab/projects/IVIV_scRNA/alevin_fry_WTK1_to_WTK6/"

# GLOBALS ##############################################################################################################
#

# Import Quant Data ###############

file.dir <- list.files(path = al.fry.root.dir, )

wtk6.dirs <- file.dir[grep("WTK6", file.dir)]

sce.wtk6.1 <- loadFry(fryDir = paste0(al.fry.root.dir, wtk6.dirs[1], "/", wtk6.dirs[1], "_alevin", "/", "countMatrix"), outputFormat = "raw")
sce.wtk6.2 <- loadFry(fryDir = paste0(al.fry.root.dir, wtk6.dirs[2], "/", wtk6.dirs[2], "_alevin", "/", "countMatrix"), outputFormat = "raw")
sce.wtk6.3 <- loadFry(fryDir = paste0(al.fry.root.dir, wtk6.dirs[3], "/", wtk6.dirs[3], "_alevin", "/", "countMatrix"), outputFormat = "raw")
sce.wtk6.4 <- loadFry(fryDir = paste0(al.fry.root.dir, wtk6.dirs[4], "/", wtk6.dirs[4], "_alevin", "/", "countMatrix"), outputFormat = "raw")
sce.wtk6.5 <- loadFry(fryDir = paste0(al.fry.root.dir, wtk6.dirs[5], "/", wtk6.dirs[5], "_alevin", "/", "countMatrix"), outputFormat = "raw")
sce.wtk6.6 <- loadFry(fryDir = paste0(al.fry.root.dir, wtk6.dirs[6], "/", wtk6.dirs[6], "_alevin", "/", "countMatrix"), outputFormat = "raw")
sce.wtk6.7 <- loadFry(fryDir = paste0(al.fry.root.dir, wtk6.dirs[7], "/", wtk6.dirs[7], "_alevin", "/", "countMatrix"), outputFormat = "raw")
sce.wtk6.8 <- loadFry(fryDir = paste0(al.fry.root.dir, wtk6.dirs[8], "/", wtk6.dirs[8], "_alevin", "/", "countMatrix"), outputFormat = "raw")


df.genes.wtk6 <- tibble(ensg = rownames(assay(sce.wtk6.1, i = "spliced")),
                          spliced_counts = rowSums(assay(sce.wtk6.1, i = "spliced")) +
                                            rowSums(assay(sce.wtk6.2, i = "spliced")) +
                                            rowSums(assay(sce.wtk6.3, i = "spliced")) +
                                            rowSums(assay(sce.wtk6.4, i = "spliced")) +
                                            rowSums(assay(sce.wtk6.5, i = "spliced")) +
                                            rowSums(assay(sce.wtk6.6, i = "spliced")) +
                                            rowSums(assay(sce.wtk6.7, i = "spliced")) +
                                            rowSums(assay(sce.wtk6.8, i = "spliced")),
                          unspliced_counts = rowSums(assay(sce.wtk6.1, i = "unspliced")) +
                                          rowSums(assay(sce.wtk6.2, i = "unspliced")) +
                                          rowSums(assay(sce.wtk6.3, i = "unspliced")) +
                                          rowSums(assay(sce.wtk6.4, i = "unspliced")) +
                                          rowSums(assay(sce.wtk6.5, i = "unspliced")) +
                                          rowSums(assay(sce.wtk6.6, i = "unspliced")) +
                                          rowSums(assay(sce.wtk6.7, i = "unspliced")) +
                                          rowSums(assay(sce.wtk6.8, i = "unspliced")),
                          ambiguous_counts = rowSums(assay(sce.wtk6.1, i = "ambiguous")) +
                              rowSums(assay(sce.wtk6.2, i = "ambiguous")) +
                              rowSums(assay(sce.wtk6.3, i = "ambiguous")) +
                              rowSums(assay(sce.wtk6.4, i = "ambiguous")) +
                              rowSums(assay(sce.wtk6.5, i = "ambiguous")) +
                              rowSums(assay(sce.wtk6.6, i = "ambiguous")) +
                              rowSums(assay(sce.wtk6.7, i = "ambiguous")) +
                              rowSums(assay(sce.wtk6.8, i = "ambiguous")))

df.genes.wtk6 %>%
    filter(spliced_counts + unspliced_counts + ambiguous_counts > 0)

df.genes.wtk6 %<>%
    filter(spliced_counts + unspliced_counts + ambiguous_counts > 0) %>%
    mutate(total_counts = spliced_counts + unspliced_counts + ambiguous_counts) %>%
    mutate(spliced_over_total = spliced_counts / total_counts,
           unspliced_over_total = unspliced_counts / total_counts,
           ambiguous_over_total = ambiguous_counts / total_counts)

df.genes.wtk6 %<>%
    mutate(spliced_over_unspliced = spliced_counts / unspliced_counts)

#pdf()
df.genes.wtk6 %>%
    ggplot(aes(x = log10(total_counts))) +
    geom_histogram()

df.genes.wtk6 %>%
    ggplot(aes(x = total_counts)) +
    geom_histogram()

df.genes.wtk6 %>%
    ggplot(aes(x = log10(spliced_counts), y = log10(unspliced_counts))) +
    geom_point() +
    geom_smooth(method = "lm")

df.genes.wtk6 %>%
    ggplot(aes(x = spliced_over_total, y = unspliced_over_total)) +
    geom_point(size = 0.5)

df.genes.wtk6 %>%
    ggplot(aes(x = spliced_over_total, y = unspliced_over_total)) +
    stat_bin_hex(bins = 20)

df.genes.wtk6 %>%
    ggplot(aes(x = spliced_over_total, y = ambiguous_over_total)) +
    geom_point(size =.5)

df.genes.wtk6 %>%
    ggplot(aes(x = spliced_over_total, y = ambiguous_over_total)) +
    stat_bin_hex(bins = 20)

df.genes.wtk6 %>%
    ggplot(aes(x = spliced_counts + ambiguous_counts, y = unspliced_counts)) +
    geom_point(size =.5)

df.genes.wtk6 %>%
    ggplot(aes(x = spliced_counts + ambiguous_counts, y = unspliced_counts)) +
    stat_bin_hex()

df.genes.wtk6 %>%
    ggplot(aes(x = spliced_over_unspliced, y = total_counts)) +
    geom_point()

#dev.off()

df.genes.wtk6.melted <- as_tibble(melt(df.genes.wtk6))

df.genes.wtk6.melted %>%
    ggplot(aes(x = value, fill = variable)) +
    geom_histogram(position = "dodge")


all(rownames(sce.wtk6.1) == rownames(sce.wtk6.2))

sum(colnames(sce.wtk6.1) %in% colnames(sce.wtk6.2))
sum(colnames(sce.wtk6.1) %in% colnames(sce.wtk6.3))

sum(duplicated(colnames(sce.wtk6.1)))

allbarcodes <- c(colnames(sce.wtk6.1), colnames(sce.wtk6.2), colnames(sce.wtk6.3), colnames(sce.wtk6.4),
                 colnames(sce.wtk6.5), colnames(sce.wtk6.6), colnames(sce.wtk6.7), colnames(sce.wtk6.8))

sum(duplicated(allbarcodes))

seurat.wtk6.1 <- CreateSeuratObject(assay(sce.wtk6.1, i = "spliced"), project = "wtk6_sl1")
seurat.wtk6.2 <- CreateSeuratObject(assay(sce.wtk6.2, i = "spliced"), project = "wtk6_sl2")

seurat.wtk6.1.2 <- merge(x = seurat.wtk6.1, y = seurat.wtk6.2)

assay(sce, i = "spliced")[10:50,10:50]
assay(sce, i = "unspliced")[10:50,10:50]
assay(sce, i = "ambiguous")[10:50,10:50]

df.genes <- tibble(ensg = rownames(assay(sce, i = "spliced")),
                   spliced_counts = rowSums(assay(sce, i = "spliced")),
                   unspliced_counts = rowSums(assay(sce, i = "unspliced")),
                   ambiguous_counts = rowSums(assay(sce, i = "ambiguous")))

rowSums(assay(sce, i = "spliced"))
sum(assay(sce, i = "unspliced"))
sum(assay(sce, i = "ambiguous"))


wtk6.sl7.seurat <- CreateSeuratObject()

# labeling cells

sce.wtk6.1 <- loadFry(fryDir = paste0(al.fry.root.dir, wtk6.dirs[1], "/", wtk6.dirs[1], "_alevin", "/", "countMatrix"), outputFormat = "snRNA")

View(as.data.frame(colData(sce.wtk6.1)))

sce.wtk6.1$bc1 <- gsub("([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})", "\\3", sce.wtk6.1$barcodes,  perl=T)
sce.wtk6.1$bc2 <- gsub("([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})", "\\2", sce.wtk6.1$barcodes,  perl=T)
sce.wtk6.1$bc3 <- gsub("([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})", "\\1", sce.wtk6.1$barcodes,  perl=T)

seurat.wtk6.1 <- CreateSeuratObject(assay(sce.wtk6.1), project = "wtk6_sl1")


str(seurat.wtk6.1)

head(seurat.wtk6.1$nFeature_RNA)
head(seurat.wtk6.1$nCount_RNA)

all(colnames(seurat.wtk6.1) == rownames(colData(sce.wtk6.1)))
seurat.wtk6.1 <- AddMetaData(seurat.wtk6.1, as.data.frame(colData(sce.wtk6.1)))

seurat.wtk6.1.filt <- subset(seurat.wtk6.1, subset = nCount_RNA > 1500 & nFeature_RNA > 1000)

View(seurat.wtk6.1@meta.data)

#df.tmp <- readRDS("/proj/steinlab/projects/IVIV_scRNA/InfoVerifiedID/RObjects/WTK1Verified.rds")
