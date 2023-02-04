# PGP1 Samples from WTK5 and WTK6 for reproducibility analysis

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

library(mikelaffr)

# OUTPUT FILES #########################################################################################################
# merged sce
sce.merged.rds <- here("results/alevin_fry/20230201_PGP1_merged_sce.rds")

# merged seurat file
seurat.merged.rds <- here("results/seurat/20230201_PGP1_seurat_object.rds")

# INPUT FILES ##########################################################################################################
# Samples Data Table
df.samples.csv <- here("data/metadata/WTK1_to_WTK6_mapping_by_sample.csv")
# Barcodes Data Table
df.barcodes.csv <- here("data/metadata/WTK1_to_WTK6_mapping_by_barcode.csv")
# Sublibraries Data Table
df.sublibraries.csv <- here("data/metadata/WTK1_to_WTK6_sublibraries.csv")

# Single-cell reproducibility samples
df.repro.samples.csv <- here("data/metadata/single_cell_reproducibility_samples.csv")

# alevin-fry quantification root directory
al.fry.root.dir <- "/proj/steinlab/projects/IVIV_scRNA/alevin_fry_WTK1_to_WTK6/"

# summarized barcode/sample count data
#df.count.data.rds <- here("results/alevin_fry/WTK1_to_WTK6_summarized_count_data.rds")

# GLOBALS ##############################################################################################################
#


# Import Metadata ######
df.barcodes <- read_csv(df.barcodes.csv)
df.samples <- read_csv(df.samples.csv)

df.sublibraries <- read_csv(df.sublibraries.csv)
df.sublibraries %<>%
    mutate(Sublibrary = paste("SL", Sublibrary, sep = ""))

df.samples.pgp1 <- read_csv(df.repro.samples.csv)

df.samples.pgp1 %<>%
    left_join(df.samples, by = c("CODissoID", "WTK_ID"))

df.samples.pgp1 %<>%
    select(-randHEX_barcodes)

sum(duplicated(df.samples.pgp1$CODissoID))


# Import alevin-fry ############

# import only WTK5 and WTK6 sublibraries
df.sublibraries %>%
    filter(WTK_ID %in% df.samples.pgp1$WTK_ID) -> df.sublibraries.pgp1

printMessage(fillChar = "-")
printMessage("Loading alevin-fry", fillChar = "@")
printMessage(fillChar = "-")

# loop over alevin-fry count matrix directories
for (i in 1:nrow(df.sublibraries.pgp1)) {

    printMessage()
    printMessage(paste("Working on sublibrary", i, "of", nrow(df.sublibraries.pgp1), ":", df.sublibraries.pgp1$Sublibrary_ID[i]))
    printMessage()

    sce <- NULL
    dir.counts <- NULL
    df.coldata <- NULL
    #seurat.obj <- NULL

    dir.counts <- paste0(al.fry.root.dir, df.sublibraries.pgp1$Sublibrary_ID[i], "/", df.sublibraries.pgp1$Sublibrary_ID[i], "_alevin", "/", "countMatrix")

    # load alevin-fry count matrix as single cell experiment using "snRNA" which is U+S+A
    sce <- loadFry(fryDir = dir.counts, outputFormat = "snRNA")


    df.coldata <- as_tibble(colData(sce))


    df.coldata$bc1 <- gsub("([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})", "\\3", df.coldata$barcodes,  perl=T)
    df.coldata$bc2 <- gsub("([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})", "\\2", df.coldata$barcodes,  perl=T)
    df.coldata$bc3 <- gsub("([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})", "\\1", df.coldata$barcodes,  perl=T)

    df.coldata %<>%
        left_join(dplyr::filter(df.barcodes, WTK_ID == df.sublibraries.pgp1$WTK_ID[i], type == "T"), by = c("bc1" = "sequence"))

    df.coldata %<>%
        mutate(Sublibrary_ID = df.sublibraries.pgp1$Sublibrary_ID[i])

    # create unique cell barcode ID
    df.coldata %<>%
        mutate(Cell_Barcode_ID = paste(barcodes, Sublibrary_ID, sep = "_"))

    #df.coldata <- as.data.frame(df.coldata)
    #rownames(df.coldata) <- df.coldata$Cell_Barcode_ID

    # relabel cells with unique name
    colnames(sce) <- df.coldata$Cell_Barcode_ID

    # select only samples we want: pgp1
    sce <- sce[, df.coldata$CODissoID %in% df.samples.pgp1$CODissoID]
    #df.coldata <- df.coldata[df.coldata$CODissoID %in% df.samples.pgp1$CODissoID, ]

    # merge sce
    if (i == 1) {
        sce.merged <- sce
    } else {
        sce.merged <- cbind(sce.merged, sce)
    }

}

printMessage(fillChar = "-")
printMessage("Finished Loading alevin-fry", fillChar = "@")
printMessage(fillChar = "-")

printMessage("Saving SCE Object")

saveRDS(sce.merged, sce.merged.rds)

# Load alevin-fry ################
#sce.merged <- readRDS(sce.merged.rds)

# build out cell metadata and link with samples data
printMessage("Building out cell metadata", fillChar = " ", justify = "l")
df.cellData <- as.data.frame(colData(sce.merged))

df.cellData$cell_label <- rownames(df.cellData)
# df.cellData$Sublibrary_ID <- paste(sapply(strsplit(df.cellData$cell_label, "_"), `[`, 2),
#                                    sapply(strsplit(df.cellData$cell_label, "_"), `[`, 3),
#                                    sapply(strsplit(df.cellData$cell_label, "_"), `[`, 4),
#                                    sapply(strsplit(df.cellData$cell_label, "_"), `[`, 5),
#                                    sapply(strsplit(df.cellData$cell_label, "_"), `[`, 6), sep = "_")

df.cellData$Sublibrary_ID <- sapply(regmatches(df.cellData$cell_label, regexpr("_", df.cellData$cell_label), invert = TRUE), `[`, 2)

df.cellData <- as_tibble(df.cellData)
# split round barcodes
df.cellData %<>%
    mutate(bc1 = gsub("([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})", "\\3", df.cellData$barcodes,  perl=T),
           bc2 = gsub("([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})", "\\2", df.cellData$barcodes,  perl=T),
           bc3 = gsub("([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})", "\\1", df.cellData$barcodes,  perl=T))

# get WTK_ID for barcode matching
df.cellData %<>%
    left_join(dplyr::select(df.sublibraries, Sublibrary_ID, WTK_ID, Sublibrary), by = "Sublibrary_ID")

# match first round barcodes and WTK_ID to get sample IDs
df.cellData %<>%
    left_join(dplyr::select(df.barcodes, WTK_ID, bc1 = sequence, well, CODissoID), by = c("WTK_ID","bc1"))

# add in sample metadata
df.cellData %<>%
    left_join(dplyr::select(df.samples.pgp1, CODissoID, Site, Day, Rep, WTK_ID), by = c("WTK_ID", "CODissoID"))


# select only samples we want: pgp1
sce.merged <- sce.merged[, df.cellData$CODissoID %in% df.samples.pgp1$CODissoID]
df.cellData <- df.cellData[df.cellData$CODissoID %in% df.samples.pgp1$CODissoID, ]

stopifnot(all(colnames(sce.merged) == df.cellData$cell_label))

# create easy to read cell name
df.cellData %<>%
    mutate(cell_name = paste(barcodes, WTK_ID, Sublibrary, Site, Day, Rep, sep = "_"))

sum(duplicated(df.cellData$cell_name))

df.cellData %<>%
    dplyr::select(-cell_label)

df.cellData <- as.data.frame(df.cellData)

# label rows, label cells in sce
rownames(df.cellData) <- df.cellData$cell_name
colnames(sce.merged) <- rownames(df.cellData)

stopifnot(all(colnames(sce.merged) == rownames(df.cellData)))

# Build Seurat Object #################
printMessage("Building Seurat Object")
seurat.merged <- CreateSeuratObject(counts = counts(sce.merged), project = "PGP1_Samples", meta.data = df.cellData)

# save seurat object
printMessage("Saving Seurat Object")
saveRDS(seurat.merged, seurat.merged.rds)

printMessage("Finished")
# seurat.merged@meta.data
#
# str(seurat.merged)
#
# seurat.merged@assays$RNA@meta.features

