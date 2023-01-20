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
#

# INPUT FILES ##########################################################################################################
# Samples Data Table
df.samples.csv <- here("data/metadata/WTK1_to_WTK6_mapping_by_sample.csv")
# Barcodes Data Table
df.barcodes.csv <- here("data/metadata/WTK1_to_WTK6_mapping_by_barcode.csv")
# Sublibraries Data Table
df.sublibraries.csv <- here("data/metadata/WTK1_to_WTK6_sublibraries.csv")

# alevin-fry quantification root directory
al.fry.root.dir <- "/proj/steinlab/projects/IVIV_scRNA/alevin_fry_WTK1_to_WTK6/"

# GLOBALS ##############################################################################################################
#

# Import Metadata ######
df.barcodes <- read_csv(df.barcodes.csv)
df.samples <- read_csv(df.samples.csv)

df.sublibraries <- read_csv(df.sublibraries.csv)



# Import alevin-fry ############

df.sublibraries$count_Cell_Barcodes <- NA
df.sublibraries$count_Genes <- NA

df.count.data <- tibble()

# loop over alevin-fry count matrix directories
for (i in 1:nrow(df.sublibraries)) {

    printMessage()
    printMessage(paste("Working on sublibrary", i, "of", nrow(df.sublibraries), ":", df.sublibraries$Sublibrary_ID[i]))
    printMessage()

    sce <- NULL
    dir.counts <- NULL
    df.coldata <- NULL
    sublibrary.samples <- NULL

    dir.counts <- paste0(al.fry.root.dir, df.sublibraries$Sublibrary_ID[i], "/", df.sublibraries$Sublibrary_ID[i], "_alevin", "/", "countMatrix")

    # load alevin-fry count matrix as single cell experiment using "snRNA" which is U+S+A
    sce <- loadFry(fryDir = dir.counts, outputFormat = "snRNA")

    # sublibrary stats
    # df.sublibraries$count_Cell_Barcodes[i] <- dim(sce)[2]
    # df.sublibraries$count_Genes[i] <- dim(sce)[1]

    # cell barcode stats
    df.coldata <- as_tibble(colData(sce))

    df.coldata$bc1 <- gsub("([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})", "\\3", df.coldata$barcodes,  perl=T)
    df.coldata$bc2 <- gsub("([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})", "\\2", df.coldata$barcodes,  perl=T)
    df.coldata$bc3 <- gsub("([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})", "\\1", df.coldata$barcodes,  perl=T)

    df.coldata %<>%
        left_join(dplyr::filter(df.barcodes, WTK_ID == df.sublibraries$WTK_ID[i], type == "T"), by = c("bc1" = "sequence"))

    # check that all cell barcodes get assigned to a sample
    if (sum(is.na(df.coldata$CODissoID))) {
        printMessage("ERROR: Some cell barcodes not assigned to samples!", fillChar = "&")
    }

    # check that all samples are represented in this sublibrary
    df.samples %>%
        filter(WTK_ID == df.sublibraries$WTK_ID[i]) %>%
        pull(CODissoID) -> sublibrary.samples
    if (!all(sublibrary.samples %in% df.coldata$CODissoID)) {
        printMessage("ERROR: Some samples not represented in this sublibrary!", fillChar = "&")
    }


    df.coldata %>%
        group_by(CODissoID, WTK_ID, well, bc1) %>%
        summarise(count_cell_barcodes = n()) %>%
        mutate(count_genes = dim(sce)[1],
               Sublibrary_ID = df.sublibraries$Sublibrary_ID[i]) -> df.tmp

    df.count.data %<>%
        bind_rows(df.tmp)

}

