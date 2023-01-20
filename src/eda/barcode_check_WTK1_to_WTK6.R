# check barcode mapping to samples after alevin-fry

library(SingleCellExperiment)
library(dplyr)
library(readr)
library(magrittr)
library(here)

# OUTPUT FILES #########################################################################################################

# INPUT FILES ##########################################################################################################
# Samples Data Table
df.samples.csv <- here("data/metadata/WTK1_to_WTK6_mapping_by_sample.csv")
# Barcodes Data Table
df.barcodes.csv <- here("data/metadata/WTK1_to_WTK6_mapping_by_barcode.csv")

# alevin-fry quantification root directory
al.fry.root.dir <- "/proj/steinlab/projects/IVIV_scRNA/alevin_fry_WTK1_to_WTK6/"

# GLOBALS ##############################################################################################################

# Import alevin-fry ############

file.dir <- list.files(path = al.fry.root.dir)

wtk6.dirs <- file.dir[grep("WTK6", file.dir)][2]

sce.wtk6.1 <- loadFry(fryDir = paste0(al.fry.root.dir, wtk6.dirs[1], "/", wtk6.dirs[1], "_alevin", "/", "countMatrix"), outputFormat = "snRNA")

View(as.data.frame(colData(sce.wtk6.1)))

sce.wtk6.1$bc1 <- gsub("([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})", "\\3", sce.wtk6.1$barcodes,  perl=T)
sce.wtk6.1$bc2 <- gsub("([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})", "\\2", sce.wtk6.1$barcodes,  perl=T)
sce.wtk6.1$bc3 <- gsub("([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})", "\\1", sce.wtk6.1$barcodes,  perl=T)

