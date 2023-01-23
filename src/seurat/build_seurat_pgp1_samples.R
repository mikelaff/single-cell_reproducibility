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
sce.merged.rds <- here("results/alevin_fry/WTK5_to_WTK6_merged_sce.rds")

# merged seurat file
seurat.merged.rds <- here("results/seurat/PGP1_seurat_object.rds")

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
df.count.data.rds <- here("results/alevin_fry/WTK1_to_WTK6_summarized_count_data.rds")

# GLOBALS ##############################################################################################################
#


# Import Metadata ######
df.barcodes <- read_csv(df.barcodes.csv)
df.samples <- read_csv(df.samples.csv)

df.sublibraries <- read_csv(df.sublibraries.csv)
df.sublibraries %<>%
    mutate(Sublibrary = paste("SL", Sublibrary, sep = ""))

df.pgp1.samples <- read_csv(df.repro.samples.csv)

df.pgp1.samples %<>%
    left_join(df.samples, by = "CODissoID")

df.pgp1.samples %<>%
    select(-randHEX_barcodes)

# select only PGP1 samples
df.samples %>%
    filter(WTK_ID == "WTK5" | WTK_ID == "WTK6",
           grepl("UNC", CODissoID) | grepl("CN", CODissoID) | grepl("CHOP", CODissoID) | grepl("Zelda", CODissoID)) %>%
    pull(CODissoID) -> samples.pgp1

sum(duplicated(samples.pgp1))

# Load Count Summary Data ###########
# df.count.data <- readRDS(df.count.data.rds)
#
# # filter for PGP1 samples
# df.count.data %>%
#     filter(CODissoID %in% samples.pgp1) -> df.count.data.pgp1
#
# df.count.data.pgp1 %>%
#     group_by(CODissoID) %>%
#     summarise(cells_per_sample = sum(count_cell_barcodes)) %>%
#     ggplot(aes(x = CODissoID, y = cells_per_sample)) +
#     geom_col() +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14),
#           axis.text.y = element_text(size = 14),
#           axis.title.y = element_text(size = 14),
#           title = element_text(size = 18)) +
#     labs(y = "Cells per Sample\n(raw unfiltered)",
#          title = "Alevin-Fry Quantified Single Cells") +
#     scale_fill_manual(values = cbPalette) +
#     scale_y_continuous(labels = scales::label_comma())
#
# ggsave("alevin-fry_quant_raw_unfiltered_pgp1_samples.pdf", width = 14, height = 8)

# Import alevin-fry ############

# # import only WTK5 and WTK6 sublibraries
# df.sublibraries %>%
#     filter(WTK_ID == "WTK5" | WTK_ID == "WTK6") -> df.sublibraries.pgp1
#
# # loop over alevin-fry count matrix directories
# for (i in 1:nrow(df.sublibraries.pgp1)) {
#
#     printMessage()
#     printMessage(paste("Working on sublibrary", i, "of", nrow(df.sublibraries.pgp1), ":", df.sublibraries.pgp1$Sublibrary_ID[i]))
#     printMessage()
#
#     sce <- NULL
#     dir.counts <- NULL
#     df.coldata <- NULL
#     #seurat.obj <- NULL
#
#     dir.counts <- paste0(al.fry.root.dir, df.sublibraries.pgp1$Sublibrary_ID[i], "/", df.sublibraries.pgp1$Sublibrary_ID[i], "_alevin", "/", "countMatrix")
#
#     # load alevin-fry count matrix as single cell experiment using "snRNA" which is U+S+A
#     sce <- loadFry(fryDir = dir.counts, outputFormat = "snRNA")
#
#
#     df.coldata <- as_tibble(colData(sce))
#
#
#     df.coldata$bc1 <- gsub("([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})", "\\3", df.coldata$barcodes,  perl=T)
#     df.coldata$bc2 <- gsub("([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})", "\\2", df.coldata$barcodes,  perl=T)
#     df.coldata$bc3 <- gsub("([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})", "\\1", df.coldata$barcodes,  perl=T)
#
#     df.coldata %<>%
#         left_join(dplyr::filter(df.barcodes, WTK_ID == df.sublibraries.pgp1$WTK_ID[i], type == "T"), by = c("bc1" = "sequence"))
#
#     df.coldata %<>%
#         mutate(Sublibrary_ID = df.sublibraries.pgp1$Sublibrary_ID[i])
#
#     df.coldata %<>%
#         mutate(Cell_Barcode_ID = paste(barcodes, Sublibrary_ID, sep = "_"))
#
#     df.coldata <- as.data.frame(df.coldata)
#     rownames(df.coldata) <- df.coldata$Cell_Barcode_ID
#
#     colnames(sce) <- df.coldata$Cell_Barcode_ID
#
#     # seurat.obj <- CreateSeuratObject(assay(sce), project = df.sublibraries.pgp1$Sublibrary_ID[i])
#     # seurat.obj <- AddMetaData(seurat.obj, df.coldata)
#     #
#     # seurat.obj <- subset(seurat.obj, subset = nCount_RNA > 1500 & nFeature_RNA > 1000)
#
#     if (i == 1) {
#         sce.merged <- sce
#     } else {
#         sce.merged <- cbind(sce.merged, sce)
#     }
#
# }
#
# saveRDS(sce.merged, sce.merged.rds)

# Load alevin-fry ################
sce.merged <- readRDS(sce.merged.rds)

# build out cell metadata and link with samples data
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
    left_join(dplyr::select(df.pgp1.samples, CODissoID, Site, Day, Rep, WTK_ID), by = c("WTK_ID", "CODissoID"))


# select only samples we want: pgp1
sce.merged <- sce.merged[, df.cellData$CODissoID %in% df.pgp1.samples$CODissoID]
df.cellData <- df.cellData[df.cellData$CODissoID %in% df.pgp1.samples$CODissoID, ]

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

seurat.merged <- CreateSeuratObject(counts = counts(sce.merged), project = "PGP1_Samples", meta.data = df.cellData)

# save seurat object
saveRDS(seurat.merged, seurat.merged.rds)

# seurat.merged@meta.data
#
# str(seurat.merged)
#
# seurat.merged@assays$RNA@meta.features

