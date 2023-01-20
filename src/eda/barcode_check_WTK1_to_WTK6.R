# check barcode mapping to samples after alevin-fry

library(fishpond)
library(SingleCellExperiment)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(here)

library(mikelaffr)

# OUTPUT FILES #########################################################################################################
# summarized barcode/sample count data
df.count.data.output.rds <- here("results/alevin_fry/WTK1_to_WTK6_summarized_count_data.rds")

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

printMessage("Saving count data...")
saveRDS(df.count.data, df.count.data.output.rds)

# Load Count Summary Data ###########
df.count.data
df.count.data %>%
    group_by(Sublibrary_ID, WTK_ID) %>%
    summarise(cells_per_sublibrary = sum(count_cell_barcodes)) %>%
    ggplot(aes(x = Sublibrary_ID, y = cells_per_sublibrary, fill = WTK_ID)) +
    geom_col()



# file.dir <- list.files(path = al.fry.root.dir)
#
# wtk6.dirs <- file.dir[grep("WTK6", file.dir)]
#
# sce.wtk6.1 <- loadFry(fryDir = paste0(al.fry.root.dir, wtk6.dirs[1], "/", wtk6.dirs[1], "_alevin", "/", "countMatrix"), outputFormat = "snRNA")
#
#
# sce.wtk6.1$bc1 <- gsub("([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})", "\\3", sce.wtk6.1$barcodes,  perl=T)
# sce.wtk6.1$bc2 <- gsub("([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})", "\\2", sce.wtk6.1$barcodes,  perl=T)
# sce.wtk6.1$bc3 <- gsub("([ACTGN]{8})([ACTGN]{8})([ACTGN]{8})", "\\1", sce.wtk6.1$barcodes,  perl=T)
#
# df.coldata <- as_tibble(colData(sce.wtk6.1))
#
# # Import Barcode List #######
#
# df.barcodes %<>%
#     filter(WTK_ID == "WTK6")
#
# all(df.barcodes$sequence[df.barcodes$type == "T"] %in% df.coldata$bc1)
# all(df.coldata$bc1 %in% df.barcodes$sequence[df.barcodes$type == "T"])
#
# df.coldata %<>%
#     left_join(df.barcodes, by = c("bc1" = "sequence"))
#
# summary(as.factor(df.coldata$type))
# sum(is.na(df.coldata$CODissoID))
#
# df.coldata %>%
#     group_by(CODissoID) %>%
#     summarise(count_cells = n()) %>%
#     ggplot(aes(x = CODissoID, y = count_cells)) +
#     geom_col()
#
# df.coldata %>%
#     group_by(bc1) %>%
#     summarise(count_cells = n())
