# curate WTK metadata

library(dplyr)
library(readr)
library(magrittr)
library(here)

# OUTPUT FILES #########################################################################################################
# Samples Data Table
df.samples.csv <- here("data/metadata/WTK1_to_WTK6_mapping_by_sample.csv")
# Barcodes Data Table
df.barcodes.csv <- here("data/metadata/WTK1_to_WTK6_mapping_by_barcode.csv")
# Sublibraries Data Table
df.sublibraries.csv <- here("data/metadata/WTK1_to_WTK6_sublibraries.csv")

# INPUT FILES ##########################################################################################################
# WTK1 Mapping csv
wtk1.mapping.csv <- "/proj/steinlab/projects/IVIV_scRNA/InfoVerifiedID/WTK1Mapping.csv"
# WTK2 Mapping csv
wtk2.mapping.csv <- "/proj/steinlab/projects/IVIV_scRNA/InfoVerifiedID/WTK2Mapping.csv"
# WTK3 Mapping csv
wtk3.mapping.csv <- "/proj/steinlab/projects/IVIV_scRNA/InfoVerifiedID/WTK3Mapping.csv"
# WTK4 Mapping csv
wtk4.mapping.csv <- "/proj/steinlab/projects/IVIV_scRNA/InfoVerifiedID/WTK4Mapping.csv"
# WTK5 Mapping csv
wtk5.mapping.csv <- "/proj/steinlab/projects/IVIV_scRNA/InfoVerifiedID/WTK5Mapping.csv"
# WTK6 Mapping csv
wtk6.mapping.csv <- "/proj/steinlab/projects/IVIV_scRNA/InfoVerifiedID/WTK6Mapping.csv"

# First barcodes (48 oligo dT barcodes, 48 rand hex barcodes)
barcode.one.csv <- here("data/barcode/bc_data_v2.csv")

# Sublibrary directories root
sublibrary.dir <- "/proj/steinlab/projects/IVIV_scRNA/alevin_fry_WTK1_to_WTK6/"

# GLOBALS ##############################################################################################################
# Well A1 to D12
WELLS <- c(paste0("A", seq(1,12)), paste0("B", seq(1,12)), paste0("C", seq(1,12)), paste0("D", seq(1,12)))

# Import Well Mapping Files #############
df.samples <- tibble()

# WTK1
df.wtk1.mapping <- read_csv(wtk1.mapping.csv)

df.wtk1.mapping %<>%
    mutate(WTK_ID = "WTK1") %>%
    select(WTK_ID,
           SampleNumber,
           CODissoID,
           Wells)

# remove spaces between dashes for consistency
df.wtk1.mapping %<>%
    mutate(Wells = gsub(" ", "", Wells))

# modify the last sample of WTK1 wells to be D7-D12 instead of D7-12
df.wtk1.mapping$Wells[match("FPfSPGPG184", df.wtk1.mapping$CODissoID)] <- "D7-D12"

# add to samples table
df.samples %<>%
    bind_rows(df.wtk1.mapping)

# WTK2
df.wtk2.mapping <- read_csv(wtk2.mapping.csv)

df.wtk2.mapping %<>%
    mutate(WTK_ID = "WTK2") %>%
    select(WTK_ID,
           SampleNumber,
           CODissoID = CoDissoID,
           Wells)

# add to samples table
df.samples %<>%
    bind_rows(df.wtk2.mapping)

# WTK3
df.wtk3.mapping <- read_csv(wtk3.mapping.csv)

df.wtk3.mapping %<>%
    mutate(WTK_ID = "WTK3") %>%
    select(WTK_ID,
           SampleNumber,
           CODissoID,
           Wells)

# add to samples table
df.samples %<>%
    bind_rows(df.wtk3.mapping)

# WTK4
df.wtk4.mapping <- read_csv(wtk4.mapping.csv)

df.wtk4.mapping %<>%
    mutate(WTK_ID = "WTK4") %>%
    select(WTK_ID,
           SampleNumber,
           CODissoID,
           Wells)

# add to samples table
df.samples %<>%
    bind_rows(df.wtk4.mapping)

# WTK5
df.wtk5.mapping <- read_csv(wtk5.mapping.csv)

df.wtk5.mapping %<>%
    mutate(WTK_ID = "WTK5") %>%
    select(WTK_ID,
           SampleNumber,
           CODissoID,
           Wells)

# add to samples table
df.samples %<>%
    bind_rows(df.wtk5.mapping)

# WTK6
df.wtk6.mapping <- read_csv(wtk6.mapping.csv)

df.wtk6.mapping %<>%
    mutate(WTK_ID = "WTK6") %>%
    select(WTK_ID,
           SampleNumber,
           CODissoID,
           Wells)

# add to samples table
df.samples %<>%
    bind_rows(df.wtk6.mapping)

rm(df.wtk1.mapping, df.wtk2.mapping, df.wtk3.mapping, df.wtk4.mapping, df.wtk5.mapping, df.wtk6.mapping)

# Build Barcode Table ###############
df.barcode.one <- read_csv(barcode.one.csv)

df.barcode.one %>%
    filter(type == "T") -> df.oligoDT.barcodes

df.barcode.one %>%
    filter(type == "R") -> df.randHEX.barcodes

df.barcodes <- tibble(WTK_ID = c("WTK1", "WTK2", "WTK3", "WTK4", "WTK5", "WTK6"))

df.barcodes %<>%
    inner_join(df.barcode.one, by = character())

# loop over samples and pull barcodes from wells
df.samples$oligoDT_barcodes <- NA
df.samples$randHEX_barcodes <- NA

df.barcode.builder <- tibble()

for (i in 1:nrow(df.samples)) {

    oligoDT.barcodes <- NA
    randHEX.barcodes <- NA

    df.tmp <- NA
    # check for more than one well, get start and end well
    if (grepl("-", df.samples$Wells[i])) {
        start.well <- strsplit(df.samples$Wells[i], "-")[[1]][1]
        end.well <- strsplit(df.samples$Wells[i], "-")[[1]][2]
    } else {
        start.well <- df.samples$Wells[i]
        end.well <- df.samples$Wells[i]
    }

    # list barcodes for each sample
    oligoDT.barcodes <- df.oligoDT.barcodes$sequence[match(start.well, df.oligoDT.barcodes$well) : match(end.well, df.oligoDT.barcodes$well)]
    randHEX.barcodes <- df.randHEX.barcodes$sequence[match(start.well, df.randHEX.barcodes$well) : match(end.well, df.randHEX.barcodes$well)]

    df.samples$oligoDT_barcodes[i] <- paste(oligoDT.barcodes, collapse = ";")
    df.samples$randHEX_barcodes[i] <- paste(randHEX.barcodes, collapse = ";")

    # barcode table
    df.tmp <- tibble(sequence = c(oligoDT.barcodes, randHEX.barcodes), type = c(rep("T", length(oligoDT.barcodes)), rep("R", length(randHEX.barcodes))))
    df.tmp %<>%
        mutate(WTK_ID = df.samples$WTK_ID[i],
               SampleNumber = df.samples$SampleNumber[i],
               CODissoID = df.samples$CODissoID[i])

    df.barcode.builder %<>%
        bind_rows(df.tmp)

}

df.barcodes %<>%
    left_join(df.barcode.builder, by = c("WTK_ID", "sequence", "type"))

# Compile Sublibraries ##########
sublib.dirs <- list.dirs(sublibrary.dir, full.names = FALSE, recursive = FALSE)

df.sublibraries <- tibble(Sublibrary_ID = sublib.dirs)

df.sublibraries %<>%
    mutate(WTK_ID = sapply(strsplit(Sublibrary_ID, "_"), `[`, 1),
           Sublibrary = sapply(strsplit(Sublibrary_ID, "_"), `[`, 2))

# Save Tables ##########

write_csv(df.samples, df.samples.csv)
write_csv(df.barcodes, df.barcodes.csv)
write_csv(df.sublibraries, df.sublibraries.csv)

