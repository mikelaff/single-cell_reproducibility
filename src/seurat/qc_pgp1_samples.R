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
library(tidytext)
library(here)

library(plyranges)

library(mikelaffr)

# OUTPUT FILES #########################################################################################################

# merged seurat filtered file
seurat.merged.filtered.rds <- here("results/seurat/20230126_PGP1_filtered_seurat_object.rds")

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

seurat.merged <- RunMiQC(seurat.merged, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA",
                         posterior.cutoff = 0.75, model.slot = "flexmix_model")


df.cellData <- as_tibble(seurat.merged@meta.data)

pdf("alevin-fry_count_filtering.pdf", width = 14, height = 8)
df.cellData %>%
    slice_sample(prop = 0.1) %>%
    arrange(percent.mt) %>%
    ggplot(aes(x = nCount_RNA, y = nFeature_RNA, color = percent.mt)) +
    geom_point(size = 1, alpha = 0.5, shape = 20) +
    scale_color_viridis_c() +
    labs(title = "Alevin-fry Unfiltered Cells")

df.cellData %>%
    slice_sample(prop = 0.1) %>%
    arrange(percent.mt) %>%
    ggplot(aes(x = percent.mt, y = miQC.probability, color = miQC.keep)) +
    geom_point(size = 1) +
    scale_color_manual(values = cbPalette[c(6,7)]) +
    labs(title = "Alevin-fry Unfiltered Cells")

df.cellData %>%
    slice_sample(prop = 0.1) %>%
    arrange(miQC.probability) %>%
    ggplot(aes(x = nFeature_RNA, y = percent.mt, color = miQC.probability)) +
    geom_point(size = 1, alpha = 0.5, shape = 20) +
    scale_color_viridis_c() +
    labs(title = "Alevin-fry Unfiltered Cells")

df.cellData %>%
    slice_sample(prop = 0.1) %>%
    arrange(miQC.probability) %>%
    ggplot(aes(x = nFeature_RNA, y = percent.mt, color = miQC.keep)) +
    geom_point(size = 1, alpha = 0.5, shape = 20) +
    scale_color_manual(values = cbPalette) +
    labs(title = "Alevin-fry Unfiltered Cells")


df.cellData %>%
    mutate(SampleName = paste(Site, Day, Rep, sep = "_")) %>%
    group_by(SampleName, WTK_ID) %>%
    dplyr::summarise(cell_count = dplyr::n(),
                     num_wells = n_distinct(well)) %>%
    ggplot(aes(x = SampleName, y = cell_count, fill = factor(num_wells))) +
    geom_col(size = 1, width = 0.9) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 14)) +
    scale_y_continuous(labels = scales::label_comma()) +
    scale_fill_viridis_d() +
    labs(title = "Alevin-fry Unfiltered Cells")

df.cellData %>%
    mutate(SampleName = paste(Site, Day, Rep, sep = "_")) %>%
    group_by(SampleName, WTK_ID) %>%
    dplyr::summarise(cell_count = dplyr::n(),
                     num_wells = n_distinct(well)) %>%
    ggplot(aes(x = reorder(SampleName, cell_count), y = cell_count, fill = factor(num_wells))) +
    geom_col(size = 1, width = 0.9) +
    geom_text(data = function(x) subset(x, cell_count == min(cell_count)), mapping = aes(label = formatC(cell_count, format = "d", big.mark = ",")), vjust = 0.5, angle = 90, hjust = 0) +
    geom_text(data = function(x) subset(x, cell_count == max(cell_count)), mapping = aes(label = formatC(cell_count, format = "d", big.mark = ",")), vjust = 0.5, angle = 90, hjust = 1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 14)) +
    scale_y_continuous(labels = scales::label_comma()) +
    scale_fill_viridis_d() +
    labs(title = "Alevin-fry Unfiltered Cells")

df.cellData %>%
    filter(miQC.keep == "keep") %>%
    mutate(SampleName = paste(Site, Day, Rep, sep = "_")) %>%
    group_by(SampleName, WTK_ID) %>%
    dplyr::summarise(cell_count = dplyr::n(),
                     num_wells = n_distinct(well)) %>%
    ggplot(aes(x = reorder(SampleName, cell_count), y = cell_count, fill = factor(num_wells))) +
    geom_col(size = 1, width = 0.9) +
    geom_text(data = function(x) subset(x, cell_count == min(cell_count)), mapping = aes(label = formatC(cell_count, format = "d", big.mark = ",")), vjust = 0.5, angle = 90, hjust = 0) +
    geom_text(data = function(x) subset(x, cell_count == max(cell_count)), mapping = aes(label = formatC(cell_count, format = "d", big.mark = ",")), vjust = 0.5, angle = 90, hjust = 1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 14)) +
    scale_y_continuous(labels = scales::label_comma()) +
    scale_fill_viridis_d() +
    labs(title = "Alevin-fry miQC Filtered Cells",
         caption = "miQC post. prob < 0.75")


df.cellData %>%
    filter(miQC.keep == "keep") %>%
    mutate(SampleName = paste(Site, Day, Rep, sep = "_")) %>%
    group_by(SampleName, WTK_ID) %>%
    ggplot(aes(x = SampleName, y = nCount_RNA, color = WTK_ID)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 14)) +
    scale_y_log10(labels = scales::label_comma()) +
    scale_color_manual(values = cbPalette) +
    labs(title = "Alevin-fry miQC Filtered Cells",
         caption = "miQC post. prob < 0.75")

df.cellData %>%
    filter(miQC.keep == "keep") %>%
    mutate(SampleName = paste(Site, Day, Rep, sep = "_")) %>%
    group_by(SampleName, WTK_ID) %>%
    ggplot(aes(x = SampleName, y = nFeature_RNA, color = WTK_ID)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 14)) +
    scale_y_log10(labels = scales::label_comma()) +
    scale_color_manual(values = cbPalette) +
    labs(title = "Alevin-fry miQC Filtered Cells",
         caption = "miQC post. prob < 0.75")


# # find thresholds
# threshold.nCount <- c(0, 100, 200, 250, 500, 750, 1000, 1500)
# threshold.nFeature <- c(0, 100, 200, 250, 500, 750, 1000, 1500)
#
# df.filtered <- tibble()
#
# for (nCount in threshold.nCount) {
#     for (nFeature in threshold.nFeature) {
#
#         df.cellData %>%
#             filter(miQC.keep == "keep" & nFeature_RNA > nFeature & nCount_RNA > nCount) %>%
#             mutate(SampleName = paste(Site, Day, Rep, sep = "_"),
#                    nFeature_thresh = nFeature,
#                    nCount_thresh = nCount) %>%
#             group_by(SampleName, Site, Day, WTK_ID, nFeature_thresh, nCount_thresh) %>%
#             dplyr::summarise(cell_count = dplyr::n(),
#                              num_wells = n_distinct(well)) %>%
#             bind_rows(df.filtered) -> df.filtered
#
#     }
# }
#
# df.filtered %>%
#     group_by(nFeature_thresh, nCount_thresh) %>%
#     ggplot(aes(x = reorder_within(SampleName, cell_count, list(nFeature_thresh, nCount_thresh)), y = cell_count, fill = Site, color = Day)) +
#     geom_col(size = 0.4, width = 0.8) +
#     facet_wrap(~nFeature_thresh + nCount_thresh, scales = "free", labeller = labeller(.rows = label_both)) +
#     #geom_text(data = function(x) subset(x, cell_count == min(cell_count)), mapping = aes(label = formatC(cell_count, format = "d", big.mark = ",")), vjust = 0.5, angle = 90, hjust = 0) +
#     #geom_text(data = function(x) subset(x, cell_count == max(cell_count)), mapping = aes(label = formatC(cell_count, format = "d", big.mark = ",")), vjust = 0.5, angle = 90, hjust = 1) +
#     theme(axis.text.x = element_blank(),
#           strip.text = element_text(size = 5, margin = margin(b = 0, t = 0))) +
#     scale_y_log10() +
#     scale_x_reordered() +
#     scale_fill_manual(values = cbPalette) +
#     scale_color_manual(values = cbPalette[c(6,4,7)]) +
#     labs(title = "Alevin-fry miQC Filtered Cells",
#          caption = "miQC post. prob < 0.75")
#
# #ggsave("filter_thresholds.pdf", height = 20, width = 25)

#hi







# Filter for low reads, low features, miQC prob
#seurat.merged.filtered <- subset(seurat.merged, subset = nCount_RNA > 1500 & nFeature_RNA > 1000 & miQC.keep == "keep")

df.cellData %>%
    filter(nCount_RNA > 1500 & nFeature_RNA > 1000 & miQC.keep == "keep") -> df.cellData.filtered

# reload seurat object for filtering
rm(seurat.merged)
seurat.merged <- readRDS(seurat.merged.rds)
# filter seurat for cells to keep
seurat.merged.filtered <- seurat.merged[,df.cellData.filtered$cell_name]

dim(seurat.merged.filtered)

rm(seurat.merged)

# save seurat filtered object
saveRDS(seurat.merged.filtered, seurat.merged.filtered.rds)



df.cellData.filtered %>%
    mutate(SampleName = paste(Site, Day, Rep, sep = "_")) %>%
    group_by(SampleName, WTK_ID) %>%
    dplyr::summarise(cell_count = dplyr::n(),
                     num_wells = n_distinct(well)) %>%
    ggplot(aes(x = reorder(SampleName, cell_count), y = cell_count, fill = factor(num_wells))) +
    geom_col(size = 1, width = 0.9) +
    geom_text(data = function(x) subset(x, cell_count < max(cell_count)), mapping = aes(label = formatC(cell_count, format = "d", big.mark = ",")), vjust = 0.5, angle = 90, hjust = 0) +
    geom_text(data = function(x) subset(x, cell_count == max(cell_count)), mapping = aes(label = formatC(cell_count, format = "d", big.mark = ",")), vjust = 0.5, angle = 90, hjust = 1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 14)) +
    scale_y_continuous(labels = scales::label_comma()) +
    scale_fill_viridis_d() +
    labs(title = "Alevin-fry miQC/Count Filtered Cells",
         caption = "miQC post. prob < 0.75 & nFeature_RNA > 1000 & nCount_RNA > 1500")


df.cellData.filtered %>%
    mutate(SampleName = paste(Site, Day, Rep, sep = "_")) %>%
    group_by(SampleName, WTK_ID) %>%
    ggplot(aes(x = SampleName, y = nCount_RNA, color = WTK_ID)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 14)) +
    scale_y_log10(labels = scales::label_comma()) +
    scale_color_manual(values = cbPalette) +
    labs(title = "Alevin-fry miQC/Count Filtered Cells",
         caption = "miQC post. prob < 0.75 & nFeature_RNA > 1000 & nCount_RNA > 1500")

df.cellData.filtered %>%
    mutate(SampleName = paste(Site, Day, Rep, sep = "_")) %>%
    group_by(SampleName, WTK_ID) %>%
    ggplot(aes(x = SampleName, y = nFeature_RNA, color = WTK_ID)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 14)) +
    scale_y_log10(labels = scales::label_comma()) +
    scale_color_manual(values = cbPalette) +
    labs(title = "Alevin-fry miQC/Count Filtered Cells",
         caption = "miQC post. prob < 0.75 & nFeature_RNA > 1000 & nCount_RNA > 1500")


dev.off()





