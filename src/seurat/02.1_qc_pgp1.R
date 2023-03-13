# QC pgp1 samples


#library(fishpond)
#library(data.table)
library(Seurat)
#library(DoubletFinder)
library(miQC)
#library(SeuratWrappers)
#library(flexmix)
#library(SingleCellExperiment)
#library(Matrix)
#library(stringr)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
#library(tidytext)
library(here)

#library(plyranges)

library(mikelaffr)

# OUTPUT FILES #########################################################################################################
# merged seurat filtered file, filtered for mitochondrial read percentage, gene count, and read count
seurat.merged.filtered.rds <- here("results/seurat/20230227_PGP1_QC_filtered_seurat_object.rds")
# # filtered files by sample dir, same as above but split into in
# by.sample.dir <- here("results/seurat/PGP1_QC_split_by_sample_seurat_objects/")
# dir.create(by.sample.dir, showWarnings = FALSE, recursive = TRUE)


# INPUT FILES ##########################################################################################################
# merged seurat file, unfiltered
seurat.merged.rds <- here("results/seurat/20230201_PGP1_seurat_object.rds")

# gencode gtf file
#gencode.gtf <- here("data/refgenome/gencode/gencode.v40.annotation.gtf.gz")
# gencode gene names, types, and ensgid
df.gencode.csv <- here("data/refgenome/gencode/gencode.v40.genes.csv.gz")

# GLOBALS ##############################################################################################################
# get min PCs
get_min_pc <- function(stdev) {
    sum.stdev <- sum(stdev)
    percent.stdev <- (stdev / sum.stdev) * 100
    cumulative <- cumsum(percent.stdev)

    co1 <- which(cumulative > 90 & percent.stdev < 5)[1]
    co2 <- sort(which((percent.stdev[1:length(percent.stdev) - 1] - percent.stdev[2:length(percent.stdev)]) > 0.1), decreasing = TRUE)[1] + 1

    min.pc <- min(co1, co2)

    return(min.pc)
}



# Load GTF File #######
# gencode <- read_gff(gencode.gtf)
#
# df.gencode <- as_tibble(gencode)
#
# rm(gencode)
#
# df.gencode %<>%
#     filter(type == "gene") %>%
#     select(gene_id, gene_type, gene_name, chrom = seqnames)
df.gencode <- read_csv(df.gencode.csv)

# Load Seurat Object ########
seurat.merged <- readRDS(seurat.merged.rds)
# label samples
seurat.merged@meta.data$Sample <- paste(seurat.merged@meta.data$Site, seurat.merged@meta.data$Day, seurat.merged@meta.data$Rep, sep = "_")


# percent mt reads
df.gencode %>%
    filter(chrom == "chrM") %>%
    pull(gene_id) -> ensg.mt

seurat.merged <- PercentageFeatureSet(object = seurat.merged, features = ensg.mt, col.name = "percent.mt")

seurat.merged <- RunMiQC(seurat.merged, percent.mt = "percent.mt", nFeature_RNA = "nFeature_RNA",
                         posterior.cutoff = 0.75, model.slot = "flexmix_model")


df.cellData <- as_tibble(seurat.merged@meta.data)

pdf("20230203_alevin-fry_count_filtering_pgp1.pdf", width = 14, height = 8)
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
    group_by(Sample, WTK_ID) %>%
    dplyr::summarise(cell_count = dplyr::n(),
                     num_wells = n_distinct(well)) %>%
    ggplot(aes(x = Sample, y = cell_count, fill = factor(num_wells))) +
    geom_col(size = 1, width = 0.9) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 14)) +
    scale_y_continuous(labels = scales::label_comma()) +
    scale_fill_viridis_d() +
    labs(title = "Alevin-fry Unfiltered Cells")

df.cellData %>%
    group_by(Sample, WTK_ID) %>%
    dplyr::summarise(cell_count = dplyr::n(),
                     num_wells = n_distinct(well)) %>%
    ggplot(aes(x = reorder(Sample, cell_count), y = cell_count, fill = factor(num_wells))) +
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
    group_by(Sample, WTK_ID) %>%
    dplyr::summarise(cell_count = dplyr::n(),
                     num_wells = n_distinct(well)) %>%
    ggplot(aes(x = reorder(Sample, cell_count), y = cell_count, fill = factor(num_wells))) +
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
    group_by(Sample, WTK_ID) %>%
    ggplot(aes(x = Sample, y = nCount_RNA, color = WTK_ID)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 14)) +
    scale_y_log10(labels = scales::label_comma()) +
    scale_color_manual(values = cbPalette[c(3,1,2)]) +
    labs(title = "Alevin-fry miQC Filtered Cells",
         caption = "miQC post. prob < 0.75")

df.cellData %>%
    filter(miQC.keep == "keep") %>%
    group_by(Sample, WTK_ID) %>%
    ggplot(aes(x = Sample, y = nFeature_RNA, color = WTK_ID)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 14)) +
    scale_y_log10(labels = scales::label_comma()) +
    scale_color_manual(values = cbPalette[c(3,1,2)]) +
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
# ggsave("20230203_filter_thresholds_pgp1.pdf", height = 20, width = 25)








# Filter for low reads, low features, miQC prob
#seurat.merged.filtered <- subset(seurat.merged, subset = nCount_RNA > 1500 & nFeature_RNA > 1000 & miQC.keep == "keep")

# df.cellData %>%
#     filter(nCount_RNA > 1500 & nFeature_RNA > 1000 & miQC.keep == "keep") -> df.cellData.filtered

# Filter for low reads, low features, low percent.mt
df.cellData %>%
    filter(nCount_RNA >= 1500 & nFeature_RNA >= 1000 & percent.mt < 10) -> df.cellData.filtered

# reload seurat object before filtering
seurat.merged <- readRDS(seurat.merged.rds)
# replace metadata
df.cellData <- as.data.frame(df.cellData)
rownames(df.cellData) <- df.cellData$cell_name

all(colnames(seurat.merged) == rownames(df.cellData))
seurat.merged@meta.data <- df.cellData

# filter seurat for cells to keep
seurat.merged.filtered <- seurat.merged[,df.cellData.filtered$cell_name]

dim(seurat.merged.filtered)

rm(seurat.merged, df.cellData)


df.cellData.filtered %>%
    group_by(Sample, WTK_ID) %>%
    dplyr::summarise(cell_count = dplyr::n(),
                     num_wells = n_distinct(well)) %>%
    ggplot(aes(x = reorder(Sample, cell_count), y = cell_count, fill = factor(num_wells))) +
    geom_col(size = 1, width = 0.9) +
    geom_text(data = function(x) subset(x, cell_count < max(cell_count)), mapping = aes(label = formatC(cell_count, format = "d", big.mark = ",")), vjust = 0.5, angle = 90, hjust = 0) +
    geom_text(data = function(x) subset(x, cell_count == max(cell_count)), mapping = aes(label = formatC(cell_count, format = "d", big.mark = ",")), vjust = 0.5, angle = 90, hjust = 1, color = "white") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 14)) +
    scale_y_continuous(labels = scales::label_comma()) +
    scale_fill_viridis_d() +
    labs(title = "Alevin-fry miQC/Count Filtered Cells",
         caption = "percent.mt < 10 & nFeature_RNA >= 1000 & nCount_RNA >= 1500")


df.cellData.filtered %>%
    group_by(Sample, WTK_ID) %>%
    ggplot(aes(x = Sample, y = nCount_RNA, color = WTK_ID)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 14)) +
    scale_y_log10(labels = scales::label_comma()) +
    scale_color_manual(values = cbPalette[c(3,1,2)]) +
    labs(title = "Alevin-fry miQC/Count Filtered Cells",
         caption = "percent.mt < 10 & nFeature_RNA >= 1000 & nCount_RNA >= 1500")

df.cellData.filtered %>%
    group_by(Sample, WTK_ID) %>%
    ggplot(aes(x = Sample, y = nFeature_RNA, color = WTK_ID)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.y = element_text(size = 14)) +
    scale_y_log10(labels = scales::label_comma()) +
    scale_color_manual(values = cbPalette[c(3,1,2)]) +
    labs(title = "Alevin-fry miQC/Count Filtered Cells",
         caption = "percent.mt < 10 & nFeature_RNA >= 1000 & nCount_RNA >= 1500")


dev.off()

# Sample Filter ############
# remove samples with less than 100 cells
df.cellData.filtered %>%
    group_by(Sample) %>%
    dplyr::summarise(cell_count = dplyr::n()) %>%
    filter(cell_count > 100) %>%
    pull(Sample) -> samples.to.keep

seurat.merged.filtered <- subset(seurat.merged.filtered, subset = Sample %in% samples.to.keep)
dim(seurat.merged.filtered)


# save filtered data
saveRDS(seurat.merged.filtered, seurat.merged.filtered.rds)










# Scratch ################

# seur.iviv <- readRDS("/proj/steinlab/projects/IVIV_scRNA/IBIS/data/20230109.keep.snRNA.ppMT.0.75.nC.1500.nF.1000.QCed.rds")
#
# df.cellData.ivi <- as_tibble(seur.iviv@meta.data)
#
# df.cellData.ivi %>%
#     group_by(sampleNames, WTK) %>%
#     dplyr::summarise(cell_count = dplyr::n()) %>%
#     ggplot(aes(x = reorder(sampleNames, cell_count), y = cell_count, fill = factor(WTK))) +
#     geom_col(size = 1, width = 0.9) +
#     geom_text(data = function(x) subset(x, cell_count < max(cell_count)), mapping = aes(label = formatC(cell_count, format = "d", big.mark = ",")), vjust = 0.5, angle = 90, hjust = 0) +
#     geom_text(data = function(x) subset(x, cell_count == max(cell_count)), mapping = aes(label = formatC(cell_count, format = "d", big.mark = ",")), vjust = 0.5, angle = 90, hjust = 1, color = "white") +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 14),
#           axis.text.y = element_text(size = 14),
#           axis.title.y = element_text(size = 14)) +
#     scale_y_continuous(labels = scales::label_comma()) +
#     scale_fill_viridis_d() +
#     labs(title = "Alevin-fry miQC/Count Filtered Cells",
#          caption = "miQC post. prob < 0.75 & nFeature_RNA > 1000 & nCount_RNA > 1500")
# #
# # col.names.pgp1 <- colnames(seur.iviv)[grepl("pgp1", colnames(seur.iviv))]
#
#
# df.cellData.ivi %>%
#     ggplot(aes(x = percent.mt, fill = WTK)) +
#     geom_histogram(position = "dodge")

# for (i in 1:length(seur.list)) {
#     printMessage(paste(i, "of", length(seur.list), ":", names(seur.list)[i]))
#
#     wtk.id <- seur.list[[i]]@meta.data$WTK_ID[1]
#
#     if (wtk.id %in% c("WTK1", "WTK2", "WTK3")) {
#         exp.D.ratio <- 0.034
#     } else {
#         exp.D.ratio <- 0.0425
#     }
#
#     seur.list[[i]]@meta.data$DF.name <- "Singlet"
#
#     df.cells <- as_tibble(seur.list[[i]]@meta.data)
#     df.cells %<>%
#         arrange(nFeature_RNA) %>%
#         slice_head(prop = 0.8)
#
#     df.cells %>%
#         slice_sample(prop = exp.D.ratio) %>%
#         pull(cell_name) -> doublet.cells
#
#     seur.list[[i]]@meta.data[doublet.cells,]$DF.name <- "Doublet"
#
#
#     #VlnPlot(seur.list[[i]], features = "nFeature_RNA", group.by = "DF.name", pt.size = 0.4) + ggtitle(names(seur.list)[i])
#
#
# }
#
# pdf("doublets2.pdf")
# for (i in 1:length(seur.list)) {
#     print(VlnPlot(seur.list[[i]], features = "nFeature_RNA", group.by = "DF.name", pt.size = 0.4) + ggtitle(names(seur.list)[i]))
# }
# dev.off()
