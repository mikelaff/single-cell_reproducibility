# Get average expression for each cluster from Kreigsten

library(Seurat)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(here)
library(ComplexHeatmap)

library(mikelaffr)

# OUTPUT FILES #########################################################################################################
# Kreigstein average expression: 42 clusters, 5000 genes
kreig.avg.expression.output.rds <- here("results/seurat/20230306_kreigstein_avgExpression_SCT_scaleData.rds")

# PGP1 resolution 0.6: 8 clusters, 2000 genes
pgp1.avg.expression.output.rds <- here("results/seurat/20230306_pgp1_avgExpression_SCT_scaleData_res0.6.rds")
# INPUT FILES ##########################################################################################################
# integrated dataset
seurat.integrated.rds <- here("results/seurat/20230220_PGP1_ALL_cells_UNCR3_reference_integrated_seurat_object.rds")

# kreigstein seurat object
seur.kreig.rds <- "/proj/steinlab/projects/IVIV_scRNA/youngsook_pine/extData_annotation/kreigstein_primary_seuratObject_SCTv2/Kreigstein_primary_integrated.rds"

# Kreigstein average expression: 42 clusters, 5000 genes
kreig.avg.expression.output.rds <- here("results/seurat/20230306_kreigstein_avgExpression_SCT_scaleData.rds")
# Kreigstein (Bhaduri) cell types
kreig.cell.type.txt <- "/proj/steinlab/projects/IVIV_scRNA/youngsook_pine/extData_annotation/markerList/markerList/bhaduri_primary.txt"

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

# Read in Gencode #######
# load gene names
df.gencode <- read_csv(df.gencode.csv)

# # Import Kreigstein Seurat Object ################
#
# seur.kreig <- readRDS(seur.kreig.rds)
#
# seur.kreig$clustInfo
#
# Idents(seur.kreig) <- seur.kreig$clustInfo
#
# average.expression <- AverageExpression(seur.kreig, assays = "SCT", slot = "scale.data", group.by = "ident")
# str(average.expression)
#
# # Save Avg Expression #########
# saveRDS(average.expression, kreig.avg.expression.output.rds)

# Load Bhaduri Data ############
#gene.list.bhaduri<-readRDS('Seurat/bhaduri/5000.variable.genes.rds')
df.avgExp.bhaduri <- readRDS(kreig.avg.expression.output.rds)[[1]] %>% data.frame()
colnames(df.avgExp.bhaduri) <- sub('X', 'bhaduri_C', colnames(df.avgExp.bhaduri))

df.avgExp.bhaduri$gene_name <- rownames(df.avgExp.bhaduri)

cell.info <- data.table::fread(kreig.cell.type.txt, data.table=F, header=T)
cell.info$cluster <- paste0('bhaduri_C',cell.info$V1)

# Load PGP1 data ############
seur.pgp1 <- readRDS(seurat.integrated.rds)


# Cluster Integrated Data #########
DefaultAssay(seur.pgp1) <- "integrated"
seur.pgp1 <- RunPCA(seur.pgp1,
                    npcs = 20,
                    verbose = FALSE)


min.pc <- get_min_pc(seur.pgp1[["pca"]]@stdev)

seur.pgp1 <- RunUMAP(seur.pgp1,
                     dims = 1:min.pc,
                     verbose = TRUE)

seur.pgp1 <- FindNeighbors(seur.pgp1,
                           dims = 1:min.pc,
                           verbose = TRUE)

seur.pgp1 <- FindClusters(seur.pgp1,
                          verbose = TRUE,
                          resolution = 0.6)



# Get Avg Expression #############

average.expression.pgp1 <- AverageExpression(seur.pgp1, assays = "SCT", slot = "scale.data", group.by = "ident")
saveRDS(average.expression.pgp1, pgp1.avg.expression.output.rds)

df.avgExp.pgp1 <- average.expression.pgp1[[1]] %>% data.frame()

df.avgExp.pgp1 <- readRDS(pgp1.avg.expression.output.rds)[[1]] %>% data.frame()
colnames(df.avgExp.pgp1) <- sub('X', 'pgp1_C', colnames(df.avgExp.pgp1))

df.avgExp.pgp1$gene_id <- rownames(df.avgExp.pgp1)
df.avgExp.pgp1$gene_name <- df.gencode$gene_name[match(df.avgExp.pgp1$gene_id, df.gencode$gene_id)]



gene.list.IBIS <- readRDS(paste0('Seurat/post-processed/counts/20230109.keep.snRNA.ppMT.0.75.nC.1500.nF.1000/',DFiltering,'.merged.5000.variable.genes.rds'))


genelist <- intersect(df.avgExp.pgp1$gene_name, df.avgExp.bhaduri$gene_name)


avgExp.bhaduri.tmp <- df.avgExp.bhaduri[match(genelist, df.avgExp.bhaduri$gene_name),]
avgExp.pgp1.tmp <- df.avgExp.pgp1[match(genelist, df.avgExp.pgp1$gene_name),]
cor.res.tmp <- cor(avgExp.bhaduri.tmp[,1:42], avgExp.pgp1.tmp[,1:9]) %>% as.matrix()
print(max(cor.res.tmp))

cor.res.tmp <- cor.res.tmp[match(cell.info$cluster, rownames(cor.res.tmp)),]
cell.class <- cell.info$Class
cell.state <- cell.info$State
cell.type <- cell.info$Type
cell.subtype<-cell.info$Subtype

cell.class.color<-setNames(c('#5cc4c3','#2a4f72','#d1d1d1'),unique(cell.class))
cell.state.color<-setNames(c('#dd7373','#3b3561','#ead94c','#d1d1d1'),unique(cell.state))
cell.type.color<-setNames(c('#f75449','#2509b5','#6e88f5','#8dfab3','#23aa4e','#3fe72d','#fdfe61','#b29415','#f75449','#d1d1d1'),unique(cell.type))

#update Bhaduri's data
rownames(cor.res.tmp) <- paste(cell.class,cell.state,cell.type,cell.subtype,sep=' ')
rownames(cor.res.tmp) <- sub(' Outlier Outlier Outlier','',rownames(cor.res.tmp))

cell.pc.perDay<-readRDS(paste0('Seurat/post-processed/UMAP/20230109.keep/RefN_2_0fv8qgdnjs/10_10/',DFiltering,'.UMAP.cell.embeddings.clusterinfo_res',res,'.rds'))
cell.N <- seur.pgp1@meta.data %>% group_by(seurat_clusters) %>% summarize(n=n()) %>%pull(n)
donor.N <- seur.pgp1@meta.data %>% group_by(Site, Day, Rep, seurat_clusters) %>% summarize(n=n()) %>% group_by(seurat_clusters, Day)%>%summarize(n=n())
seur.pgp1@meta.data %>% group_by(Day) %>% summarize(n=n()) %>% left_join(seur.pgp1@meta.data %>% group_by(seurat_clusters,Day)%>% summarize(nc=n())) %>% data.frame() %>% mutate(pc = 100*nc/n) %>% arrange(seurat_clusters) -> cell.pc.perDay

#seur.pgp1@meta.data$pc <- perc

# ha = HeatmapAnnotation( 'Cell proportion' = anno_barplot(matrix(nc = 2, c(seur.pgp1@meta.data %>% filter(Day=='D14') %>% pull(pc), seur.pgp1@meta.data %>% filter(Day=='D84') %>% pull(pc)),dimnames=list(unique(seur.pgp1@meta.data$seurat_clusters),c('D14','D84'))), beside = TRUE, attach = TRUE, gp=gpar(fill = c('D14'='navy','D84'='orange'))),
#                         'Donor #' = anno_barplot(matrix(nc = 2, c(donor.N %>% filter(Day=='D14') %>% pull(n), donor.N %>% filter(Day=='D84') %>% pull(n)),dimnames=list(unique(donor.N$seurat_clusters),c('D14','D84'))), beside = TRUE, attach = TRUE, gp=gpar(fill = c('D14'='navy','D84'='orange'))),
#                         'Cell counts' = anno_barplot(cell.N,gp=gpar(fill='darkgreen')))
#
# ha = HeatmapAnnotation( 'Donor #' = anno_barplot(matrix(nc = 2, c(donor.N %>% filter(Day=='D14') %>% pull(n), donor.N %>% filter(Day=='D84') %>% pull(n)),dimnames=list(unique(donor.N$seurat_clusters),c('D14','D84'))), beside = TRUE, attach = TRUE, gp=gpar(fill = c('D14'='navy','D84'='orange'))),
#                         'Cell counts' = anno_barplot(cell.N,gp=gpar(fill='darkgreen')))

#Heatmap(cor.res.tmp)

Heatmap(cor.res.tmp,name="Pearson's correlation",
        #top_annotation=ha,
        right_annotation=rowAnnotation(Class=cell.class,state=cell.state,type=cell.type,
                                       col=list(Class=cell.class.color,state=cell.state.color,type=cell.type.color)),
        column_title = paste0("Res 0.8"),row_names_gp=gpar(fontsize=6),column_names_gp=gpar(fontsize=6),
        cell_fun = function(j, i, x, y, width, height, fill) {
            if(cor.res.tmp[i, j] > 0.3)
                grid.text(sprintf("%.2f", cor.res.tmp[i, j]), x, y, gp = gpar(fontsize = 6))
        }
)

pdf(paste0(dir.pdf, "20230306_ALL_cells_UNCR3_reference_0.4_clustered.pdf"), height = 7, width = 10)

ElbowPlot(seur.pgp1)

printMessage("Plotting PCA")
DimPlot(seur.pgp1, reduction = "pca", group.by = "Day")
DimPlot(seur.pgp1, reduction = "pca", group.by = "Site")
DimPlot(seur.pgp1, reduction = "pca", group.by = "WTK_ID")
FeaturePlot(seur.pgp1, reduction = "pca", features = c("percent.mt"))


DefaultAssay(seur.pgp1) <- "SCT"
printMessage("Plotting UMAPs")
DimPlot(seur.pgp1,
        label = TRUE,
        raster = FALSE)

Heatmap(cor.res.tmp,name="Pearson's correlation",
        #top_annotation=ha,
        right_annotation=rowAnnotation(Class=cell.class,state=cell.state,type=cell.type,
                                       col=list(Class=cell.class.color,state=cell.state.color,type=cell.type.color)),
        column_title = paste0("Res 0.4"),row_names_gp=gpar(fontsize=6),column_names_gp=gpar(fontsize=6),
        cell_fun = function(j, i, x, y, width, height, fill) {
            if(cor.res.tmp[i, j] > 0.3)
                grid.text(sprintf("%.2f", cor.res.tmp[i, j]), x, y, gp = gpar(fontsize = 6))
        }
)

DimPlot(seur.pgp1,
        group.by = "Day",
        label = FALSE, raster = FALSE)
DimPlot(seur.pgp1,
        group.by = "Site",
        label = FALSE, raster = FALSE)
DimPlot(seur.pgp1,
        group.by = "WTK_ID",
        label = FALSE, raster = FALSE)
#ggsave("iddrc_umap_by_site.pdf", height = 5, width = 7)
DimPlot(seur.pgp1,
        group.by = "Sample",
        label = FALSE, raster = FALSE)

#pdf("scaled_transformed_resolution_10kCells_featurePlots.pdf", height = 7, width = 10)

ptsize <- 0.1

FeaturePlot(seur.pgp1, reduction = "umap", features = c("percent.mt"), pt.size = ptsize, raster = FALSE)

FeaturePlot(seur.pgp1, reduction = "umap", features = c("ENSG00000181449.4"), pt.size = ptsize, raster = FALSE) + ggtitle("SOX2")
FeaturePlot(seur.pgp1, reduction = "umap", features = c("ENSG00000148773.14"), pt.size = ptsize, raster = FALSE) + ggtitle("MKI67")
FeaturePlot(seur.pgp1, reduction = "umap", features = c("ENSG00000007372.25"), pt.size = ptsize, raster = FALSE) + ggtitle("PAX6")

FeaturePlot(seur.pgp1, reduction = "umap", features = c("ENSG00000136535.15"), pt.size = ptsize, raster = FALSE) + ggtitle("TBR1")
FeaturePlot(seur.pgp1, reduction = "umap", features = c("ENSG00000163508.13"), pt.size = ptsize, raster = FALSE) + ggtitle("EOMES")

FeaturePlot(seur.pgp1, reduction = "umap", features = c("ENSG00000171476.22"), pt.size = ptsize, raster = FALSE) + ggtitle("HOPX")
FeaturePlot(seur.pgp1, reduction = "umap", features = c("ENSG00000134853.12"), pt.size = ptsize, raster = FALSE) + ggtitle("PDGFRA")

FeaturePlot(seur.pgp1, reduction = "umap", features = c("ENSG00000077279.20"), pt.size = ptsize, raster = FALSE) + ggtitle("DCX")
FeaturePlot(seur.pgp1, reduction = "umap", features = c("ENSG00000258947.8"), pt.size = ptsize, raster = FALSE) + ggtitle("TUBB3")
FeaturePlot(seur.pgp1, reduction = "umap", features = c("ENSG00000164600.7"), pt.size = ptsize, raster = FALSE) + ggtitle("NEUROD6")

FeaturePlot(seur.pgp1, reduction = "umap", features = c("ENSG00000127152.18"), pt.size = ptsize, raster = FALSE) + ggtitle("BCL11B")
FeaturePlot(seur.pgp1, reduction = "umap", features = c("ENSG00000119042.17"), pt.size = ptsize, raster = FALSE) + ggtitle("SATB2")

FeaturePlot(seur.pgp1, reduction = "umap", features = c("ENSG00000204531.21"), pt.size = ptsize, raster = FALSE) + ggtitle("OCT4")
FeaturePlot(seur.pgp1, reduction = "umap", features = c("ENSG00000144355.15"), pt.size = ptsize, raster = FALSE) + ggtitle("DLX1")
FeaturePlot(seur.pgp1, reduction = "umap", features = c("ENSG00000115844.11"), pt.size = ptsize, raster = FALSE) + ggtitle("DLX2")
FeaturePlot(seur.pgp1, reduction = "umap", features = c("ENSG00000231764.11"), pt.size = ptsize, raster = FALSE) + ggtitle("DLX6-AS1")

FeaturePlot(seur.pgp1, reduction = "umap", features = c("ENSG00000118271.12"), pt.size = ptsize, raster = FALSE) + ggtitle("TTR")
FeaturePlot(seur.pgp1, reduction = "umap", features = c("ENSG00000176165.12"), pt.size = ptsize, raster = FALSE) + ggtitle("FOXG1")
FeaturePlot(seur.pgp1, reduction = "umap", features = c("ENSG00000180613.11"), pt.size = ptsize, raster = FALSE) + ggtitle("GSX2")
FeaturePlot(seur.pgp1, reduction = "umap", features = c("ENSG00000106852.16"), pt.size = ptsize, raster = FALSE) + ggtitle("LHX6")





dev.off()



# Scratch

df.cellData <- as_tibble(seur.pgp1@meta.data)

df.cellData %>%
    group_by(seurat_clusters) %>%
    summarise(count_cells = n(),
              percent_cells = 100*n()/111894) -> df.tmp


summary(df.tmp)


# neuron: c0, c7, c3, c1
# proj: c5, c4, c2, c8

(23071+7396+14267+21912)/111894 # 0.596
(10564+11581+14335+151)/111894 # 0.327

seur.pgp1 <- readRDS(seurat.merged.filtered.rds)
df.cellData <- as_tibble(seur.pgp1@meta.data)

df.cellData %>%
    group_by(Sample) %>%
    summarise(count_cells = n()) -> df.tmp

summary(df.tmp)

summary(df.cellData)

