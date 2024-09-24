m(list=ls())
library(ArchR)
library(Seurat)
library(grid)
library(ggplot2)

threads = 40
addArchRThreads(threads = threads)

addArchRGenome("mm10")

inputFiles <- './FFPE_S1.fragments.sort.bed.gz'
sampleNames <- 'S1_ATAC'

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = sampleNames,
  minTSS = 1, #Dont set this too high because you can always increase later
  filterFrags = NULL,
  minFrags = 10,
  maxFrags = 1e+07,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  offsetPlus = 0,
  offsetMinus = 0,
  TileMatParams = list(tileSize = 5000), force =  TRUE
)
ArrowFiles

projHeme1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = sampleNames,
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)
projHeme1

## Prepare meta data
meta.data <- as.data.frame(getCellColData(ArchRProj = projHeme1))
meta.data['cellID_archr'] <- row.names(meta.data)
new_row_names <- row.names(meta.data)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub(".*#","", x)))
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(meta.data) <- new_row_names

# qc_10x <- read.csv('./singlecell.csv', header = TRUE)
# qc_10x['Proportion_of_TSS_fragments'] <- qc_10x$TSS_fragments / qc_10x$passed_filters
# qc_10x['Proportion_of_mito_fragments'] <- qc_10x$mitochondrial / qc_10x$total
# qc_10x['FRiP'] <- qc_10x$peak_region_fragments / qc_10x$passed_filters
# 
# summary(qc_10x$Proportion_of_TSS_fragments)

# qc_10x <- qc_10x[-1,]
# row.names(qc_10x) <- qc_10x$barcode
# new_row_names <- row.names(qc_10x)
# new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
# row.names(qc_10x) <- new_row_names
# qc_10x <- qc_10x[, c('Proportion_of_TSS_fragments', 'Proportion_of_mito_fragments'), drop =FALSE]
# qc_10x <- qc_10x[, c('Proportion_of_TSS_fragments', 'Proportion_of_mito_fragments', 'FRiP'), drop =FALSE]
# 
# meta.data <- merge(meta.data, qc_10x, by=0, all=TRUE)
# row.names(meta.data) <- meta.data$Row.names
# meta.data <- meta.data[-1]


data.dir <- getwd()
assay = "Spatial"
filter.matrix = TRUE
slice = "slice1"

image <- Read10X_Image(image.dir = file.path(data.dir, "spatial"), 
                       filter.matrix = filter.matrix)

# export meta data
meta.data.spatial <- meta.data[row.names(image@coordinates), ]
summary(meta.data.spatial$nFrags)
summary(meta.data.spatial$Proportion_of_TSS_fragments)
summary(meta.data.spatial$Proportion_of_mito_fragments)
summary(meta.data.spatial$FRiP)

write.table(meta.data.spatial, paste0(sampleNames, '_meta_data.txt'), row.names = TRUE, col.names = TRUE)


projCUTA <- projHeme1[meta.data.spatial$cellID_archr, ]
projCUTA

p <- plotFragmentSizes(ArchRProj = projCUTA)
png(filename = 'FragmentSizes.png', width = 1200, height = 1200, res = 300)
p
dev.off()

projCUTA <- addIterativeLSI(
  ArchRProj = projCUTA,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,
  force = TRUE
)

projCUTA <- addClusters(
  input = projCUTA,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 1, # default:0.5
  force = TRUE
)

projCUTA <- addUMAP(
  ArchRProj = projCUTA, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)

projCUTA <- addImputeWeights(projCUTA)

#start from here
#saveArchRProject(ArchRProj = projCUTA, outputDirectory = paste0("Save-inTissue-", sampleNames), load = FALSE)
projCUTA <- loadArchRProject(path = paste0("Save-inTissue-", sampleNames), force = FALSE, showLogo = TRUE)
projCUTA

# QC plot
df <- getCellColData(projCUTA, select = c("log10(nFrags)", "TSSEnrichment"))
df

p <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), 5.5),
  ylim = c(1, 8),
  baseSize = 12
)
p
#dev.off()

plotPDF(p, name = "TSS-vs-Frags.pdf", ArchRProj = projCUTA, addDOC = FALSE)

#identifying marker genes
markersGS <- getMarkerFeatures(
  ArchRProj = projCUTA, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  testMethod = "wilcoxon"
)

markerList_pos <- getMarkers(markersGS, cutOff = "FDR < 0.05 & Log2FC >= 0.25")

markerGenes <- list()
for (i in seq_len(length(markerList_pos))) {
  markerGenes <- c(markerGenes, markerList_pos[[i]]$name)
}
markerGenes <- unlist(markerGenes)
markerGenes <- getFeatures(ArchRProj = projCUTA, useMatrix = "GeneScoreMatrix")

## Spatial plots
library(ggplot2)
library(patchwork)
library(dplyr)

source('./scripts/getGeneScore_ArchR.R')
source('./scripts/SpatialPlot_new.R')
#source('SpatialPlot_new_K27me3.R')

## Prepare meta data
meta.data <- as.data.frame(getCellColData(ArchRProj = projCUTA))
meta.data['cellID_archr'] <- row.names(meta.data)
new_row_names <- row.names(meta.data)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub(".*#","", x)))
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(meta.data) <- new_row_names

projCUTA <- addImputeWeights(projCUTA)
gene_score <- getGeneScore_ArchR(ArchRProj = projCUTA, name = markerGenes, imputeWeights = getImputeWeights(projCUTA))


data.dir <- getwd()
assay = "Spatial"
filter.matrix = TRUE
slice = "slice1"

object <- CreateSeuratObject(counts = gene_score, assay = assay, meta.data = meta.data)

image <- Read10X_Image(image.dir = file.path(data.dir, "spatial"), filter.matrix = filter.matrix)
image <- image[Cells(x = object)]
DefaultAssay(object = image) <- assay
object[[slice]] <- image

spatial.obj <- object

#p <- VlnPlot(spatial.obj, features = "nFrags", pt.size = 0.1, log = TRUE) + NoLegend()
median(meta.data$nFrags)
# png(filename = 'nFrags.png', width = 1200, height = 1200, res = 300)
# p
# dev.off()

# SpatialPlot(spatial.obj, features = "nFrags",  pt.size.factor = 4, min.cutoff = "q10", max.cutoff = "q90", image.alpha = 0) + 
#      theme(legend.position = "right")

p <- SpatialPlot_new(spatial.obj, features = "nFrags",  pt.size.factor = 4.5, min.cutoff = "q10", max.cutoff = "q90", image.alpha = 0, stroke = 0) + 
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
p
png(filename = paste0(sampleNames, '_nFrags_spatial.png'), width = 2400, height = 2400, res = 300)
p
dev.off()

#SpatialDimPlot(spatial.obj, label = FALSE, label.size = 3, group.by = 'Clusters', pt.size.factor = 4)
#Idents(spatial.obj) <- 'Clusters'
#SpatialDimPlot(spatial.obj, cells.highlight = CellsByIdentities(object = spatial.obj, idents = c('C1', 'C2')), facet.highlight = TRUE, ncol = 2, pt.size.factor = 3)

n_clusters <- length(unique(projCUTA$Clusters))
cols <- ArchRPalettes$stallion[as.character(seq_len(n_clusters))]
names(cols) <- paste0('C', seq_len(n_clusters))
#names(cols) <- paste0('C', n_clusters-seq_len(n_clusters)+1)
cols
p1 <- SpatialPlot(spatial.obj, label = FALSE, label.size = 3, group.by = 'Clusters', pt.size.factor = 4.5, cols = cols, image.alpha = 0, stroke = 0)
p1$layers[[1]]$aes_params <- c(p1$layers[[1]]$aes_params, shape=22)
p1

png(filename = paste0(sampleNames, '_clusters_spatial.png'), width = 3600, height = 3000, res = 300)
plotPDF(p1, name = paste0(sampleNames, '_clusters_spatial.pdf'), addDOC = FALSE, width = 5, height = 5)
p1
dev.off()

p2 <- plotEmbedding(ArchRProj = projCUTA, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", size = 1.5, pal = cols)
p2

plotPDF(p2, name = paste0(sampleNames, '_clusters_umap.pdf'), ArchRProj = projCUTA, addDOC = FALSE, width = 5, height = 5)

png(filename = paste0(sampleNames, '_clusters_umap.png'), width = 3600, height = 3000, res = 300)
p2
dev.off()

markers_cluster <- as.data.frame(markerList_pos$C2)
markers_cluster

features_spatial <- markerGenes
feature <- features_spatial[4]
feature <- 'Bcl11b'

p <- SpatialPlot_new(spatial.obj, features = feature, slot = 'counts', pt.size.factor = 4.5, image.alpha = 0, stroke = 0, min.cutoff = "q10", max.cutoff = "q99") +
  theme(legend.position = "right", legend.text=element_text(size=15), legend.title=element_text(size=15))
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
p

png(filename = paste0('./markers_plot_liver/', feature, '_upReg_spatial.png'), width = 1200, height = 1200, res = 300)
p
dev.off()


p_track <- plotBrowserTrack(
  ArchRProj = projCUTA,
  groupBy = "Clusters",
  geneSymbol = feature,
  upstream = 50000, #50000,
  downstream = 50000, #50000,
  normMethod = 'ReadsInTSS',
  minCells = 0,
  borderWidth = 0
)

# p_track <- plotBrowserTrack(
#   ArchRProj = projCUTA,
#   plotSummary = c("bulkTrack", "scTrack", "geneTrack"),
#   sizes = c(6, 6, 2),
#   groupBy = "Clusters",
#   geneSymbol = features_spatial,
#   upstream = 50000,
#   downstream = 50000,
#   scCellsMax = 100,
#   normMethod = 'nFrags'
# )

grid::grid.newpage()
png(filename = 'Bcl11b_upReg_track.png', width = 1200, height = 1200, res = 300)
grid::grid.draw(p_track$Bcl11b)
dev.off()
#
plotPDF(plotList = p_track,
        name = paste0(sampleNames, "_Plot-Tracks-Marker-Genes-10kb.pdf"),
        ArchRProj = projCUTA,
        addDOC = FALSE, width = 5, height = 5)
