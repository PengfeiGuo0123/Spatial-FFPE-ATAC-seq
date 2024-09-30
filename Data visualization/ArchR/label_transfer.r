library(zellkonverter)
library(Seurat)
library(SingleCellExperiment)
library(ArchR)

############ process ref data
seurat_obj <- readRDS('./label_transfer/Human_thymus_all.rds')
ref <- seurat_obj
ref <- NormalizeData(ref)
ref <- FindVariableFeatures(ref)
ref <- ScaleData(ref)
ref <- RunPCA(ref)
ref <- FindNeighbors(ref, dims = 1:30)
ref <- FindClusters(ref)
ref <- RunUMAP(ref, dims = 1:30)

genes.use <- VariableFeatures(ref)


######### process query data
threads = 40
addArchRThreads(threads = threads)

addArchRGenome("hg38")

inputFiles <- './SpFFPEHThymus04.fragments.sort.bed.gz'
sampleNames <- 'Human_Thymus_ATAC'

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = sampleNames,
  minTSS = 0, #Dont set this too high because you can always increase later
  filterFrags = NULL,
  minFrags = 0,
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

saveArchRProject(ArchRProj = projCUTA, outputDirectory = paste0("Save-inTissue-", sampleNames), load = FALSE)
projCUTA <- loadArchRProject(path = paste0("Save-inTissue-", sampleNames), force = FALSE, showLogo = TRUE)

## Spatial plots
library(ggplot2)
library(patchwork)
library(dplyr)

source('./scripts/getGeneScore_ArchR.R')
source('./scripts/SpatialPlot_new.R')

## Prepare meta data
meta.data <- as.data.frame(getCellColData(ArchRProj = projCUTA))
meta.data['cellID_archr'] <- row.names(meta.data)
new_row_names <- row.names(meta.data)
new_row_names <- unlist(lapply(new_row_names, function(x) gsub(".*#","", x)))
new_row_names <- unlist(lapply(new_row_names, function(x) gsub("-.*","", x)))
row.names(meta.data) <- new_row_names

projCUTA <- addImputeWeights(projCUTA)

markerGenes <- getFeatures(ArchRProj = projCUTA, useMatrix = "GeneScoreMatrix")

gene_score <- getGeneScore_ArchR(ArchRProj = projCUTA, name = markerGenes, imputeWeights = getImputeWeights(projCUTA))
saveRDS(gene_score, './lm_gene_score_ATAC_all_gene.rds')


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

p1 <- SpatialPlot(spatial.obj, label = FALSE, label.size = 3, group.by = 'Clusters', pt.size.factor = 4.5, image.alpha = 0, stroke = 0)
p1$layers[[1]]$aes_params <- c(p1$layers[[1]]$aes_params, shape=22)
p1

png(filename = paste0(sampleNames, '_test.png'), width = 2400, height = 2400, res = 300)
p1
dev.off()


p2 <- plotEmbedding(ArchRProj = projCUTA, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", size = 1.5)
p2

plotPDF(p2, name = paste0(sampleNames, '_clusters_umap.pdf'), ArchRProj = projCUTA, addDOC = FALSE, width = 5, height = 5)


getAvailableMatrices(projCUTA)
genes_CUTA <- getFeatures(
  ArchRProj = projCUTA,
  useMatrix = "GeneScoreMatrix",
  select = NULL,
  ignoreCase = TRUE
)

genes.use <- genes.use[which(genes.use %in% genes_CUTA)]

projCUTA <- addImputeWeights(projCUTA)
imputation <- getGeneScore_ArchR(ArchRProj = projCUTA, name = genes.use, imputeWeights = getImputeWeights(projCUTA)) #name = gsub("\\-", "_", genes.use)
imputation <- log(imputation + 1) #use natural log

object_CUTA <- CreateSeuratObject(counts = imputation, assay = 'GeneScore', project = 'CUTA')
DefaultAssay(object_CUTA) <- "GeneScore"
object_CUTA <- Seurat::ScaleData(object_CUTA)


query <- object_CUTA
query <-  FindVariableFeatures(query)
query <- RunPCA(query)
query <- FindNeighbors(query, dims = 1:30)
query <- FindClusters(query)
query <- RunUMAP(query, dims = 1:30)

########## start integration
anchors <- FindTransferAnchors(reference = ref, query = query,
                               reference.assay = "originalexp", query.assay = "GeneScore", reduction = "cca")


prediction_class <- TransferData(anchorset = anchors, refdata = ref$Anno_level_1, dims = 1:30, k.weight = 50, weight.reduction='cca' )
query <- AddMetaData(query, metadata = prediction_class
                     #, col.name = 'prediction_class'
)


# add cols
n_clusters <- length(unique(ref$Anno_level_1))
cols <- ArchRPalettes$stallion[as.character(seq_len(n_clusters))]
names(cols) <-unique(ref$Anno_level_1)

# Unimodal UMAP Projection
ref <- RunUMAP(ref, dims = 1:30, reduction = "pca", return.model = TRUE)
query <- MapQuery(anchorset = anchors, reference = ref, query = query,
                  refdata = list(celltype = "Anno_level_1"), reference.reduction = "pca", reduction.model = "umap")

p1 <- DimPlot(ref, reduction = "umap", group.by = "Anno_level_1", label = TRUE, label.size = 3,
              repel = TRUE, raster=FALSE, cols = cols) + NoLegend() + ggtitle("Reference annotations")
p2 <- DimPlot(query, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE, cols = cols) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2

setwd('label_transfer')
pdf(paste0("./dim_plot_projection_with_ref_region.pdf"), width = 12, height = 6)
p1 + p2
dev.off()


################ overlied umap
# Extract UMAP coordinates from reference and query
ref_umap <- Embeddings(ref, "umap")
query_umap <- Embeddings(query, "ref.umap")

# Combine the UMAP coordinates into a single data frame for plotting
ref_umap_df <- as.data.frame(ref_umap)
ref_umap_df$group <- "Reference"
ref_umap_df$Region <- ref$Anno_level_1  # Assuming 'Region' is in the metadata


query_umap_df <- as.data.frame(query_umap)
query_umap_df$group <- "Query"
query_umap_df$Region <- query$predicted.celltype  # Assuming 'prediction_region' is in the metadata


# Combine the data frames
colnames(query_umap_df) <- c( 'UMAP_1',    'UMAP_2',     'group', 'Region')
umap_combined_df <- rbind(ref_umap_df, query_umap_df)

# Plot using ggplot2 to overlay UMAPs
library(ggplot2)


p <- ggplot() +
  # Plot the reference data first (background) with color by region
  geom_point(data = subset(umap_combined_df, group == "Reference"),
             aes(x = UMAP_1, y = UMAP_2, color = Region),
             size = 0.3) +
  # Plot the query data on top with black color
  geom_point(data = subset(umap_combined_df, group == "Query"),
             aes(x = UMAP_1, y = UMAP_2),
             color = "black",
             size = 0.3) +
  theme_minimal() +
  scale_color_manual(values = cols) +
  ggtitle("Overlay of Reference and Query UMAPs") +
  theme(legend.position = "right",
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),  # Remove minor gridlines
        panel.background = element_blank())  # Remove panel background


pdf(paste0("./dim_plot_co_projection.pdf"), width = 8, height = 6)
p
dev.off()



p <- ggplot() +
  # Plot the reference data first (background) with color by region
  geom_point(data = subset(umap_combined_df, group == "Reference"),
             aes(x = UMAP_1, y = UMAP_2),
             color = "grey",
             size = 0.3) +
  # Plot the query data on top with black color
  geom_point(data = subset(umap_combined_df, group == "Query"),
             aes(x = UMAP_1, y = UMAP_2, color = Region),
             size = 0.3) +
  theme_minimal() +
  scale_color_manual(values = cols) +
  ggtitle("Overlay of Reference and Query UMAPs") +
  theme(legend.position = "right",
        panel.grid.major = element_blank(),  # Remove major gridlines
        panel.grid.minor = element_blank(),  # Remove minor gridlines
        panel.background = element_blank())  # Remove panel background


pdf(paste0("./dim_plot_co_projection_opp.pdf"), width = 8, height = 6)
p
dev.off()

save.image(file = "lm_thymus_label_transfer.RData")




# spatial label plot
data.dir <- '/mnt/HDD1/Users/liran/05_FFPE/FFPEThymus'
assay = "Spatial"
filter.matrix = TRUE
slice = "slice1"

image <- Read10X_Image(image.dir = file.path(data.dir, "spatial"), filter.matrix = filter.matrix)
image <- image[Cells(x = query)]
DefaultAssay(query = image) <- assay
query[[slice]] <- image

p <- SpatialDimPlot(query, label = FALSE, label.size = 3, group.by = 'predicted.celltype', pt.size.factor = 4, image.alpha = 0, cols = cols)
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
p

pdf('label_spatial_image_alpha0.pdf',  width = 6, height = 6)
p
dev.off()

p <- SpatialDimPlot(query, label = FALSE, label.size = 3, group.by = 'predicted.celltype', pt.size.factor = 4, image.alpha = 0.6, cols = cols)
p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
p

pdf('label_spatial_image_alpha0.6.pdf',  width = 6, height = 6)
p
dev.off()

source('/mnt/HDD1/Users/liran/nano_review/72nanobody/SpnbRNA72/SpnbRNA72/output/align_ori/Solo.out/Gene/GeneFull/scripts/SpatialDimPlot_new.R')
Idents(query) <- 'predicted.id'

table(query$predicted.id)
ids.highlight <- names(table(query$predicted.id))
ids.highlight


# plot list
features_spatial <- names(table(query$predicted.id))


plot_features <- function(feature){
  p <- SpatialDimPlot_new(query, cells.highlight = CellsByIdentities(object = query, idents = feature),
                          facet.highlight = TRUE, pt.size.factor = 3, alpha = c(1,0.05), stroke = 0, image.alpha = 0)
  p$layers[[1]]$aes_params <- c(p$layers[[1]]$aes_params, shape=22) # set spots to square shape
  p
}

ggList <- lapply(features_spatial, plot_features)
ggList[[1]]

pdf('atac_MOCA_label_highlight.pdf')
ggList
dev.off()





