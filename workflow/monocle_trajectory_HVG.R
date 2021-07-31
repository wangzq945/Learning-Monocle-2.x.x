# Monocle ----

# package
library(monocle)
library(Seurat)
library(ggplot2)

# path
c("Monocle", "trajectory", "Seurat") %>% walk(function(x){if(!dir.exists(x)) {dir.create(x)}})

## Load data ----

exprs <- readRDS("data/packer_embryo/packer_embryo_expression.rds")
phenoData <- readRDS("data/packer_embryo/packer_embryo_colData.rds")
featureData <- readRDS("data/packer_embryo/packer_embryo_rowData.rds")

## Seurat workflow ----

seu <- CreateSeuratObject(counts = exprs, project = "packer_embryo", min.cells = 3, min.features = 200, meta.data = phenoData)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
ElbowPlot(seu, ndims = 50)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.5)
seu <- RunUMAP(seu, dims = 1:30)
saveRDS(seu, file = "Seurat/seu.rds")

## Generate a cell_data_set ----

if(F){
  exprs <- seu@assays$RNA@counts # sparse matrix from seurat object
  phenoData <- seu@meta.data
  featureData <- data.frame(symbol = rownames(exprs)) %>% mutate(gene_short_name = symbol) %>% column_to_rownames(var = "symbol")
}

pd <- new("AnnotatedDataFrame", data = phenoData)
fd <- new("AnnotatedDataFrame", data = featureData)
cds <- newCellDataSet(
  exprs, 
  phenoData = pd, 
  featureData = fd, 
  expressionFamily=negbinomial.size()
)

## Estimate size factors and dispersions ----

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

## Choose high variable features  ----

HVG <- seu@assays$RNA@var.features

ordering_genes <- HVG

cds <- setOrderingFilter(cds, ordering_genes)

plot_ordering_genes(cds)

## Reduce data dimensionality ----

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')

## Order cells along the trajectory ----

cds <- orderCells(cds)

plot_cell_trajectory(cds, color_by = "Pseudotime")
plot_cell_trajectory(cds, color_by = "State")
plot_cell_trajectory(cds, color_by = "cell.type")
plot_cell_trajectory(cds, color_by = "embryo.time")
plot_cell_trajectory(cds, color_by = "embryo.time.bin")

## Order cells along the trajectory with a root ----

GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$embryo.time.bin)[,"130-170"]
    return(as.numeric(names(T0_counts)[which(T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

cds <- orderCells(cds, root_state = GM_state(cds))

plot_cell_trajectory(cds, color_by = "Pseudotime")
plot_cell_trajectory(cds, color_by = "State")
plot_cell_trajectory(cds, color_by = "cell.type")
plot_cell_trajectory(cds, color_by = "embryo.time")
plot_cell_trajectory(cds, color_by = "embryo.time.bin")

## Save ----

saveRDS(cds, file = str_c("Monocle/cds.rds"))
