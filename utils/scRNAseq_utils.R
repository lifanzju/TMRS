library(Seurat)
library(DoubletFinder)
library(tidyverse)
library(ggplotify)

QC <- function(data, min_genes, max_fold, max_percent_mt, npcs, res) {
  data$percent.mt <- PercentageFeatureSet(data, pattern = "^MT-") ## Compute mitochondrial gene percentage
  # Filter cells based on QC
  max_genes <- max_fold * median(data$nFeature_RNA)
  data <- subset(data, subset = nFeature_RNA >= min_genes & nFeature_RNA <= max_genes & percent.mt < max_percent_mt)
  #### Pre-process Seurat object
  data <- NormalizeData(data)
  data <- ScaleData(data, verbose = T)
  data <- FindVariableFeatures(data)
  data <- RunPCA(data, verbose = T)
  data <- RunUMAP(data, reduction = "pca", dims = 1:npcs)
  data <- FindNeighbors(data, reduction = "pca", dims = 1:npcs)
  data <- FindClusters(data, resolution = res)
  ### pK Identification
  sweep.res.list <- paramSweep_v3(data, PCs = 1:npcs, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  ### Homotypic Doublet Proportion Estimate
  annotations <- data@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  DoubletRate <- ncol(data) * 8 * 1e-6 
  nExp_poi <- round(DoubletRate * length(data$seurat_clusters))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop)) 
  #### Run DoubletFinder
  data <- doubletFinder_v3(data, PCs = 1:npcs, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = F)
  #### Filter out predicted doublets
  Idents(data) <- select(data@meta.data, contains("DF.classifications"))
  newdata <- subset(data, idents = "Singlet")
  cat("Task", data@project.name, "completed.\n")
  return(newdata)
}

cal_pathway_auc <- function(data, genesets) {
  exp_matrix <- data@assays$RNA@counts
  cells_rankings <- AUCell_buildRankings(exp_matrix)
  cells_AUC <- AUCell_calcAUC(geneSets = genesets, rankings = cells_rankings, aucMaxRank = nrow(cells_rankings) * 0.1)
  auc <- getAUC(cells_AUC)
  return(auc)
}
