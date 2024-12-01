library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(harmony)
library(DoubletFinder)
library(stringr)
library(cowplot)
library(ggsci)
library(ggpubr)
library(data.table)
library(qs)

rm(list = ls())
Sample.dir <- "data/GSE156632"
intermediate_data.dir <- "results/step_1_hdwgcna_wgcna/scRNAseq_preprocessing/intermediate_data"
results.dir <- "results/step_1_hdwgcna_wgcna/scRNAseq_preprocessing/results/qc_results"
if (!dir.exists(results.dir)) {
  dir.create(results.dir, recursive = T)
}
source("scripts/utils/scRNAseq_utils.R")
source("scripts/utils/plot_settings.R")
Samples <- list.files(Sample.dir)
sceList <- lapply(Samples, function(x) {
  path <- file.path(Sample.dir, x)
  project_name <- substr(x, 12, 16)
  data <- fread(path, header = T) %>% data.frame(row.names = 1)
  sce <- CreateSeuratObject(
    counts = data,
    project = project_name
  )
  return(sce)
})
n <- length(sceList)

# quality control sample by sample
min_genes <- 500
max_fold <- 2
max_percent_mt <- 10
npcs <- 30
res <- 0.5
sceList_new <- lapply(sceList, QC,
  min_genes = min_genes,
  max_fold = max_fold,
  max_percent_mt = max_percent_mt,
  npcs = npcs,
  res = res
)

# merge samples
sce_merged <- merge(sceList_new[[1]], y = sceList_new[c(2:length(sceList))], add.cell.ids = substr(Samples, 12, 16))
sce_merged$group <- if_else(substr(sce_merged$orig.ident, 5, 5) == "t", "Tumor", "Normal")
qsave(sce_merged, paste0(intermediate_data.dir, "/merged_data.qs"))

# check qc results
Idents(sce_merged) <- "orig.ident"

plot1 <- FeatureScatter(sce_merged, feature1 = "nCount_RNA", feature2 = "percent.mt", plot.cor = F, group.by = "orig.ident", shuffle = T)
plot2 <- FeatureScatter(sce_merged, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", plot.cor = F, group.by = "orig.ident", shuffle = T)
plot3 <- FeatureScatter(sce_merged, feature1 = "nFeature_RNA", feature2 = "percent.mt", plot.cor = F, group.by = "orig.ident", shuffle = T)
qc_plot <- (plot1 | plot2 | plot3) + plot_layout(guides = "collect") & my_plot_theme(x.text.angle = 45) & scale_color_manual(values = adaptive_palette(12))
qsave(qc_plot, paste0(results.dir, "/qc_plot.qs"))
ggsave(paste0(results.dir, "/qc_plot.tiff"), qc_plot, width = 8, height = 3.5)

# preprocessing merged data
sce_merged <- NormalizeData(sce_merged, normalization.method = "LogNormalize", scale.factor = 1e4)
sce_merged <- FindVariableFeatures(sce_merged)
sce_merged <- ScaleData(sce_merged)
sce_merged <- RunPCA(sce_merged, features = VariableFeatures(object = sce_merged))
DimPlot(sce_merged, reduction = "pca")
# remove batch effect by harmony
sce_merged <- RunHarmony(sce_merged, group.by.vars = "orig.ident", plot_convergence = T)
# determine the pcs
ElbowPlot(sce_merged, ndims = 50)
# clustering
sce_merged <- FindNeighbors(sce_merged, dims = 1:30, reduction = "harmony")
sce_merged <- FindClusters(sce_merged, resolution = 0.5)
# umap dimensionality reduction
sce_merged <- RunUMAP(sce_merged, dims = 1:30, reduction = "harmony")
DimPlot(sce_merged, reduction = "umap", label = TRUE)
qsave(sce_merged, file.path(intermediate_data.dir, "preprocessed_data.qs"))
# check batch effect
plot1 <- DimPlot(sce_merged, shuffle = T, reduction = "umap", label = F, group.by = "seurat_clusters") + 
  ggtitle("Clusters") + 
  labs(x = "UMAP1", y = "UMAP2") +
  scale_color_manual(values = adaptive_palette(28))
plot2 <- DimPlot(sce_merged, shuffle = T, reduction = "umap", group.by = "orig.ident") + 
  ggtitle("Samples") + 
  scale_color_manual(values = adaptive_palette(12)) +
  labs(x = "UMAP1", y = "UMAP2")
batch_effect_plot <- (plot1 + plot2) & my_sc_plot_theme() & theme(legend.key.size = unit(0.35, "cm"))
qsave(batch_effect_plot, paste0(results.dir, "/all_cells_batch_effect.qs"))
ggsave(filename = paste0(results.dir, "/all_cells_batch_effect.tiff"), width = 8, height = 3.5)
