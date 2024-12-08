library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyverse)
library(AUCell)
library(qs)
library(scales)

rm(list = ls())
intermediate_data.dir <- "results/step_1_hdwgcna_wgcna/scRNAseq_preprocessing/intermediate_data"
results.dir <- "results/step_1_hdwgcna_wgcna/scRNAseq_preprocessing/results/cell_annotation"
if (!dir.exists(results.dir)) dir.create(results.dir, recursive = T)
source("scripts/utils/plot_settings.R")
sce_merged <- qread("results/step_1_hdwgcna_wgcna/scRNAseq_preprocessing/intermediate_data/preprocessed_data.qs")

# round 1 cell annotation
all_cells_result_dir <- file.path(results.dir, "all_cells")
if (!dir.exists(all_cells_result_dir)) dir.create(all_cells_result_dir, recursive = T)
all_cells_annotation_clusters <- DotPlot(sce_merged, group.by = "seurat_clusters",
  features = c(
    "CD3D", "IL7R", "CD3E", "CD2", # T cells
    "NKG7", "GNLY", "KLRF1", "KLRD1", # NK cells
    "MS4A1", "CD79A", "IGKC", "MZB1", # B cells
    "C1QB", "CD14", "CSF1R", "CD68", # Myeloid cells
    "PDGFRB", "ACTA2", "TAGLN", "MYL9", # Fibroblasts
    "PECAM1", "ENG", "CDH5", "VWF", # Endothelial cells
    "EPCAM", "KRT8", "FXYD2", "KRT18" # Epithelial cells
  )
) +
  scale_size_continuous(limits = c(0,100),breaks = c(25, 50, 75),range = c(0,5)) +
  scale_color_gradient(
    low = "gold",
    high = "firebrick3",
    limits = c(-1, 1),
    breaks = c(-1,-0.5,0,0.5,1),
    oob = scales::squish
    )+
  guides(
    size = guide_legend(title = "Percent\nExpressed"),
    color = guide_colorbar(title = "Average\nExpression")
  )+
  labs(x = NULL, y = "Clusters") +
  my_plot_theme(x.text.angle = 90) +
  theme(legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(l = 0))
qsave(
  all_cells_annotation_clusters,
  paste0(all_cells_result_dir, "/all_cells_annotation_clusters.qs")
)
ggsave(file.path(all_cells_result_dir, "all_cells_annotation_clusters.pdf"), 
       width = 3.25, height = 3)
t_nk_cells <- 5
b_cells <- 22
myeloid_cells <- c(1, 6, 11, 18, 24, 25)
fibroblast_cells <- c(4, 12, 17, 20, 21, 26)
endothelial_cells <- c(2, 8, 9, 14, 16)
epithelial_cells <- c(0, 3, 7, 10, 13, 15, 19, 23, 27)
new.cluster.ids <- vector(length = length(table(sce_merged$seurat_clusters)))
new.cluster.ids[epithelial_cells + 1] <- "Epithelial Cells"
new.cluster.ids[t_nk_cells + 1] <- "T Cells & NK Cells"
new.cluster.ids[myeloid_cells + 1] <- "Myeloid Cells"
new.cluster.ids[endothelial_cells + 1] <- "Endothelial Cells"
new.cluster.ids[fibroblast_cells + 1] <- "Fibroblast Cells"
new.cluster.ids[b_cells + 1] <- "B Cells"
names(new.cluster.ids) <- c(0:(length(table(sce_merged$seurat_clusters)) - 1))
sce_merged <- RenameIdents(sce_merged, new.cluster.ids)
sce_merged$CellType <- Idents(sce_merged)
sce_merged$CellType <- factor(sce_merged$CellType, levels = rev(c("T Cells & NK Cells", "B Cells", "Myeloid Cells", "Fibroblast Cells", "Endothelial Cells", "Epithelial Cells")), ordered = T)
Idents(sce_merged) <- "CellType"
qsave(sce_merged, file.path(intermediate_data.dir, "all_cells_annotated.qs"))
all_cells_annotation_celltype <- DotPlot(sce_merged,
  cols = c("gold", "firebrick3"),
  features = c(
    "CD3D", "IL7R", "CD3E", "CD2", # T cells
    "NKG7", "GNLY", "KLRF1", "KLRD1", # NK cells
    "MS4A1", "CD79A", "IGKC", "MZB1", # B cells
    "C1QB", "CD14", "CSF1R", "CD68", # Myeloid cells
    "PDGFRB", "ACTA2", "TAGLN", "MYL9", # Fibroblasts
    "PECAM1", "ENG", "CDH5", "VWF", # Endothelial cells
    "EPCAM", "KRT8", "FXYD2", "KRT18" # Epithelial cells
  )
) + labs(x = NULL, y = NULL) +
  scale_size_continuous(limits = c(0,100),breaks = c(25, 50, 75),range = c(0,5))  +
  scale_color_gradient(
    low = "gold",
    high = "firebrick3",
    limits = c(-1, 1),
    breaks = c(-1,-0.5,0,0.5,1),
    oob = scales::squish
  )+
  guides(
    size = guide_legend(title = "Percent\nExpressed"),
    color = guide_colorbar(title = "Average\nExpression")
  ) +
  my_plot_theme(x.text.angle = 90) +
  theme(legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(l = 0))

qsave(all_cells_annotation_celltype, 
      file.path(all_cells_result_dir, "all_cells_annotation_celltype.qs"))
ggsave(filename = file.path(all_cells_result_dir, "all_cells_annotation_celltype.tiff"), 
       width = 3.25, height = 2.67)
all_cells_celltype_umap <- DimPlot(sce_merged)+
  scale_color_npg()+
  labs(x="UMAP1", y="UMAP2")+
  my_sc_plot_theme()+
  theme(legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(l = 0))
qsave(all_cells_celltype_umap, file.path(all_cells_result_dir, "all_cells_celltype_umap.qs"))
ggsave(file.path(all_cells_result_dir, "all_cells_celltype_umap.tiff"), 
       all_cells_celltype_umap, width = 3.25, height = 2.67)

# round 2 cell annotation
## round 2 cell preprocessing
epi_cells_result_dir <- file.path(results.dir, "epi_cells")
if (!dir.exists(epi_cells_result_dir)) dir.create(epi_cells_result_dir, recursive = T)
sce_epi <- subset(sce_merged, idents = "Epithelial Cells")
sce_epi <- NormalizeData(sce_epi, normalization.method = "LogNormalize", scale.factor = 1e4)
sce_epi <- FindVariableFeatures(sce_epi)
sce_epi <- ScaleData(sce_epi)
sce_epi <- RunPCA(sce_epi, features = VariableFeatures(object = sce_epi))
sce_epi <- RunHarmony(sce_epi, group.by.vars = "orig.ident", plot_convergence = T)
ElbowPlot(sce_epi, ndims = 50)
sce_epi <- FindNeighbors(sce_epi, dims = 1:30, reduction = "harmony")
sce_epi <- FindClusters(sce_epi, resolution = 0.5)
sce_epi <- RunUMAP(sce_epi, dims = 1:30, reduction = "harmony")
plot1 <- DimPlot(sce_epi, reduction = "umap", group.by = "seurat_clusters", label = F, shuffle = T) +
  scale_color_manual(values = adaptive_palette(18)) + 
  ggtitle("Clusters") +
  labs(x = "UMAP1", y = "UMAP2") +
  guides(color = guide_legend(ncol = 2, override.aes = list(size = 3)))+
  my_sc_plot_theme()
plot2 <- DimPlot(sce_epi, reduction = "umap", group.by = "orig.ident", shuffle = T) +
  scale_color_manual(values = adaptive_palette(12)) + 
  ggtitle("Samples")+
  labs(x = "UMAP1", y = "UMAP2") +
  my_sc_plot_theme()
epi_cells_batch_effect <- (plot1 | plot2) & theme(legend.key.height = unit(0.3, "cm"),
                                                  legend.key.width = unit(0.1, "cm"),
                                                  legend.margin = margin(l = 0))
ggsave(filename = paste0(epi_cells_result_dir, "/epi_cells_batch_effect.tiff"), width = 6.25, height = 2.67)
qsave(epi_cells_batch_effect, paste0(epi_cells_result_dir, "/epi_cells_batch_effect.qs"))
## round 2 annotation
epi_cells_annotation_clusters <- DotPlot(sce_epi,
  features = c(
    "ALDOB", "SLC22A8", "GLYAT", "MIOX", # Proximal Tubule Cells
    "CA9", "NDUFA4L2", "ANGPTL4", "NNMT", # RCC Cells
    "CALB1", "KCNJ1", "SLC8A1", "CLDN8", # Distal Tubule Cells
    "SLC12A1", "UMOD", "CLDN16", "CLDN10", # Loop of Henle Cells
    "FOXI1", "SLC4A1", "SLC26A4", "ATP6V0D2" # Collecting Tubule Cells
  ), group.by = "seurat_clusters"
)  + labs(x = NULL, y = NULL) +
  scale_size_continuous(limits = c(0,100),breaks = c(25, 50, 75),range = c(0,5)) +
  scale_color_gradient(
    low = "gold",
    high = "firebrick3",
    limits = c(-1, 1),
    breaks = c(-1,-0.5,0,0.5,1),
    oob = scales::squish
  )+
  guides(
    size = guide_legend(title = "Percent\nExpressed"),
    color = guide_colorbar(title = "Average\nExpression")
  ) +
  my_plot_theme(x.text.angle = 90) +
  theme(legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(l = 0))
qsave(
  epi_cells_annotation_clusters,
  paste0(epi_cells_result_dir, "/epi_cells_annotation_clusters.qs")
)
ggsave(
  filename = paste0(epi_cells_result_dir, "/epi_cells_annotation_clusters.tiff"),
  width = 3.25,
  height = 2.67
)
normal_malignat_proximal_tubule_cells <- c(0, 2, 4, 15, 16, 1, 3, 5, 6, 7, 9, 11, 12, 13, 17)
distal_tubule_cells <- c(10)
loop_of_henle_cells <- c(14)
collecting_tubule_cells <- c(8)
new.cluster.ids <- vector(length = length(table(sce_epi$seurat_clusters)))
new.cluster.ids[normal_malignat_proximal_tubule_cells + 1] <- "Normal Proximal Tubule Cells & Malignant Cells"
new.cluster.ids[loop_of_henle_cells + 1] <- "Loop of Henle Cells"
new.cluster.ids[distal_tubule_cells + 1] <- "Distal Tubule Cells"
new.cluster.ids[collecting_tubule_cells + 1] <- "Collecting Tubule Cells"
names(new.cluster.ids) <- c(0:(length(table(sce_epi$seurat_clusters)) - 1))
sce_epi <- RenameIdents(sce_epi, new.cluster.ids)
sce_epi$CellType <- Idents(sce_epi)
sce_epi$CellType <- factor(sce_epi$CellType,
  levels = rev(c(
    "Normal Proximal Tubule Cells & Malignant Cells",
    "Distal Tubule Cells",
    "Loop of Henle Cells",
    "Collecting Tubule Cells"
  )),
  ordered = T
)
Idents(sce_epi) <- "CellType"
epi_cells_annotation_celltypes <- DotPlot(sce_epi,
  features = c(
    "ALDOB", "SLC22A8", "GLYAT", "MIOX", # Proximal Tubule Cells
    "CA9", "NDUFA4L2", "ANGPTL4", "NNMT", # RCC Cells
    "CALB1", "KCNJ1", "SLC8A1", "CLDN8", # Distal Tubule Cells
    "SLC12A1", "UMOD", "CLDN16", "CLDN10", # Loop of Henle Cells
    "FOXI1", "SLC4A1", "SLC26A4", "ATP6V0D2" # Collecting Tubule Cells
  ), group.by = "CellType"
) +
  labs(x = NULL, y = NULL) +
  scale_y_discrete(labels = rev(c(
    "Normal Proximal Tubule Cells & Malignant Cells",
    "Distal Tubule Cells",
    "Loop of Henle Cells",
    "Collecting Tubule Cells"
  )) %>% str_wrap(width = 15)) +
  scale_size_continuous(limits = c(0,100),breaks = c(25, 50, 75),range = c(0,5)) +
  scale_color_gradient(
    low = "gold",
    high = "firebrick3",
    limits = c(-1, 1),
    breaks = c(-1,-0.5,0,0.5,1),
    oob = scales::squish
  )+
  guides(
    size = guide_legend(title = "Percent\nExpressed"),
    color = guide_colorbar(title = "Average\nExpression")
  ) +
  my_plot_theme(x.text.angle = 90) +
  theme(legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(l = 0))
qsave(
  epi_cells_annotation_celltypes ,
  file.path(epi_cells_result_dir, "epi_cells_annotation_celltypes.qs")
)
ggsave(
  filename = file.path(epi_cells_result_dir, "epi_cells_annotation_celltypes.tiff"),
  width = 3.25,
  height = 2.67
)
qsave(sce_epi, file.path(intermediate_data.dir, "epi_cells_annotated.qs"))

sce_epi_celltype_umap <- DimPlot(sce_epi, reduction = "umap",label = FALSE) +
  scale_color_manual(values = rev(pal_npg()(4)),
                     labels = str_wrap(levels(sce_epi$CellType),
                                       width = 25)) +
  labs(x="UMAP1", y="UMAP2")+
  my_sc_plot_theme()+
  theme(legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(l = 0))
qsave(sce_epi_celltype_umap, file.path(epi_cells_result_dir, "epi_cells_celltype_umap.qs"))
ggsave(file.path(epi_cells_result_dir, "epi_cells_celltype_umap.tiff"), 
       sce_epi_celltype_umap, width = 3.25, height = 2.67)


# round 3 annotation
mali_pro_result_dir <- file.path(results.dir, "mali_pro")
if (!dir.exists(mali_pro_result_dir)) dir.create(mali_pro_result_dir, recursive = T)
sce_mali_pro <- subset(sce_epi, idents = "Normal Proximal Tubule Cells & Malignant Cells")
sce_mali_pro <- NormalizeData(sce_mali_pro, normalization.method = "LogNormalize", scale.factor = 1e4)
sce_mali_pro <- FindVariableFeatures(sce_mali_pro)
sce_mali_pro <- ScaleData(sce_mali_pro)
sce_mali_pro <- RunPCA(sce_mali_pro, features = VariableFeatures(object = sce_mali_pro))
sce_mali_pro <- RunHarmony(sce_mali_pro, group.by.vars = "orig.ident", plot_convergence = T)
ElbowPlot(sce_mali_pro, ndims = 50)
sce_mali_pro <- FindNeighbors(sce_mali_pro, dims = 1:30, reduction = "harmony")
sce_mali_pro <- FindClusters(sce_mali_pro, resolution = 1)
sce_mali_pro <- RunUMAP(sce_mali_pro, dims = 1:30, reduction = "harmony")
plot1 <- DimPlot(sce_mali_pro, reduction = "umap", label = F, group.by = "seurat_clusters") +
  scale_color_manual(values = adaptive_palette(22)) +
  ggtitle("Clusters")+
  labs(x = "UMAP1", y = "UMAP2")+
  my_sc_plot_theme()

plot2 <- DimPlot(sce_mali_pro, reduction = "umap", group.by = "orig.ident", shuffle = T) +
  scale_color_manual(values = adaptive_palette(12)) +
  ggtitle("Samples")+
  labs(x = "UMAP1", y = "UMAP2")+
  my_sc_plot_theme()
malignant_pro_batch_effect <- (plot1 | plot2) & theme(legend.key.height = unit(0.3, "cm"),
                                                      legend.key.width = unit(0.1, "cm"),
                                                      legend.margin = margin(l = 0))
ggsave(filename = file.path(mali_pro_result_dir, "mali_pro_batch_effect.tiff"), width = 6.5, height = 2.67)
qsave(malignant_pro_batch_effect, file.path(mali_pro_result_dir, "mali_pro_batch_effect.qs"))
mali_pro_cells_annotation_clusters <- DotPlot(sce_mali_pro,
  features = c(
    "ALDOB", "SLC22A8", "GLYAT", "MIOX",
    "CA9", "NDUFA4L2", "ANGPTL4", "NNMT"
  ), group.by = "seurat_clusters"
) + labs(x = NULL, y = NULL) +
  scale_size_continuous(limits = c(0,100),breaks = c(25, 50, 75),range = c(0,5)) +
  scale_color_gradient(
    low = "gold",
    high = "firebrick3",
    limits = c(-1, 1),
    breaks = c(-1,-0.5,0,0.5,1),
    oob = scales::squish
  )+
  guides(
    size = guide_legend(title = "Percent\nExpressed"),
    color = guide_colorbar(title = "Average\nExpression")
  ) +
  my_plot_theme(x.text.angle = 90) +
  theme(legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(l = 0))
qsave(
  mali_pro_cells_annotation_clusters,
  file.path(mali_pro_result_dir, "mali_pro_cells_annotation_clusters.qs")
)
ggsave(
  filename = file.path(mali_pro_result_dir, "mali_pro_cells_annotation_clusters.tiff"),
  width = 3.25,
  height = 2.67
)
proximal_tubule_cells <- c(0, 1, 2, 6, 7, 8, 15, 17, 20)
malignant_cells <- c(3, 4, 5, 6, 9, 12, 14)
transitional_cells <- c(6, 10, 11, 13, 16, 18, 19, 21)
new.cluster.ids <- vector(length = length(table(sce_mali_pro$seurat_clusters)))
new.cluster.ids[malignant_cells + 1] <- "Malignant Cells"
new.cluster.ids[proximal_tubule_cells + 1] <- "Proximal Tubule Cells"
new.cluster.ids[transitional_cells + 1] <- "Transitional Cells"
names(new.cluster.ids) <- c(0:(length(table(sce_mali_pro$seurat_clusters)) - 1))
sce_mali_pro <- RenameIdents(sce_mali_pro, new.cluster.ids)
sce_mali_pro$CellType <- Idents(sce_mali_pro)
sce_mali_pro$CellType <- factor(sce_mali_pro$CellType, levels = c("Malignant Cells", "Transitional Cells", "Proximal Tubule Cells"), ordered = T)
mali_pro_cells_annotation_celltypes <- DotPlot(sce_mali_pro,
  features = c(
    "ALDOB", "SLC22A8", "GLYAT", "MIOX",
    "CA9", "NDUFA4L2", "ANGPTL4", "NNMT"
  ), group.by = "CellType",
  cols = c("gold", "firebrick3")
) + labs(x = NULL, y = NULL) +
  scale_size_continuous(limits = c(0,100),breaks = c(25, 50, 75),range = c(0,5)) +
  scale_color_gradient(
    low = "gold",
    high = "firebrick3",
    limits = c(-1, 1),
    breaks = c(-1,-0.5,0,0.5,1),
    oob = scales::squish
  )+
  guides(
    size = guide_legend(title = "Percent\nExpressed"),
    color = guide_colorbar(title = "Average\nExpression")
  ) +
  my_plot_theme(x.text.angle = 90) +
  theme(legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(l = 0))
qsave(
  mali_pro_cells_annotation_celltypes,
  file.path(mali_pro_result_dir, "mali_pro_cells_annotation_celltypes.qs")
)
ggsave(
  filename = file.path(mali_pro_result_dir, "mali_pro_cells_annotation_celltypes.tiff"),
  width = 3.25,
  height = 2.67
)
qsave(sce_mali_pro, file.path(intermediate_data.dir, "mali_pro_cells_annotated.qs"))

sce_mali_pro$CellType <- factor(sce_mali_pro$CellType, levels = c("Proximal Tubule Cells", "Transitional Cells", "Malignant Cells"), ordered = T)
mali_pro_celltype_umap <- DimPlot(sce_mali_pro,reduction = "umap",label = FALSE) +
  scale_color_manual(values = c("Malignant Cells"="#E64B35FF",
                                "Transitional Cells"="#4DBBD5FF",
                                "Proximal Tubule Cells"="#00A087FF"))+
  labs(x="UMAP1", y="UMAP2")+
  my_sc_plot_theme()+
  theme(legend.key.height = unit(0.3, "cm"),
        legend.margin = margin(l = 0))
qsave(mali_pro_celltype_umap, file.path(mali_pro_result_dir, "mali_pro_celltype_umap.qs"))
ggsave(file.path(mali_pro_result_dir, "mali_pro_celltype_umap.tiff"), 
       mali_pro_celltype_umap, width = 3.25, height = 2.67)






