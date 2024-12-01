library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(AUCell)

rm(list = ls())

source("scripts/utils/plot_settings.R")
source("scripts/utils/scRNAseq_utils.R")
all_cells <- qread("results/step_1_hdwgcna_wgcna/scRNAseq_preprocessing/intermediate_data/all_cells_annotated.qs")
epi_cells <- qread("results/step_1_hdwgcna_wgcna/scRNAseq_preprocessing/intermediate_data/epi_cells_annotated.qs")
mali_pro_cells <- qread("results/step_1_hdwgcna_wgcna/scRNAseq_preprocessing/intermediate_data/mali_pro_cells_annotated.qs")
results.dir <- "results/step_1_hdwgcna_wgcna/pathway_enrichment_score"
tryptophan_metabolism_geneset <- read.delim("data/genesets/tryptophan_metabolism.txt", header = F)
genesets <- list()
genesets[["Tryptophan_Metabolism"]] <- tryptophan_metabolism_geneset[,1]
# Calculate AUCell scores
all_cells$AUC <- cal_pathway_auc(all_cells, genesets)[1,]
epi_cells$AUC <- cal_pathway_auc(epi_cells, genesets)[1,]
mali_pro_cells$AUC <- cal_pathway_auc(mali_pro_cells, genesets)[1,]

# draw figures
# all cells
all_cells_results_dir <- file.path(results.dir, "all_cells")

all_cells_auc_umap <- FeaturePlot(all_cells, features = "AUC") + 
  ggtitle("Tryptophan Metabolism") +
  viridis::scale_color_viridis(option = "A") +
  labs(x = "UMAP1", y = "UMAP2",color = "AUC") +
  my_sc_plot_theme()+
  theme(legend.key.height = unit(0.3, "cm"))

qsave(all_cells_auc_umap, file.path(all_cells_results_dir, "all_cells_auc_umap.qs"))
ggsave(file.path(all_cells_results_dir, "all_cells_auc_umap.tiff"), all_cells_auc_umap, width = 3.25, height = 2.67)

all_cells_auc_violin<- ggplot(all_cells@meta.data, aes(x = CellType, y = AUC, fill = CellType)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  scale_x_discrete(labels = str_wrap(levels(all_cells$CellType),width = 10)) +
  scale_fill_npg(alpha = 0.8) +
  labs(y = "Tryptophan Metabolism (AUC)", x = "", title = "") +
  my_plot_theme(x.text.angle = 0)+
  theme(legend.position = "none")

qsave(all_cells_auc_violin, file.path(all_cells_results_dir, "all_cells_auc_violin.qs"))
ggsave(file.path(all_cells_results_dir, "all_cells_auc_violin.tiff"), all_cells_auc_violin, width = 3.25, height = 2.67)

qsave(all_cells, file.path(all_cells_results_dir, "all_cells_auc_calculated.qs"))
# epi cells
epi_cells_results_dir <- file.path(results.dir, "epi_cells")
epi_cells_auc_umap <- FeaturePlot(epi_cells, features = "AUC") + 
  ggtitle("Tryptophan Metabolism") +
  viridis::scale_color_viridis(option = "A") +
  labs(x = "UMAP1", y = "UMAP2",color = "AUC") +
  my_sc_plot_theme()+
  theme(legend.key.height = unit(0.3, "cm"))
qsave(epi_cells_auc_umap, file.path(epi_cells_results_dir, "epi_cells_auc_umap.qs"))
ggsave(file.path(epi_cells_results_dir, "epi_cells_auc_umap.tiff"), epi_cells_auc_umap, width = 3.25, height = 2.67)
      
epi_cells_auc_violin <- ggplot(epi_cells@meta.data, aes(x = CellType, y = AUC, fill = CellType)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  my_plot_theme(x.text.angle = 45) +
  scale_x_discrete(labels = str_wrap(levels(epi_cells$CellType),width = 15)) +
  scale_fill_npg(alpha = 0.8) +
  labs(y = "Tryptophan Metabolism (AUC)", x = "", title = "") +
  my_plot_theme(x.text.angle = 0)+
  theme(legend.position = "none")
qsave(epi_cells_auc_violin, file.path(epi_cells_results_dir, "epi_cells_auc_violin.qs"))
ggsave(file.path(epi_cells_results_dir, "epi_cells_auc_violin.tiff"), epi_cells_auc_violin, width = 3.25, height = 2.67)

qsave(epi_cells, file.path(epi_cells_results_dir, "epi_cells_auc_calculated.qs"))

# mali_pro cells
malig_pro_results_dir <- file.path(results.dir, "mali_pro_cells")
mali_pro_cells_auc_umap <- FeaturePlot(mali_pro_cells, features = "AUC") + 
  ggtitle("Tryptophan Metabolism") +
  viridis::scale_color_viridis(option = "A") +
  labs(x = "UMAP1", y = "UMAP2",color = "AUC") +
  my_sc_plot_theme() +
  theme(legend.key.height = unit(0.3, "cm"))
qsave(mali_pro_cells_auc_umap, file.path(malig_pro_results_dir, "mali_pro_cells_auc_umap.qs"))
ggsave(file.path(malig_pro_results_dir, "mali_pro_cells_auc_umap.tiff"), mali_pro_cells_auc_umap, width = 3.25, height = 2.67)

comparasions <- list(c("Malignant Cells", "Transitional Cells"), c("Transitional Cells", "Proximal Tubule Cells"), c("Malignant Cells", "Proximal Tubule Cells"))

mali_pro_cells_auc_violin <- ggplot(mali_pro_cells@meta.data, aes(x = CellType, y = AUC, fill = CellType)) +
  geom_violin() +
  geom_boxplot(width = 0.1) +
  my_plot_theme(x.text.angle = 45) +
  scale_x_discrete(labels = str_wrap(levels(mali_pro_cells$CellType),width = 10)) +
  scale_fill_npg(alpha = 0.8) +
  scale_y_continuous(limits = c(0, max(mali_pro_cells$AUC) + 0.2)) +
  stat_compare_means(method = "anova", label.x.npc = 0.38, label.y = max(mali_pro_cells$AUC) + 0.15,size = 2.5) +
  stat_compare_means(comparisons = comparasions, method = "t.test", hide.ns = F, label = "p.signif")+
  labs(y = "Tryptophan Metabolism (AUC)", x = "", title = "") +
  my_plot_theme(x.text.angle = 0)+
  theme(legend.position = "none")
qsave(mali_pro_cells_auc_violin, file.path(malig_pro_results_dir, "mali_pro_cells_auc_violin.qs"))
ggsave(file.path(malig_pro_results_dir, "mali_pro_cells_auc_violin.tiff"), mali_pro_cells_auc_violin, width = 6, height = 5)

qsave(mali_pro_cells, file.path(malig_pro_results_dir, "mali_pro_cells_auc_calculated.qs"))
