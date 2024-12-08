library(Seurat)
library(ggplot2)
library(dplyr)
library(qs)
library(monocle)
rm(list = ls())
source("scripts/utils/plot_settings.R")

mali_pro_cells <- qread("results/step_1_hdwgcna_wgcna/pathway_enrichment_score/all_cells/all_cells_auc_calculated.qs")
mali_pro_exp <- mali_pro_cells@assays$RNA@counts
p_data <- mali_pro_cells@meta.data
p_data <- select(p_data, !contains("DF.classifications"))
p_data <- select(p_data, !contains("pANN"))
f_data <- data.frame(gene_short_name = row.names(mali_pro_cells), row.names = row.names(mali_pro_cells))
pd <- new("AnnotatedDataFrame", data = p_data)
fd <- new("AnnotatedDataFrame", data = f_data)
cds <- newCellDataSet(mali_pro_exp, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10)) 
print(head(fData(cds)))
diff <- differentialGeneTest(cds[expressed_genes, ], fullModelFormulaStr = "~CellType", cores = 6) 
deg <- subset(diff, qval < 0.001)
deg <- arrange(deg, qval)
deg <- deg[1:2000, ]
head(deg)
write.table(deg, file = "results/step_1_hdwgcna_wgcna/trajectory_inference/train.monocle.DEG.txt", col.names = T, row.names = T, sep = "\t", quote = F)
ordergene <- rownames(deg)
cds <- setOrderingFilter(cds, ordergene)
plot_ordering_genes(cds)
cds <- reduceDimension(cds, max_components = 2, method = "DDRTree")
cds <- orderCells(cds)
cds$CellType <- factor(cds$CellType, levels = c("Malignant Cells", "Transitional Cells", "Proximal Tubule Cells"), ordered = T)
qsave(cds, "results/step_1_hdwgcna_wgcna/trajectory_inference/cds.qs")

pseudotime_plot <- plot_cell_trajectory(cds = cds, color_by = "Pseudotime", size = 1, show_backbone = T) +
  viridis::scale_color_viridis(option = "D") +
  guides(
    color = guide_colorbar(barwidth = 2.5, barheight = 0.5, direction = "horizontal")
  )+
  my_sc_plot_theme()+
  theme(legend.position = c(0.7,0.95))+
  ggtitle("Pseudotime Trajectory")
celltype_plot <- plot_cell_trajectory(cds, color_by = "CellType", size = 0.1, show_backbone = T) +
  scale_color_manual(values = c("Malignant Cells"="#E64B35FF",
                                "Transitional Cells"="#4DBBD5FF",
                                "Proximal Tubule Cells"="#00A087FF")) +
  labs(color = "") + 
  my_sc_plot_theme() +
  theme(legend.position =c(0.7,1),
        legend.key.height  = unit(0.3, "cm")) +
  ggtitle("Cell Type")
tryptophan_metabolism_plot <- plot_cell_trajectory(cds, color_by = "AUC", cell_size = 1) +
  viridis::scale_color_viridis(option = "B", begin = 0, end = 1) +
  labs(color = "AUC")+
  guides(
    color = guide_colorbar(barwidth = 2.5, barheight = 0.5, direction = "horizontal")
  )+
  my_sc_plot_theme()+
  theme(legend.position =c(0.7,0.95)) +
  ggtitle("Tryptophan Metabolism")

qsave(pseudotime_plot, "results/step_1_hdwgcna_wgcna/trajectory_inference/trajectory_inference.qs")
qsave(celltype_plot, "results/step_1_hdwgcna_wgcna/trajectory_inference/celltype.qs")
qsave(tryptophan_metabolism_plot, "results/step_1_hdwgcna_wgcna/trajectory_inference/tryptophan_metabolism.qs")



