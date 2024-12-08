library(Seurat)
library(hdWGCNA)
library(ggplotify)
library(qs)
library(gridGraphics)
library(grid)
library(extrafont)
library(magick)
library(gtools)
rm(list = ls())
source("scripts/utils/plot_settings.R")
sce_mali_pro <- qread("results/step_1_hdwgcna_wgcna/pathway_enrichment_score/mali_pro_cells/mali_pro_cells_auc_calculated.qs")
results_dir <- "results/step_1_hdwgcna_wgcna/hdWGCNA"
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}
# Set up Seurat object for hdWGCNA
sce_mali_pro_WGCNA <- SetupForWGCNA(
  sce_mali_pro,
  gene_select = "fraction",
  fraction = 0.03,
  wgcna_name = "sce_mali_pro_WGCNA"
)
# Construct metacells
sce_mali_pro_WGCNA <- MetacellsByGroups(
  seurat_obj = sce_mali_pro_WGCNA,
  group.by = c("CellType", "orig.ident"),
  reduction = "harmony",
  k = 25,
  max_shared = 10,
  ident.group = "CellType"
)
# Normalize metacell expression matrix
sce_mali_pro_WGCNA <- NormalizeMetacells(sce_mali_pro_WGCNA)
# Determine the gene expression matrix for hdWGCNA
sce_mali_pro_WGCNA <- SetDatExpr(
  group_name = names(table(sce_mali_pro_WGCNA$CellType)),
  sce_mali_pro_WGCNA,
  group.by = "CellType",
  assay = "RNA",
  slot = "data"
)
# Select soft-power threshold
sce_mali_pro_WGCNA <- TestSoftPowers(
  sce_mali_pro_WGCNA,
  networkType = "signed"
)
power_table <- GetPowerTable(sce_mali_pro_WGCNA)
selected_soft_power <- power_table$Power[power_table$SFT.R.sq >= 0.85][1]
write.csv(power_table, file = paste0(results_dir, "/Soft_Power_Table.csv"), row.names = F)
soft_power_selection_fig_list <- PlotSoftPowers(sce_mali_pro_WGCNA, 
                                                selected_power = selected_soft_power,
                                                text_size = 2)
soft_power_selection_fig <- (soft_power_selection_fig_list[[1]] +
  soft_power_selection_fig_list[[2]]) &
  my_plot_theme()
qsave(soft_power_selection_fig, file = paste0(results_dir, "/Soft_Power_Selection_Fig.qs"))
ggsave(
  paste0(results_dir, "/Soft_Power_Selection_Fig.tiff"),
  soft_power_selection_fig,
  width = 6.5,
  height = 2.67
)
  
# Construct co-expression network
sce_mali_pro_WGCNA <- ConstructNetwork(
  sce_mali_pro_WGCNA,
  soft_power = selected_soft_power,
  minModuleSize = 50,
  mergeCutHeight = 0.2,
  setDatExpr = FALSE,
  tom_outdir = results_dir,
  overwrite_tom = T,
  tom_name = "sce_mali_pro_WGCNA_tom"
)
font_import(pattern = "Arial", prompt = FALSE) 
loadfonts(device = "pdf") 
pdf(file.path(results_dir, "hdWGCNA_Dendrogram.pdf"),width = 3.25, height = 2.67,
    pointsize = 7, family = "Arial")
PlotDendrogram(sce_mali_pro_WGCNA, main = "hdWGCNA Dendrogram")
dev.off()
hdwgcna_dendrogram <- image_read(file.path(results_dir, "hdWGCNA_Dendrogram.pdf"), density = 600)
grob_image <- rasterGrob(hdwgcna_dendrogram, interpolate = TRUE) 
hdwgcna_dendrogram <- ggplot() +
  annotation_custom(grob_image, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
  theme_void()
qsave(hdwgcna_dendrogram, file = paste0(results_dir, "/hdWGCNA_Dendrogram.qs"))
ggsave(
  paste0(results_dir, "/hdWGCNA_Dendrogram.tiff"),
  hdwgcna_dendrogram,
  width = 4,
  height = 3.5
)
# Module eigengenes and connectivity
sce_mali_pro_WGCNA <- ScaleData(sce_mali_pro_WGCNA, features = VariableFeatures(sce_mali_pro_WGCNA))
sce_mali_pro_WGCNA <- ModuleEigengenes(
  sce_mali_pro_WGCNA,
  group.by.vars = "orig.ident"
)
hMEs <- GetMEs(sce_mali_pro_WGCNA)
MEs <- GetMEs(sce_mali_pro_WGCNA, harmonized = FALSE)
# compute eigengene-based connectivity (kME):
sce_mali_pro_WGCNA <- ModuleConnectivity(
  sce_mali_pro_WGCNA,
  group.by = "CellType",
  group_name = names(table(sce_mali_pro_WGCNA$CellType)),
  corOptions = "use='p',method = 'spearman'"
)
qsave(sce_mali_pro_WGCNA, file = file.path(results_dir, "/sce_mali_pro_WGCNA.qs"))
# get the module assignment table
modules <- GetModules(sce_mali_pro_WGCNA)
write.csv(modules, file = file.path(results_dir, "module_information.csv"), row.names = F)
# get the module-trait relationship
traitData <- data.frame(KEGG_Tryptophan_Metabolism = sce_mali_pro_WGCNA$AUC)
moduleTraitCor <- cor(hMEs, traitData, method = "spearman")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, ncol(sce_mali_pro_WGCNA))
# plot the correlation heatmap
textMatrix <- paste("R = ",signif(moduleTraitCor, 2), "\nP=",
                    signif(moduleTraitPvalue, 1),
                    sep = ""
)
dim(textMatrix) <- dim(moduleTraitCor)
pdf(file.path(results_dir, "hdWGCNA_Traits_Correlation.pdf"),width = 1.625, height = 2,
    pointsize = 6, family = "Arial")
par(mar = c(2, 7, 2, 0), family = "Arial",font.main=2)
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = "Tryptophan Metabolism",
  yLabels = paste0("ME", colnames(hMEs)),
  ySymbols = paste0("ME", colnames(hMEs)),
  xLabelsAngle = 0,
  xLabelsAdj = c(0.5),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 1,
  zlim = c(-1, 1)
)
title("hdWGCNA",adj=0.35,cex = 1.2)
dev.off()
hdwgcna_traits_cor <- image_read(file.path(results_dir, "hdWGCNA_Traits_Correlation.pdf"), density = 600)
grob_image <- rasterGrob(hdwgcna_traits_cor, interpolate = TRUE) 
hdwgcna_traits_cor <- ggplot() +
  annotation_custom(grob_image, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
  theme_void()
qsave(hdwgcna_traits_cor, file = paste0(results_dir, "/hdWGCNA_Traits_Correlation.qs"))
ggsave(
  paste0(results_dir, "/hdWGCNA_Traits_Correlation.tiff"),
  hdwgcna_traits_cor,
  width = 3.25,
  height = 2.67
)

# calculate the correlation between gene exp and module
selected_trait <- as.data.frame(traitData$KEGG_Tryptophan_Metabolism)
datExpr <- as.matrix(sce_mali_pro_WGCNA@assays$RNA@data)[row.names(modules), row.names(hMEs)]
datExpr <- t(datExpr)
names(selected_trait) <- "KEGG_Tryptophan_Metabolism"
moduleNames <- names(hMEs)
geneModuleMembership <- as.data.frame(cor(datExpr, hMEs, method = "spearman"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), ncol(sce_mali_pro_WGCNA))) 
names(geneModuleMembership) <- paste("MM", moduleNames, sep = "")
names(MMPvalue) <- paste("p.MM", moduleNames, sep = "")
geneTraitSignificance <- as.data.frame(cor(datExpr, selected_trait, method = "spearman")) 
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), ncol(sce_mali_pro_WGCNA)))
names(geneTraitSignificance) <- paste("GS.", names(selected_trait), sep = "")
names(GSPvalue) <- paste("p.GS.", names(selected_trait), sep = "")
module <- "turquoise"
selected_trait_name <- "Tryptophan Metabolism"
column <- match(module, moduleNames)
moduleGenes <- rownames(modules)[modules$color == module]

gene_module_correlation <- data.frame(
  module_membership = abs(geneModuleMembership[moduleGenes, column]),
  gene_significance = abs(geneTraitSignificance[moduleGenes, 1])
)

mm_gs<- ggscatter(gene_module_correlation, x = "module_membership", y = "gene_significance",
          color = module,
          size = 1.5,
          cor.coef.size= 3,
          conf.int = TRUE,                 
          cor.coef = TRUE,                 
          cor.method = "pearson",          
          xlab = "Module Membership", 
          ylab = "Gene Significance") +
  ggtitle("hdWGCNA")+
  geom_vline(xintercept = quantile(abs(geneModuleMembership[moduleGenes, column]), 0.5), color = "black", linetype = "dashed")+
  geom_hline(yintercept = quantile(abs(geneTraitSignificance[moduleGenes, 1]), 0.5), color = "black", linetype = "dashed")+
  my_plot_theme()
qsave(mm_gs, file = paste0(results_dir, "/Module_Membership_Gene_Significance.qs"))
ggsave(
  paste0(results_dir, "/Module_Membership_Gene_Significance.tiff"),
  mm_gs,
  width = 3.25,
  height = 2.67
)
mm_gs$layers[[1]]$aes_params$size <- 1
# save the results
geneInfo0 <- base::merge(modules, geneTraitSignificance, by = "row.names", all = T)[, -1]
rownames(geneInfo0) <- geneInfo0$gene_name
geneInfo1 <- base::merge(geneInfo0, GSPvalue, by = "row.names", all = T)[, -1]
rownames(geneInfo1) <- geneInfo1$gene_name
geneInfo <- base::merge(geneInfo1, geneModuleMembership, by = "row.names", all = T)[, -1]
turquoise_gene <- geneInfo[geneInfo$color == "turquoise", ]
selected_gene <- turquoise_gene[abs(turquoise_gene$GS.KEGG_Tryptophan_Metabolism) >= quantile(abs(turquoise_gene$GS.KEGG_Tryptophan_Metabolism), 0.5) &
                                  abs(turquoise_gene$MMturquoise) >= quantile(abs(turquoise_gene$MMturquoise), 0.5), ]
write.csv(geneInfo, file = file.path(results_dir, "gene_information.csv"), row.names = F, quote = F)
write.csv(selected_gene, file = file.path(results_dir, "selected_gene_information.csv"), row.names = F, quote = F)
qsave(list(geneModuleMembership, geneTraitSignificance), file = file.path(results_dir, "mm_gs_info.qs"))

