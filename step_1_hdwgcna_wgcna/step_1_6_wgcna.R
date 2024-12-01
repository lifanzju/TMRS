library(WGCNA)
library(GSVA)
library(dplyr)
library(qs)
library(extrafont)
library(magick)
library(gridGraphics)
library(grid)
library(ggplotify)

rm(list = ls())
data <- read.delim("data/TCGA/TCGA.txt", row.names = 1)
protein_coding_gene <- read.delim(file = "data/protein_coding_gene.txt", header = F, row.names = NULL)[, 1]
genelist <- read.delim(file = "data/genesets/tryptophan_metabolism.txt", header = F, row.names = NULL)[, 1]
results_dir <- "results/step_1_hdwgcna_wgcna/WGCNA"
source("scripts/utils/plot_settings.R")
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

data <- data[protein_coding_gene, ]
data <- data.matrix(data)

# ssGSEA
ssgsea_results <- gsva(data, list(genelist), mx.diff = FALSE, verbose = FALSE, method = "ssgsea") %>%
  t() %>%
  data.frame()
colnames(ssgsea_results) <- "tryptophan_metabolism"
write.csv(ssgsea_results, file.path(results_dir, "ssgsea_results.csv"))

# WGCNA
datExpr <- as.data.frame(t(data))
# quality control
gsg <- goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) {
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ", ")))
  }
  if (sum(!gsg$goodSamples) > 0) {
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ", ")))
  }
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}
rm(data)
gc()
# prepare trait data
traitData <- ssgsea_results
traitData <- traitData[match(row.names(datExpr), row.names(traitData)), ]
dim(traitData) <- dim(ssgsea_results)
colnames(traitData) <- colnames(ssgsea_results)
row.names(traitData) <- row.names(datExpr)
traitData <- as.data.frame(traitData)
# formal analysis
options(stringsAsFactors = FALSE)
powers <- c(c(1:10), seq(from = 12, to = 30, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, blockSize = ncol(datExpr))
cex1 <- 0.85
soft_power <- sft$powerEstimate
pt <- sft$fitIndices
pt$text_color <- ifelse(pt$Power == soft_power, "white", "black")
sft_r <- pt$SFT.R.sq[pt$Power == soft_power]
scale_free_plot <- ggplot(data = pt, aes(x = Power, y = SFT.R.sq)) +
  geom_rect(
    aes(
      xmin = -Inf, xmax = Inf,
      ymin = -Inf, ymax = 0.8
    ),
    fill = "grey80", alpha = 0.8,
    color = NA
  ) +
  geom_hline(
    yintercept = sft_r,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = soft_power,
    linetype = "dashed"
  ) +
  geom_point(
    data = pt[pt$Power ==
      soft_power, c("Power", "SFT.R.sq")], aes(
      x = Power,
      y = SFT.R.sq
    ), inherit.aes = FALSE, color = "black",
    size = 5
  ) +
  geom_text(
    label = pt$Power,
    color = pt$text_color,
    size = 2
  ) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  ylab("Scale-free Topology Model Fit") +
  xlab("Soft Power Threshold") +
  my_plot_theme()

mean_k <- as.numeric(pt[pt$Power == soft_power, "mean.k."])
median_connectivity_plot <- ggplot(pt, aes(x = Power, y = mean.k.)) +
  geom_hline(yintercept = mean_k, linetype = "dashed") +
  geom_vline(xintercept = soft_power, linetype = "dashed") +
  geom_point(
    data = pt[
      pt$Power == soft_power,
      c("Power", "mean.k.")
    ], aes(x = Power, y = mean.k.),
    inherit.aes = FALSE, color = "black", size = 5
  ) +
  geom_text(
    label = pt$Power, color = pt$text_color,
    size = 2
  ) +
  scale_y_continuous(labels = scales::comma) +
  ylab("Mean Connectivity") +
  xlab("Soft Power Threshold") +
  my_plot_theme()




soft_power_selection_fig <- scale_free_plot + median_connectivity_plot
qsave(soft_power_selection_fig, file = paste0(results_dir, "/Soft_Power_Selection_Fig.qs"))
ggsave(
  paste0(results_dir, "/Soft_Power_Selection_Fig.tiff"),
  soft_power_selection_fig,
  width = 6.5,
  height = 2.67
)

qsave(scale_free_plot, file = paste0(results_dir, "/Scale_Free_Plot.qs"))
qsave(median_connectivity_plot, file = paste0(results_dir, "/Median_Connectivity_Plot.qs"))
power_selected <- sft$powerEstimate
adjacency <- adjacency(datExpr, power = power_selected)
TOM <- TOMsimilarity(adjacency)
rm(adjacency)
gc()
dissTOM <- 1 - TOM
qsave(TOM, file = file.path(results_dir, "TOM.qs"))
rm(TOM)
gc()
geneTree <- hclust(as.dist(dissTOM), method = "average")
plot(geneTree,
  xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity",
  labels = FALSE, hang = 0.04
)
minModuleSize <- 100
dynamicMods <- cutreeDynamic(
  dendro = geneTree, distM = dissTOM,
  deepSplit = 2, pamRespectsDendro = FALSE,
  minClusterSize = minModuleSize
)
table(dynamicMods)
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
sizeGrWindow(8, 6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05,
  main = "Gene dendrogram and module colors"
)
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss <- 1 - cor(MEs)
# Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average")
# Plot the result
MEDissThres <- 0.25
# Call an automatic merging function
merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors <- merge$colors
# Eigengenes of the new merged modules
mergedMEs <- merge$newMEs
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
  c("Dynamic Tree Cut", "Merged dynamic"),
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05
)
font_import(pattern = "Arial", prompt = FALSE)
loadfonts(device = "pdf")
pdf(file.path(results_dir, "WGCNA_Dendrogram.pdf"), width = 3.25, height = 2.67,
    pointsize = 7, family = "Arial")
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
  c("Dynamic\nTree Cut", "Merged dynamic"),
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05,
  main = "WGCNA Dendrogram"
)
dev.off()
wgcna_dendrogram <- image_read(file.path(results_dir, "WGCNA_Dendrogram.pdf"), density = 600)
grob_image <- rasterGrob(wgcna_dendrogram, interpolate = TRUE)
wgcna_dendrogram <- ggplot() +
  annotation_custom(grob_image, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
  theme_void()
qsave(wgcna_dendrogram, file = paste0(results_dir, "/WGCNA_Dendrogram.qs"))
ggsave(
  file.path(results_dir, "WGCNA_Dendrogram.tiff"),
  wgcna_dendrogram,
  width = 3.25,
  height = 2.67
)

# Rename to moduleColors
moduleColors <- mergedColors
# Construct numerical labels corresponding to the colors
colorOrder <- c("grey", standardColors(50))
moduleLabels <- match(moduleColors, colorOrder) - 1
MEs <- mergedMEs
# Save module colors and labels for use in subsequent parts
qsave(list(MEs, moduleLabels, moduleColors, geneTree), file = file.path(results_dir, "WGCNA.qs"))
# correlation between module eigengenes and traits
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
moduleTraitCor <- cor(MEs, traitData, method = "spearman")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
textMatrix <- paste("R = ", signif(moduleTraitCor, 2), ", P=",
  signif(moduleTraitPvalue, 1),
  sep = ""
)
dim(textMatrix) <- dim(moduleTraitCor)
pdf(file.path(results_dir, "WGCNA_Traits_Correlation.pdf"),width = 2.17, height = 4,
    pointsize = 6, family = "Arial")
par(mar = c(2, 8.5, 1, 0), family = "Arial",font.main=2)
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = "Tryptophan Metabolism",
  yLabels = names(MEs),
  ySymbols = names(MEs),
  xLabelsAngle = 0,
  xLabelsAdj = c(0.5),
  colorLabels = TRUE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 1,
  zlim = c(-1, 1),
)
title("WGCNA",adj=0.35,cex = 1.2)
dev.off()
wgcna_traits_cor <- image_read(file.path(results_dir, "WGCNA_Traits_Correlation.pdf"), density = 600)
grob_image <- rasterGrob(wgcna_traits_cor, interpolate = TRUE) 
wgcna_traits_cor <- ggplot() +
  annotation_custom(grob_image, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf)+
  theme_void()
qsave(wgcna_traits_cor, file = file.path(results_dir, "WGCNA_Traits_Correlation.qs"))
ggsave(
  file.path(results_dir, "WGCNA_Traits_Correlation.tiff"),
  wgcna_traits_cor,
  width = 3.25,
  height = 2.67*2
)
# calculate the correlation between gene exp and module
selected_trait <- as.data.frame(traitData$tryptophan_metabolism)
names(selected_trait) <- "Tryptophan Metabolism"
moduleNames <- substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, method = "spearman"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) <- paste("MM", moduleNames, sep = "")
names(MMPvalue) <- paste("p.MM", moduleNames, sep = "")
geneTraitSignificance <- as.data.frame(cor(datExpr, selected_trait, method = "spearman"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) <- paste("GS.", names(selected_trait), sep = "")
names(GSPvalue) <- paste("p.GS.", names(selected_trait), sep = "")
selected_trait_name <- "Tryptophan Metabolism"
module_list <- unique(moduleColors)
module <- "purple"
column <- match(module, moduleNames)
moduleGenes <- moduleColors == module
gene_module_correlation <- data.frame(
  module_membership = abs(geneModuleMembership[moduleGenes, column]),
  gene_significance = abs(geneTraitSignificance[moduleGenes, 1])
)
mm_gs <- ggscatter(gene_module_correlation,
  x = "module_membership", y = "gene_significance",
  color = module,
  size = 1.5,
  cor.coef.size= 3,
  conf.int = TRUE,
  cor.coef = TRUE,
  cor.method = "pearson",
  xlab = "Module Membership",
  ylab = "Gene Significance") +
  ggtitle("hdWGCNA")+
  geom_vline(xintercept = quantile(abs(geneModuleMembership[moduleGenes, column]), 0.5), color = "black", linetype = "dashed") +
  geom_hline(yintercept = quantile(abs(geneTraitSignificance[moduleGenes, 1]), 0.5), color = "black", linetype = "dashed") +
  my_plot_theme()
qsave(mm_gs, file = paste0(results_dir, "/Module_Membership_Gene_Significance.qs"))
ggsave(
  paste0(results_dir, "/Module_Membership_Gene_Significance.tiff"),
  mm_gs,
  width = 3.25,
  height = 2.67
)

# save the results
geneInfo0 <- data.frame(
  Gene = colnames(datExpr),
  moduleColor = moduleColors,
  geneTraitSignificance,
  GSPvalue
)
geneInfo <- base::merge(geneInfo0, geneModuleMembership, by = "row.names", all = T)[, -1]
all_gene_dir <- file.path(results_dir, "GeneOfEachModules/ALL")
selected_gene_dir <- file.path(results_dir, "GeneOfEachModules/Selected")
if (!dir.exists(all_gene_dir)) {
  dir.create(all_gene_dir, recursive = T)
}
if (!dir.exists(selected_gene_dir)) {
  dir.create(selected_gene_dir, recursive = T)
}
for (module in module_list) {
  genes <- geneInfo[geneInfo$moduleColor == module, ]
  module_name <- paste0("MM", module)
  selected_gene <- genes[abs(genes$GS.Tryptophan.Metabolism) >= median(abs(genes$GS.Tryptophan.Metabolism)) &
    abs(genes[, module_name]) >= median(abs(genes[, module_name])), ]
  write.csv(genes, file = file.path(all_gene_dir, paste0("Selected_Gene_Information_of_", module_name, ".csv")), row.names = F, quote = F)
  write.csv(selected_gene, file = file.path(selected_gene_dir, paste0("Selected_Gene_Information_of_", module_name, ".csv")), row.names = F, quote = F)
}
write.csv(geneInfo, file = file.path(results_dir, "Gene_Information.csv"), row.names = F, quote = F)
qsave(list(geneModuleMembership, geneTraitSignificance, geneInfo, selected_gene), file = file.path(results_dir, "MM&GS_Infor.qs"))
