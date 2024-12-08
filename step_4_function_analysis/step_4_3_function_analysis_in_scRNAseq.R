library(Seurat)
library(Scissor)
library(AUCell)
library(clusterProfiler)
library(org.Hs.eg.db)
library(linkET)
library(monocle)
library(harmony)
library(ggsci)
library(qs)
library(gbm)
library(stringr)

rm(list = ls())
source("scripts/utils/plot_settings.R")

output_dir <- "results/step_4_function_analysis/scRNAseq"
graph_output_dir <- paste0(output_dir, "/graph")
data_output_dir <- paste0(output_dir, "/data")
if (!dir.exists(graph_output_dir)) {
  dir.create(graph_output_dir, recursive = TRUE)
}
if (!dir.exists(data_output_dir)) {
  dir.create(data_output_dir, recursive = TRUE)
}

# TMRS in proximal tubule cells
sce_mali_pro <- qread("results/step_1_hdwgcna_wgcna/pathway_enrichment_score/mali_pro_cells/mali_pro_cells_auc_calculated.qs")
sce_mali_pro$CellType <- factor(sce_mali_pro$CellType, levels = c("Proximal Tubule Cells", "Transitional Cells", "Malignant Cells"), ordered = T)
p1 <- DimPlot(sce_mali_pro) +
  ggtitle("Cell Type") +
  labs(x = "UMAP1", y = "UMAP2") +
  scale_color_manual(values = c(
    rgb(114, 188, 213, maxColorValue = 255),
    rgb(231, 098, 084, maxColorValue = 255),
    rgb(255, 208, 111, maxColorValue = 255)
  ), labels = str_wrap(levels(risk_score_data$CellType), width = 15)) +
  my_sc_plot_theme() +
  theme(
    legend.key.size = unit(0.2, "cm"),
    legend.position = c(0.75, 0.1)
  )
ggsave(plot = p1, filename = paste0(graph_output_dir, "/cell_type.pdf"), width = 1.625, height = 2)
final_model <- qread(file = "results/step_2_model_construction/machine_learning/final_model.qs")
genes_in_model <- final_model$var.names

sce_exp <- as.matrix(sce_mali_pro@assays$RNA@data[genes_in_model, ]) %>% t()
sce_scaled_exp <- scale(sce_exp) %>% data.frame()
sce_risk_score <- as.numeric(predict(final_model, newdata = sce_scaled_exp, type = "link"))
names(sce_risk_score) <- rownames(sce_scaled_exp)
sce_mali_pro$risk_score <- sce_risk_score
p2 <- FeaturePlot(sce_mali_pro, features = "risk_score", min.cutoff = -1, max.cutoff = "q95") +
  viridis::scale_color_viridis(option = "A") +
  ggtitle("TMRS") +
  my_sc_plot_theme() +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(
    legend.key.size = unit(0.2, "cm"),
    legend.position = c(0.75, 0.1)
  )
ggsave(plot = p2, filename = paste0(graph_output_dir, "/risk_score.pdf"), width = 1.625, height = 2)
risk_score_data <- sce_mali_pro@meta.data
risk_score_data$CellType <- factor(risk_score_data$CellType, levels = c("Proximal Tubule Cells", "Transitional Cells", "Malignant Cells"), ordered = T)
comparasions <- list(
  c("Proximal Tubule Cells", "Transitional Cells"),
  c("Transitional Cells", "Malignant Cells"),
  c("Proximal Tubule Cells", "Malignant Cells")
)
p3 <- ggplot(risk_score_data, aes(x = CellType, y = risk_score, fill = CellType)) +
  geom_violin() +
  scale_fill_manual(values = c(
    rgb(114, 188, 213, maxColorValue = 255),
    rgb(255, 208, 111, maxColorValue = 255),
    rgb(231, 098, 084, maxColorValue = 255)
  )) +
  geom_boxplot(width = 0.1) +
  labs(x = NULL, y = "Risk Score", title = "") +
  scale_x_discrete(labels = str_wrap(levels(risk_score_data$CellType), width = 10)) +
  stat_compare_means(
    method = "anova",
    label.y = max(risk_score_data$risk_score) + 1.3,
    size = 2
  ) +
  stat_compare_means(
    comparisons = comparasions,
    method = "t.test",
    hide.ns = F,
    label = "p.signif"
  ) +
  my_plot_theme(legend = "none")
ggsave(plot = p3, filename = paste0(graph_output_dir, "/risk_score_compare.pdf"), width = 1.625, height = 2)


# TMRS in malignant cells
sce_mali <- subset(sce_mali_pro, idents = c("Malignant Cells"))
sce_mali <- NormalizeData(sce_mali, normalization.method = "LogNormalize", scale.factor = 1e4)
sce_mali <- FindVariableFeatures(sce_mali)
sce_mali <- ScaleData(sce_mali)
sce_mali <- RunPCA(sce_mali, features = VariableFeatures(object = sce_mali))
sce_mali <- RunHarmony(sce_mali, group.by.vars = "orig.ident", plot_convergence = T)
ElbowPlot(sce_mali, ndims = 50)
sce_mali <- FindNeighbors(sce_mali, dims = 1:20, reduction = "harmony")
sce_mali <- FindClusters(sce_mali, resolution = 1)
sce_mali <- RunUMAP(sce_mali, dims = 1:20, reduction = "harmony")

sce_mali_exp <- as.matrix(sce_mali@assays$RNA@data[genes_in_model, ]) %>% t()
sce_mali_scaled_exp <- scale(sce_mali_exp) %>% data.frame()
sce_mali_risk_score <- as.numeric(predict(final_model, newdata = sce_mali_scaled_exp, type = "link"))
names(sce_mali_risk_score) <- rownames(sce_mali_scaled_exp)
sce_mali$risk_score <- sce_mali_risk_score
sce_mali$rs_group <- ifelse(sce_mali$risk_score >= median(sce_mali$risk_score), "High", "Low")
Idents(sce_mali) <- "rs_group"
high_rs_cells <- rownames(sce_mali@meta.data)[sce_mali$rs_group == "High"]
low_rs_cells <- rownames(sce_mali@meta.data)[sce_mali$rs_group == "Low"]
sce_mali_cluster <- DimPlot(sce_mali,
  group.by = "seurat_clusters",
  cols = c(pal_npg(alpha = 0.8)(7), pal_jama(alpha = 0.8)(7))
) +
  labs(
    title = "Malignant Cell Cluster",
    x = "UMAP1", y = "UMAP2"
  ) +
  my_sc_plot_theme() +
  guides(color = guide_legend(ncol = 2)) +
  theme(
    legend.key.height = unit(0.2, "cm"),
    legend.key.width = unit(0.2, "cm"),
    legend.position = c(0.75, 0.1)
  )
ggsave(filename = paste0(graph_output_dir, "/Malignant Cells Clusters.pdf"), width = 1.625, height = 2)
sce_mali_rs_group <- DimPlot(sce_mali, cols = pal_npg(alpha = 0.8)(2)) +
  labs(
    title = "TMRS Group in\nMalignant Cells",
    x = "UMAP1",
    y = "UMAP2"
  ) +
  my_sc_plot_theme() +
  theme(legend.position = c(0.75, 0.1))
ggsave(filename = paste0(graph_output_dir, "/Malignant Cells RS Group.pdf"), width = 1.625, height = 2)
qsave(sce_mali, file = paste0(data_output_dir, "/sce_mali.qs"))
sce_mali_rs <- FeaturePlot(sce_mali, features = "risk_score", min.cutoff = -0.8, max.cutoff = 0.8) +
  viridis::scale_color_viridis(option = "A") +
  labs(
    title = "TMRS in\nMalignant Cells",
    x = "UMAP1",
    y = "UMAP2"
  ) +
  my_sc_plot_theme() +
  theme(
    legend.position = c(0.75, 0.1),
    legend.key.size = unit(0.2, "cm")
  )
ggsave(
  plot = sce_mali_rs, filename = paste0(graph_output_dir, "/Malignant Cells RS.pdf"),
  width = 1.625, height = 2
)


# scissor
bulk_exp <- fread(file = "data/Datasets/proprocessed_exp/TCGA.csv") %>%
  data.frame(row.names = 1) %>%
  as.matrix()
bulk_exp <- 2^bulk_exp - 0.001 # return to TPM
bulk_survival <- qread(file = "results/step_3_model_evaluation_and_nomogram/processed_clinical_data/TCGA.qs")
bulk_exp <- bulk_exp[, bulk_survival$CaseID]
phenotype <- dplyr::select(bulk_survival, c("OS", "OS.Censor"))
colnames(phenotype) <- c("time", "status")
rownames(phenotype) <- bulk_survival$CaseID
scissor_results <- Scissor(bulk_exp,
  sce_mali,
  phenotype,
  family = "cox",
  Save_file = file.path(data_output_dir, "scissor_input.RData")
)
scissor_select <- rep("Background cells", ncol(sce_mali))
names(scissor_select) <- colnames(sce_mali)
scissor_select[scissor_results$Scissor_pos] <- "scissor+ cells"
scissor_select[scissor_results$Scissor_neg] <- "scissor- cells"
sce_mali <- AddMetaData(sce_mali, metadata = scissor_select, col.name = "scissor")
qsave(sce_mali, file = file.path(data_output_dir, "sce_mali_scissor.qs"))
scissor_plot <- DimPlot(sce_mali,
  group.by = "scissor",
  cols = c("grey", "royalblue", "indianred1")
) +
  ggtitle("Scissor Analysis") +
  theme(plot.title = element_text(hjust = 0.5)) +
  my_sc_plot_theme() +
  labs(x = "UMAP1", y = "UMAP2") +
  theme(
    legend.position = c(0.7, 0.1),
    legend.key.size = unit(0.2, "cm")
  )
ggsave(scissor_plot, filename = paste0(graph_output_dir, "/Scissor Analysis.pdf"), width = 1.625, height = 2)
cell_count <- table(sce_mali$rs_group, sce_mali$scissor)
proportion_scissor_group <- prop.table(cell_count, margin = 2) %>% data.frame()
proportion_scissor_group$Var2 <- factor(proportion_scissor_group$Var2,
  levels = c("scissor- cells", "Background cells", "scissor+ cells"),
  labels = c("Scissor-\nCells", "Background\nCells", "Scissor+\nCells")
)
stack_plot <- ggplot(proportion_scissor_group) +
  geom_bar(aes(x = Var2, y = Freq, fill = Var1), stat = "identity", width = 0.7, size = 0.5, colour = "#222222") +
  scale_fill_manual(values = pal_npg(alpha = 0.8)(2)) +
  theme_classic() +
  labs(x = "", y = "Proportion", fill = "TMRS") +
  my_plot_theme(legend = "top") +
  theme(
    legend.key.size = unit(0.2, "cm"),
    legend.box.margin = margin(0, 0, -15, 0)
  )
ggsave(stack_plot, filename = paste0(graph_output_dir, "/Scissor Proportion.pdf"), width = 1.625, height = 2)

risk_score_data <- sce_mali@meta.data
risk_score_data$scissor <- factor(risk_score_data$scissor,
  levels = c("scissor- cells", "Background cells", "scissor+ cells"),
  labels = c("Scissor-\nCells", "Background\nCells", "Scissor+\nCells")
)
comparasions <- list(c("Scissor-\nCells", "Background\nCells"), c("Background\nCells", "Scissor+\nCells"), c("Scissor+\nCells", "Scissor-\nCells"))
scissor_rs_vlnplot <- ggplot(risk_score_data, aes(x = scissor, y = risk_score, fill = scissor)) +
  geom_violin() +
  scale_fill_manual(values = c(rgb(114, 188, 213, maxColorValue = 255), rgb(255, 208, 111, maxColorValue = 255), rgb(231, 098, 084, maxColorValue = 255))) +
  geom_boxplot(width = 0.1) +
  labs(x = NULL, y = "Risk Score") +
  stat_compare_means(
    method = "anova",
    label.x.npc = 0.38,
    label.y = max(risk_score_data$risk_score) + 1.2,
    size = 2
  ) +
  stat_compare_means(comparisons = comparasions, method = "t.test", hide.ns = F, label = "p.signif") +
  my_plot_theme(legend = "none")
ggsave(plot = scissor_rs_vlnplot, filename = paste0(graph_output_dir, "/Scissor RS Violinplot.pdf"), width = 1.625, height = 2)


# function analysis
exp_matrix <- sce_mali@assays$RNA@counts
cells_rankings <- AUCell_buildRankings(exp_matrix)
HALLMARK_genesets <- read.delim(file = "data/genesets/other_genesets/h.all.v2023.1.Hs.symbols.txt", row.names = 1, header = F)
KEGG_genesets <- read.delim(file = "data/genesets/other_genesets/c2.cp.kegg.v2023.1.Hs.symbols.txt", row.names = 1, header = F)

genesets <- list(HALLMARK_genesets, KEGG_genesets)
names(genesets) <- c("HALLMARK", "KEGG")
for (i in 1:length(genesets)) {
  geneset <- genesets[[i]]
  geneset[geneset == ""] <- NA
  geneset_list <- c()
  nrow <- nrow(geneset)
  for (n in 1:nrow) {
    a <- geneset[n, ]
    a <- a[!is.na(a)]
    a <- list(as.character(as.matrix(a))) 
    geneset_list <- c(geneset_list, a)
  }
  names(geneset_list) <- rownames(geneset)
  cells_AUC <- AUCell_calcAUC(geneSets = geneset_list, rankings = cells_rankings)
  AUC_results <- cells_AUC@assays@data$AUC
  risk_score <- sce_mali@meta.data$risk_score
  names(risk_score) <- rownames(sce_mali@meta.data)
  risk_score <- risk_score[colnames(AUC_results)]
  cor_analysis_table <- rbind(RS = risk_score, AUC_results)
  cor_analysis_table <- t(cor_analysis_table) %>%
    data.frame() %>%
    dplyr::arrange(RS) %>%
    t()
  AUC_RS_Cor <- cor_fun(cor_analysis_table, method = "spearman")
  AUC_RS_Cor <- cbind(RS = "Risk Score", Terms = rownames(AUC_RS_Cor), AUC_RS_Cor)
  AUC_RS_Cor$logP <- -log10(AUC_RS_Cor$pvalue)
  AUC_RS_Cor$Sig <- ifelse(AUC_RS_Cor$pvalue < 0.05, "Sig", "Not_Sig")
  write.csv(AUC_results, file = paste0(data_output_dir, "/", names(genesets)[i], "_AUC.csv"))
  write.csv(AUC_RS_Cor, file = paste0(data_output_dir, "/", names(genesets)[i], "_AUC_RS_Cor.csv"))
}
Proper_name <- function(data, geneset) {
  data <- str_remove_all(data, pattern = paste0(geneset, "_"))
  data <- str_replace_all(data, pattern = "_", replacement = " ")
  data <- str_to_title(data)
  data <- str_replace_all(data, pattern = "Nod like", replacement = "NOD-like")
  data <- str_replace_all(data, pattern = "Mapk", replacement = "MAPK")
  data <- str_replace_all(data, pattern = "To", replacement = "to")
  data <- str_replace_all(data, pattern = "Atp", replacement = "ATP")
  data <- str_replace_all(data, pattern = "Tca", replacement = "TCA")
  data <- str_replace_all(data, pattern = "Of", replacement = "of")
  data <- str_replace_all(data, pattern = "On", replacement = "on")
  data <- str_replace_all(data, pattern = "Or", replacement = "or")
  data <- str_replace_all(data, pattern = "Nad P H", replacement = "NADPH")
  data <- str_replace_all(data, pattern = "Nad", replacement = "NAD")
  data <- str_replace_all(data, pattern = "NADh", replacement = "NADH")
  data <- str_replace_all(data, pattern = "As", replacement = "as")
  data <- str_replace_all(data, pattern = "Il", replacement = "IL")
  data <- str_replace_all(data, pattern = "Nfkb", replacement = "NFKB")
  data <- str_replace_all(data, pattern = "E2f", replacement = "E2F")
  data <- str_replace_all(data, pattern = "G2m", replacement = "G2M")
  data <- str_replace_all(data, pattern = "Tgf", replacement = "TGF")
  data <- str_replace_all(data, pattern = "Kras", replacement = "KRAS")
  data <- str_replace_all(data, pattern = "Uv", replacement = "UV")
  data <- str_replace_all(data, pattern = "Jak", replacement = "JAK")
  data <- str_replace_all(data, pattern = "Stat", replacement = "STAT")
  data <- str_replace_all(data, pattern = "Mtorc", replacement = "mTORC")
  data <- str_replace_all(data, pattern = "Pi3k", replacement = "PI3K")
  data <- str_replace_all(data, pattern = "Mtor", replacement = "mTOR")
  data <- str_replace_all(data, pattern = "Dna", replacement = "DNA")
  data <- str_replace_all(data, pattern = "Via", replacement = "via")
  data <- str_replace_all(data, pattern = "Tnfa", replacement = "TNFA")
  data <- str_replace_all(data, pattern = "Rrna", replacement = "rRNA")
  data <- str_replace_all(data, pattern = "Mrna", replacement = "mRNA")
  data <- str_replace_all(data, pattern = "Utr", replacement = "UTR")
  data <- str_replace_all(data, pattern = "Ecm", replacement = "ECM")
  data <- str_replace_all(data, pattern = "assembly", replacement = "Assembly")
  data <- str_replace_all(data, pattern = "And", replacement = "and")
  data <- str_replace_all(data, pattern = "organic", replacement = "Organic")
}
# bubble plot
for (geneset in names(genesets)) {
  AUC_RS_Cor <- read.csv(file = paste0(data_output_dir, "/", geneset, "_AUC_RS_Cor.csv"), row.names = 1)
  AUC_RS_Cor <- arrange(AUC_RS_Cor, desc(cor)) %>% 
    subset(Sig == "Sig")
  AUC_RS_Cor$Terms <- Proper_name(AUC_RS_Cor$Terms, geneset) %>% 
    str_wrap(width = 20)
  pos_15 <- AUC_RS_Cor[rev(1:15), ]
  neg_15 <- AUC_RS_Cor[(nrow(AUC_RS_Cor) - 14):nrow(AUC_RS_Cor), ]
  pos_15$Terms <- factor(pos_15$Terms, levels = pos_15$Terms)
  neg_15$Terms <- factor(neg_15$Terms, levels = neg_15$Terms)
  
  up <- ggplot(pos_15,aes(x=logP, y= Terms, size= logP, col=cor,fill=cor)) + geom_point(shape=21)+ 
    scale_color_gradient(low = rgb(255,208,111,maxColorValue = 255,alpha=200), high = "firebrick3")+
    scale_fill_gradient(low = rgb(255,208,111,maxColorValue = 255,alpha = 200), high = "firebrick3")+
    scale_size_continuous(range = c(3,6))+
    labs(y=NULL,x="-log10 P-value",title =paste0(geneset," - Up"),fill="Cor",color="Cor")+
    my_plot_theme()+
    theme(legend.position = c(0.8,0.4),
          legend.key.height = unit(0.4, "cm"),
          legend.spacing.y = unit(-0.3, "cm"))
  ggsave(plot = up, filename = paste0(graph_output_dir, "/", geneset, "_bubble_plot_up.pdf"), width = 3.25, height = 3.5)
  
  down <- ggplot(neg_15,aes(x=logP, y= Terms, size= logP, col=cor,fill=cor)) + geom_point(shape=21)+ 
    scale_color_gradient(low = "dodgerblue4", high = rgb(114,188,213,maxColorValue = 255,alpha=200))+
    scale_fill_gradient(low = "dodgerblue4", high = rgb(114,188,213,maxColorValue = 255,alpha=200))+
    scale_size_continuous(range = c(3,6))+
    labs(y=NULL,x="-log10 P-value",title =paste0(geneset," - Down"),fill="Cor",color="Cor")+
    my_plot_theme()+
    theme(legend.position = c(0.8,0.4),
          legend.key.height = unit(0.4, "cm"),
          legend.spacing.y = unit(-0.3, "cm"))
  ggsave(plot = down, filename = paste0(graph_output_dir, "/", geneset, "_bubble_plot_down.pdf"), width = 3.25, height = 3.5)
}






