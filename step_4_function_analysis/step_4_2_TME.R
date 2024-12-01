library(GSVA)
library(stringr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(CIBERSORT)
library(dplyr)
library(data.table)
library(utils)
library(estimate)
library(ggstatsplot)
library(survival)
library(survminer)
library(ComplexHeatmap)
library(circlize)
library(linkET)
library(data.table)
library(qs)

rm(list = ls())
source("scripts/utils/plot_settings.R")

clinical_dir <- "results/step_3_model_evaluation_and_nomogram/processed_clinical_data"
exp_dir <- "data/Datasets/proprocessed_exp"
datasets <- c("TCGA", "E_MTAB_1980", "CPTAC")
protein_coding_gene <- read.delim(file = "data/protein_coding_gene.txt", header = F, row.names = NULL)[, 1]
output_dir <- "results/step_4_function_analysis/TME"
if (!dir.exists(output_dir)) dir.create(output_dir)
# ssGSEA
til_genelist <- read.delim("data/genesets/other_genesets/28 tumor-infiltrating lymphocytes.txt", row.names = 1, header = F)
ssgsea_output_dir <- paste0(output_dir, "/ssGSEA")
ssgsea_score_dir <- paste0(ssgsea_output_dir, "/ssGSEA_scores")
ssgsea_cor_dir <- paste0(ssgsea_output_dir, "/correlation")
ssgsea_graph_dir <- paste0(ssgsea_output_dir, "/graph")
if (!dir.exists(ssgsea_output_dir)) dir.create(ssgsea_output_dir)
if (!dir.exists(ssgsea_score_dir)) dir.create(ssgsea_score_dir)
if (!dir.exists(ssgsea_cor_dir)) dir.create(ssgsea_cor_dir)
if (!dir.exists(ssgsea_graph_dir)) dir.create(ssgsea_graph_dir)
ssgsea_fun <- function(exp, genesets) {
  genesets <- as.matrix(genesets)
  genesets <- apply(genesets, 1, function(row) {
    cleaned_row <- row[!is.na(row) & row != ""]
    as.character(cleaned_row)
  })
  exp <- data.matrix(exp)
  result <- gsva(exp, genesets, mx.diff = FALSE, verbose = T, method = "ssgsea")
  result <- data.frame(result)
  return(result)
}
cor_fun <- function(data, method) {
  data <- as.matrix(data)
  ncol <- ncol(data)
  nrow <- nrow(data)
  x <- data[1, ]
  x <- as.numeric(x)
  pvalue <- c()
  cor <- c()
  for (i in 2:nrow) {
    y <- data[i, ]
    y <- as.numeric(y)
    p <- cor.test(x, y, method = method)$p.value
    pvalue <- c(pvalue, p)
    co <- cor(x, y, method = method)
    cor <- c(cor, co)
  }
  correlation <- cbind(cor, pvalue)
  row.names(correlation) <- row.names(data)[-1]
  correlation <- data.frame(correlation)
  return(correlation)
}
for (dataset in datasets) {
  exp <- fread(file = file.path(exp_dir, paste0(dataset, ".csv"))) %>%
    data.frame(row.names = 1)
  clinical <- qread(file = paste0(clinical_dir, "/", dataset, ".qs"))
  exp <- exp[protein_coding_gene, clinical$CaseID]
  ssgsea_result <- ssgsea_fun(exp = exp, genesets = til_genelist)
  write.csv(ssgsea_result, file = file.path(ssgsea_score_dir, paste0(dataset, "_ssGSEA_scores.csv")))
  cor_analysis_table <- rbind(risk_score = clinical$risk_score, ssgsea_result)
  rs_ssgsea_cor <- cor_fun(data = cor_analysis_table, method = "pearson")
  rs_ssgsea_cor <- arrange(rs_ssgsea_cor, desc(cor))
  write.csv(rs_ssgsea_cor, file = paste0(ssgsea_cor_dir, "/", dataset, "_ssGSEA_cor.csv"))
}

# CIBERSORT
data(LM22)
# update gene symbols
HGNC_data <- fread("data/HGNC.txt")
table(rownames(LM22) %in% HGNC_data$C)
gene_symbols <- rownames(LM22)
not_found_gene <- gene_symbols[!gene_symbols %in% HGNC_data$C]
current_symbol <- c()
for (gene in not_found_gene) {
  judgment <- str_detect(HGNC_data$`Previous symbols`, fixed(gene))
  current_symbol[gene] <- ifelse(T %in% judgment, HGNC_data$C[judgment], "Not")
}
newLM22 <- data.frame(LM22)
newLM22$symbol <- row.names(newLM22)
for (gene in not_found_gene) {
  newLM22$symbol[newLM22$symbol == gene] <- current_symbol[gene]
}
newLM22 <- subset(newLM22, symbol != "Not")
table(newLM22$symbol %in% HGNC_data$C)
rownames(newLM22) <- newLM22$symbol
table(rownames(newLM22) %in% HGNC_data$C)
newLM22 <- newLM22[, -23] %>% as.matrix()

cibersort_output_dir <- paste0(output_dir, "/cibersort")
cibersort_score_dir <- paste0(cibersort_output_dir, "/cibersort_scores")
cibersort_graph_dir <- paste0(cibersort_output_dir, "/graph")
cibersort_all_score_dir <- paste0(cibersort_score_dir, "/cibersort_all_scores")
cibersort_sig_score_dir <- paste0(cibersort_score_dir, "/cibersort_sig_scores")
if (!dir.exists(cibersort_output_dir)) dir.create(cibersort_output_dir)
if (!dir.exists(cibersort_score_dir)) dir.create(cibersort_score_dir)
if (!dir.exists(cibersort_cor_dir)) dir.create(cibersort_cor_dir)
if (!dir.exists(cibersort_graph_dir)) dir.create(cibersort_graph_dir)
if (!dir.exists(cibersort_all_score_dir)) dir.create(cibersort_all_score_dir)
if (!dir.exists(cibersort_sig_score_dir)) dir.create(cibersort_sig_score_dir)

for (dataset in datasets) {
  exp <- fread(file = paste0(exp_dir, "/", dataset, ".csv"), header = T) %>%
    data.frame(row.names = 1)
  exp <- as.matrix(exp)
  exp <- 2^exp - 0.001
  result <- cibersort(sig_matrix = newLM22, mixture_file = exp, perm = 1000, QN = T)
  write.csv(result, file = paste0(
    cibersort_all_score_dir,
    "/", dataset, "_cibersort_scores_all.csv"
  ))
  result <- data.frame(result)
  result <- subset(result, result$P.value < 0.05)[, -c(24:26)]
  write.csv(result, file = paste0(
    cibersort_sig_score_dir, "/", dataset,
    "_cibersort_scores_sig.csv"
  ))
}

# Estimate
estimate_dir <- paste0(output_dir, "/estimate")
estimate_intermediate_data <- paste0(estimate_dir, "/intermediate")
gct_file_dir <- paste0(estimate_dir, "/gct_files")
estimate_score_dir <- paste0(estimate_dir, "/estimate_score")
estimate_cor_dir <- paste0(estimate_dir, "/correlation")
estimate_analysis_table_dir <- paste0(estimate_dir, "/analysis_table")
estimate_graph_dir <- paste0(estimate_dir, "/graph")
if (!dir.exists(estimate_dir)) dir.create(estimate_dir)
if (!dir.exists(estimate_intermediate_data)) dir.create(estimate_intermediate_data)
if (!dir.exists(gct_file_dir)) dir.create(gct_file_dir)
if (!dir.exists(estimate_score_dir)) dir.create(estimate_score_dir)
if (!dir.exists(estimate_cor_dir)) dir.create(estimate_cor_dir)
if (!dir.exists(estimate_analysis_table_dir)) dir.create(estimate_analysis_table_dir)
if (!dir.exists(estimate_graph_dir)) dir.create(estimate_graph_dir)
for (dataset in datasets) {
  exp <- fread(file = paste0(exp_dir, "/", dataset, ".csv"), header = T) %>%
    data.frame(row.names = 1)
  write.table(exp, file = paste0(estimate_intermediate_data, "/", dataset, ".txt"), quote = F, sep = "\t")
  gct_dir <- paste0(gct_file_dir, "/", dataset, ".gct")
  scores_dir <- paste0(estimate_score_dir, "/", dataset, "_estimate_scores.txt")
  filterCommonGenes(
    input.f = paste0(estimate_intermediate_data, "/", dataset, ".txt"),
    output.f = gct_dir
  )
  estimateScore(input.ds = gct_dir, output.ds = scores_dir)
  clinical <- qread(file = paste0(clinical_dir, "/", dataset, ".qs"))
  clinical <- arrange(clinical, risk_score)
  estimate_score <- read.delim(file = scores_dir, skip = 2, row.names = 1)
  estimate_score <- estimate_score[, -1]
  estimate_score <- estimate_score[, clinical$CaseID]
  correlation_table <- rbind(RS = clinical$risk_score, estimate_score)
  correlation <- cor_fun(correlation_table, method = "pearson")
  write.csv(correlation, file = paste0(estimate_cor_dir, "/", dataset, "_cor.csv"))
  write.csv(correlation_table, file = paste0(estimate_analysis_table_dir, "/", dataset, "_analysis_table.csv"))
}

# CIBERSORT visualization
for (dataset in datasets) {
  cibersort_score <- read.csv(file = paste0(cibersort_sig_score_dir, "/", dataset, "_cibersort_scores_sig.csv"), row.names = 1)
  clinical <- qread(file = paste0(clinical_dir, "/", dataset, ".qs"))
  risk_score_group <- data.frame(clinical$risk_score_group)
  rownames(risk_score_group) <- clinical$CaseID
  analysis_table <- merge(cibersort_score, risk_score_group, by = "row.names")
  colnames(analysis_table)[1] <- "CaseID"
  colnames(analysis_table)[ncol(analysis_table)] <- "Group"
  analysis_table <- melt(analysis_table, id = c("CaseID", "Group"))
  analysis_table$Group <- factor(analysis_table$Group, levels = c("Low", "High"), ordered = T)
  colnames(analysis_table)[c(3:4)] <- c("CellType", "Fraction")
  labels <- str_replace_all(levels(analysis_table$CellType), pattern = fixed("."), replacement = " ")
  ggplot(analysis_table, aes(x = CellType, y = Fraction, fill = Group)) +
    geom_boxplot(outlier.shape = 21, color = "gray20", outlier.size = 1.2) +
    theme_bw() +
    labs(x = NULL, title = dataset) +
    my_plot_theme(legend = "right", x.text.angle = 45) +
    scale_y_continuous(limits = c(0, max(analysis_table$Fraction) + 0.1)) +
    scale_fill_manual(values = c(
      "High" = rgb(242, 104, 36, maxColorValue = 255),
      "Low" = rgb(016, 139, 150, maxColorValue = 255)
    )) +
    stat_compare_means(aes(group = Group, label = ..p.signif.., ),
      method = "t.test", size = 2,
      label.y = max(analysis_table$Fraction) + 0.05
    ) +
    scale_x_discrete(labels = labels) +
    theme(legend.key.size = unit(0.3, "cm"))
  ggsave(filename = paste0(cibersort_graph_dir, "/", dataset, "_cibersort.pdf"), width = 6.5, height = 2, device = "pdf")
}

# ESTIMATE visualization
ImmuneScore_graph_dir <- paste0(estimate_graph_dir, "/ImmuneScore")
StromalScore_graph_dir <- paste0(estimate_graph_dir, "/StromalScore")
EstimateScore_graph_dir <- paste0(estimate_graph_dir, "/EstimateScore")
TumorPurity_graph_dir <- paste0(estimate_graph_dir, "/TumorPurity")
if (!dir.exists(ImmuneScore_graph_dir)) dir.create(ImmuneScore_graph_dir)
if (!dir.exists(StromalScore_graph_dir)) dir.create(StromalScore_graph_dir)
if (!dir.exists(EstimateScore_graph_dir)) dir.create(EstimateScore_graph_dir)
if (!dir.exists(TumorPurity_graph_dir)) dir.create(TumorPurity_graph_dir)
for (dataset in datasets) {
  estimate_result <- read.csv(file = paste0(estimate_analysis_table_dir, "/", dataset, "_Analysis_table.csv"), header = T, row.names = 1) %>%
    t() %>%
    data.frame()
  for (score in colnames(estimate_result)[2:5]) {
    cor_result <- cor.test(estimate_result$RS, estimate_result[, score])
    cor_coeff <- cor_result$estimate
    p_value <- cor_result$p.value
    p_value_to_show <- ifelse(p_value < 1e-16,
      "P-value < 1e-16",
      paste0("P-value = ", base::format.pval(p_value, digits = 2))
    )
    scatter_plot <- ggplot(estimate_result, aes(x = RS, y = get(score))) +
      geom_point(size = 1) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "#4DBBD5FF", size = 1) +
      labs(
        x = "Risk Score",
        y = gsub("([a-z])([A-Z])", "\\1 \\2", score),
        title = paste0(dataset),
        subtitle = paste("Pearson's r = ", round(cor_coeff, 2), "\n", p_value_to_show)
      ) +
      my_plot_theme() +
      theme(plot.subtitle = element_text(hjust = 0.5))
    ggsave(filename = paste0(estimate_graph_dir, "/", score, "/", dataset, ".pdf"), width = 1.625, height = 2, device = "pdf")
  }
}

# ssGSEA visualization
for (dataset in datasets) {
  ssgsea_score <- read.csv(file = paste0(ssgsea_score_dir, "/", dataset, "_ssgsea_scores.csv"), row.names = 1) %>% t()
  clinical <- qread(file = paste0(clinical_dir, "/", dataset, ".qs")) %>% arrange(risk_score)
  ssgsea_score <- ssgsea_score[clinical$CaseID, ] %>% data.frame()
  ssgsea_score_cor <- correlate(ssgsea_score)
  rs_ssgsea_cor <- read.csv(file = paste0(ssgsea_cor_dir, "/", dataset, "_ssGSEA_cor.csv"), row.names = NULL)
  rs_ssgsea_cor <- cbind(RS = "Risk Score", rs_ssgsea_cor)
  colnames(rs_ssgsea_cor)[2] <- "Cell Type"
  rs_ssgsea_cor$`Cell Type` <- str_replace_all(rs_ssgsea_cor$`Cell Type`, pattern = " ", replacement = ".")
  rs_ssgsea_cor$logP <- -log10(rs_ssgsea_cor$pvalue)
  rs_ssgsea_cor$sig <- ifelse(rs_ssgsea_cor$pvalue < 0.05, "Sig", "Not_Sig")
  qcorrplot(correlate(ssgsea_score), type = "lower") +
    geom_square() +
    scale_fill_gradientn(
      colors = c("dodgerblue2", "white", "orange", "firebrick2"),
      values = scales::rescale(c(-1, 0, 0.5, 1)),
      limits = c(-1, 1)
    ) +
    geom_couple(
      data = rs_ssgsea_cor,
      aes(color = cor, size = logP, linetype = sig),
      curvature = nice_curvature(),
      label.size = 2, label.family = "Arial"
    ) +
    scale_color_gradient2(
      low = "dodgerblue2",
      high = "firebrick2",
      mid = "white",
      midpoint = 0,
    ) +
    scale_size_continuous(limits = c(0, 20), range = c(0.1, 1.5), breaks = c(0, 5, 10)) +
    scale_linetype_manual(
      values = c("Sig" = "solid", "Not_Sig" = "dashed"),
      labels = c("Sig" = "p<0.05", "Not_Sig" = "p>0.05")
    ) +
    annotate("text", x = 23, y = 27, label = dataset, size = 3, family = "Arial", fontface = "bold") +
    labs(
      x = "",
      y = "",
      fill = "Heatmap\nCor",
      color = "TMRS\nCor",
      size = "P(-log10)", linetype = "Sig"
    ) +
    my_plot_theme(x.text.angle = 90) +
    theme(
      legend.key.height = unit(0.3, "cm"),
      legend.key.width = unit(0.4, "cm"),
      legend.spacing.x = unit(-0.3, "cm"),
      legend.position = "bottom",
      legend.box.margin = margin(-25, 0, -20, -10),
      plot.margin = margin(-10, 0, 0, -10)
    ) +
    guides(
      linetype = guide_legend(direction = "vertical", keywidth = unit(0.6, "cm"), order = 4),
      size = guide_legend(direction = "vertical", order = 3),
      fill = guide_colorbar(barheight = unit(1, "cm"), direction = "vertical", order = 2),
      color = guide_colorbar(barheight = unit(1, "cm"), direction = "vertical", order = 1)
    )
  ggsave(filename = paste0(ssgsea_graph_dir, "/", dataset, "_ssGSEA_cor_graph.pdf"), width = 3.25, height = 4, device = "pdf")
}


# TIP
TIP_output_dir <- file.path(output_dir, "cancer_immunity_cycle")
if (!dir.exists(TIP_output_dir)) dir.create(TIP_output_dir)
TIP <- read.delim(file = "data/ssGSEA.normalized.score.txt", row.names = 1)
TCGA_sample_type <- read.delim(file = "data/Datasets/SampleType/TCGA_SampleType.txt", row.names = NULL)
tumor_sample <- TCGA_sample_type$ID[TCGA_sample_type$Type == "Tumor"]
TIP <- TIP[, tumor_sample]
colnames(TIP) <- substr(colnames(TIP), 1, 12)
clinical <- qread(file = file.path(clinical_dir, "TCGA.qs"))
TIP <- TIP[, clinical$CaseID]
TIP_rs <- rbind(TMRS = clinical$risk_score, TIP)
TIP_rs_cor <- cor_fun(TIP_rs, method = "pearson")
TIP_rs_cor$logP <- -log10(TIP_rs_cor$pvalue)
TIP_rs_cor <- cbind(rs = "TMRS", steps = rownames(TIP_rs_cor), TIP_rs_cor)
rownames(TIP_rs_cor) <- NULL
TIP_rs_cor$steps <- str_replace_all(TIP_rs_cor$steps, pattern = fixed(" "), replacement = ".")
TIP_rs_cor$sig <- ifelse(TIP_rs_cor$pvalue < 0.05, "Sig", "Not Sig")
TIP <- t(TIP) %>% data.frame()
# correlation heatmap
qcorrplot(correlate(TIP), type = "lower") +
  geom_square() +
  scale_fill_gradientn(colors = c("dodgerblue2", "white", "orange", "firebrick2"),
                       values = scales::rescale(c(-1, 0, 0.5, 1)), 
                       limits = c(-1, 1)) +
  geom_couple(data = TIP_rs_cor, 
              aes(color = cor, size = logP, linetype = sig),
              label.size = 2, label.family = "Arial",
              curvature = nice_curvature()) +
  scale_color_gradient2(low = "dodgerblue2", high = "firebrick2", mid = "white", midpoint = 0, limits = c(-0.4, 0.4)) +
  scale_size_continuous(limits = c(0, 20), range = c(0.5, 2), breaks = c(0, 5, 10)) +
  scale_linetype_manual(
    values = c("Sig" = "solid", "Not Sig" = "dashed"),
    labels = c("Sig" = "p<0.05", "Not Sig" = "p>0.05")
  ) +
  labs(
    fill = "Heatmap\nCor",
    color = "TMRS\nCor",
    size = "P(-log10)", linetype = "Sig",
    x = NULL, y = NULL
  ) +
  my_plot_theme(x.text.angle = 90) +
  theme(
    legend.key.height = unit(0.3, "cm"),
    legend.key.width = unit(0.4, "cm"),
    legend.spacing.x = unit(-0.3, "cm"),
    legend.position = "bottom",
    legend.box.margin = margin(-50, 0, 0, -10),
    plot.margin = margin(0, 0, 0, 0)
  ) +
  guides(
    linetype = guide_legend(direction = "vertical", keywidth = unit(0.6, "cm"), order = 4),
    size = guide_legend(direction = "vertical", order = 3),
    fill = guide_colorbar(barheight = unit(1, "cm"), direction = "vertical", order = 2),
    color = guide_colorbar(barheight = unit(1, "cm"), direction = "vertical", order = 1)
  )
ggsave(filename = paste0(TIP_output_dir, "/", "TCGA_TIP_cor_graph.pdf"), width = 5.5, height = 5.5, device = "pdf")
