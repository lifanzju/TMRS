library(umap)
library(Rtsne)
library(ggforce)
library(cowplot)
library(stringr)
library(dplyr)
library(GSVA)
library(data.table)
library(qs)
library(scales)

rm(list = ls())
source("scripts/utils/plot_settings.R")

clinical_dir <- "results/step_3_model_evaluation_and_nomogram/processed_clinical_data"
exp_dir <- "data/Datasets/proprocessed_exp"
datasets <- str_remove_all(list.files(exp_dir), pattern = fixed(".csv"))
protein_coding_gene <- read.delim(file = "data/protein_coding_gene.txt", header = F, row.names = NULL)[, 1]
output_dir <- "results/step_4_function_analysis/pathway_analysis"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = T)
HALLMARK_genesets <- read.delim("data/genesets/other_genesets/h.all.v2023.1.Hs.symbols.txt",
  row.names = 1, header = F
)

# umap analysis
umap_dir <- file.path(output_dir, "umap")
if (!dir.exists(umap_dir)) dir.create(umap_dir)
for (dataset in datasets) {
  exp <- fread(file = file.path(exp_dir, paste0(dataset, ".csv"))) %>%
    data.frame(row.names = 1)
  clinical <- qread(file = paste0(clinical_dir, "/", dataset, ".qs"))
  exp <- exp[protein_coding_gene, clinical$CaseID] %>%
    na.omit() %>%
    t()
  umap_results <- umap(exp, random_state = 527)
  umap_df <- as.data.frame(umap_results$layout)
  umap_df$group <- factor(clinical$risk_score_group)
  scatterplot <- ggplot(umap_df, aes(x = V1, y = V2, color = group, shape = group)) +
    geom_point(alpha = 0.7, size = 1) +
    scale_color_manual(values = c(
      "High" = rgb(242, 104, 36, maxColorValue = 255),
      "Low" = rgb(016, 139, 150, maxColorValue = 255)
    )) +
    labs(x = "UMAP1", y = "UMAP2", title = NULL, color = "TMRS", shape = "TMRS") +
    my_plot_theme() +
    theme(legend.key.size = unit(0.1, "cm"))
  density_x <- ggplot(umap_df, aes(x = V1, fill = group)) +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(
      "High" = rgb(242, 104, 36, maxColorValue = 255),
      "Low" = rgb(016, 139, 150, maxColorValue = 255)
    )) +
    scale_y_continuous() +
    labs(x = NULL, title = dataset, fill = "TMRS") +
    my_plot_theme() +
    theme(
      legend.position = "none", axis.ticks.x = element_blank(), axis.text.x = element_blank(),
      panel.grid = element_blank(), plot.title = element_text(hjust = 0.5)
    )
  density_y <- ggplot(umap_df, aes(x = V2, fill = group)) +
    geom_density(alpha = 0.7) +
    scale_fill_manual(values = c(
      "High" = rgb(242, 104, 36, maxColorValue = 255),
      "Low" = rgb(016, 139, 150, maxColorValue = 255)
    )) +
    labs(x = NULL) +
    coord_flip() +
    my_plot_theme() +
    scale_y_continuous() +
    theme(
      legend.position = "none", axis.ticks.y = element_blank(), axis.text.y = element_blank(),
      panel.grid = element_blank()
    )
  plot_legend <- get_legend(scatterplot)
  scatterplot <- scatterplot + theme(legend.position = "none")
  main_plot <- plot_grid(
    density_x, NULL, scatterplot, density_y,
    ncol = 2, rel_widths = c(4, 1), rel_heights = c(1, 4)
  )

  final_plot <- ggdraw() +
    draw_plot(main_plot, 0, 0, 1, 1) + 
    draw_plot(plot_legend, 0.8, 0.8, 0.2, 0.2)
  ggsave(filename = paste0(umap_dir, "/", dataset, ".pdf"), plot = final_plot, device = "pdf", width = 1.625, height = 1.9)
}

# GSVA
gsva_result_dir <- file.path(output_dir, "gsva_result")
gsva_score_dir <- paste0(gsva_result_dir, "/gsva_scores")
correlation_dir <- paste0(gsva_result_dir, "/cor")
sig_correlation_dir <- paste0(gsva_result_dir, "/sig_cor")

if (!dir.exists(gsva_result_dir)) dir.create(gsva_result_dir)
if (!dir.exists(gsva_score_dir)) dir.create(gsva_score_dir)
if (!dir.exists(correlation_dir)) dir.create(correlation_dir)
if (!dir.exists(sig_correlation_dir)) dir.create(sig_correlation_dir)

gsva_fun <- function(exp, genesets) {
  genesets <- as.matrix(genesets)
  genesets <- apply(genesets, 1, function(row) {
    cleaned_row <- row[!is.na(row) & row != ""]
    as.character(cleaned_row)
  })
  exp <- data.matrix(exp)
  result <- gsva(exp, genesets, mx.diff = FALSE, verbose = T)
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

for (dataset in c("TCGA","E_MTAB_1980","CPTAC")) {
  exp <- fread(file = file.path(exp_dir, paste0(dataset, ".csv"))) %>%
    data.frame(row.names = 1)
  clinical <- qread(file = paste0(clinical_dir, "/", dataset, ".qs"))
  exp <- exp[protein_coding_gene, clinical$CaseID]
  gsva_result <- gsva_fun(exp = exp, genesets = HALLMARK_genesets)
  cor_analysis_table <- rbind(risk_score = clinical$risk_score, gsva_result)
  rs_gsva_cor <- cor_fun(data = cor_analysis_table, method = "pearson")
  rs_gsva_cor <- arrange(rs_gsva_cor, desc(cor))
  write.csv(gsva_result, file = paste0(gsva_score_dir, "/", dataset, ".csv"))
  write.csv(rs_gsva_cor, file = paste0(correlation_dir, "/", dataset, ".csv"))
}

# take intersection of significant pathways
cor_pos <- list()
cor_neg <- list()
for (dataset in c("TCGA","E_MTAB_1980","CPTAC")) {
  all_cor <- read.csv(file = paste0(correlation_dir, "/", dataset, ".csv"), row.names = 1)
  postive_correlated <- subset(all_cor, cor > 0 & pvalue < 0.05)
  negative_correlated <- subset(all_cor, cor < 0 & pvalue < 0.05)
  cor_pos[[dataset]] <- rownames(postive_correlated)
  cor_neg[[dataset]] <- rownames(negative_correlated)
}
intersects <- function(x) {
  Reduce(intersect, x)
}
intersection_pos <- intersects(cor_pos)
intersection_neg <- intersects(cor_neg)
terms_selected <- c(intersection_pos, intersection_neg)

# draw bubble plot
proper_name <- function(data) {
  data <- str_remove_all(data, pattern = "HALLMARK_")
  data <- str_replace_all(data, pattern = "_", replacement = " ")
  data <- str_to_title(data)
}
selected_pathways_all <- data.frame()
for (dataset in datasets) {
  cor <- read.csv(file = paste0(correlation_dir, "/", dataset, ".csv"), row.names = 1)
  cor <- cor[terms_selected, ]
  cor$term <- terms_selected
  cor$dataset <- dataset
  selected_pathways_all <- rbind(selected_pathways_all, cor)
}
selected_pathways_all$dataset <- str_replace_all(selected_pathways_all$dataset,
  pattern = "_", replacement = "-"
)
selected_pathways_all$logP <- -log10(selected_pathways_all$pvalue)
selected_pathways_all$term <- proper_name(selected_pathways_all$term)
selected_pathways_all$term <- str_replace_all(selected_pathways_all$term, "Il6 Jak Stat3", replacement = "IL6-JAK-STAT3")
selected_pathways_all$term <- str_replace_all(selected_pathways_all$term, "Epithelial Mesenchymal Transition", replacement = "EMT")
selected_pathways_all$term <- str_replace_all(selected_pathways_all$term, "G2m", replacement = "G2M")
selected_pathways_all$term <- str_replace_all(selected_pathways_all$term, "E2f", replacement = "E2F")
selected_pathways_all$term <- factor(selected_pathways_all$term, levels = unique(selected_pathways_all$term), ordered = T)

ggplot(selected_pathways_all, aes(x = term, y = dataset, fill = cor, size = logP)) +
  geom_point(shape = 21, color = "black") +
  scale_fill_gradientn(
    name = "Correlation\n(Pearson's)",
    colours = c(
      "#1663A9", rgb(033, 158, 188, maxColorValue = 255),
      "white", rgb(255, 183, 005, maxColorValue = 255), "#B9181A"
    ),
    values = rescale(c(-0.5, -0.25, 0, 0.25, 0.5)),
    breaks = c(-0.5, 0, 0.5)
  ) +
  scale_size_continuous(
    name = "-log10\np-value",
    limits = c(-0.001, max(selected_pathways_all$logP)),
    breaks = c(
      round(min(selected_pathways_all$logP)),
      mean(c(
        round(min(selected_pathways_all$logP)),
        round(max(selected_pathways_all$logP))
      )),
      floor(max(selected_pathways_all$logP))
    ),
    range = c(0.5, 5)
  ) +
  geom_vline(xintercept = c(length(intersection_pos) + 0.5)) +
  labs(
    x = NULL,
    y = NULL,
    title = "HALLMARK Signatures"
  ) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 20))+
  my_plot_theme(x.text.angle = 60,legend = "top") +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_line(colour = "gray"),
    legend.key.size = unit(0.2, "cm"),
    legend.spacing = unit(0, "cm"),
    legend.box.margin = margin(0, 0, -10, -10),
  )
ggsave(filename = file.path(gsva_result_dir,"bubble_plot.pdf"), device = "pdf", height = 2, width = 3.25)
