library(maftools)
library(tidyverse)
library(stringr)
library(forestplot)
library(ggpubr)
library(mclust)
library(qs)

rm(list = ls())
maf_dir <- "data/genomic/SNV"
clinical_dir <- "results/step_3_model_evaluation_and_nomogram/processed_clinical_data"
output_dir <- "results/step_5_genomic"
datasets <- c("TCGA", "CPTAC")
source("scripts/utils/plot_settings.R")
mutation_data_output <- paste0(output_dir, "/", "mutation_data")
mutation_graph_output <- paste0(output_dir, "/", "mutation_graph")
tmb_graph_output <- paste0(output_dir, "/", "tmb_graph")
math_graph_output <- paste0(output_dir, "/", "math_graph")
somatic_interaction_graph_output <- paste0(output_dir, "/", "somatic_interacrion_graph")
if (!dir.exists(mutation_data_output)) dir.create(mutation_data_output, recursive = T)
if (!dir.exists(mutation_graph_output)) dir.create(mutation_graph_output)
if (!dir.exists(tmb_graph_output)) dir.create(tmb_graph_output)
if (!dir.exists(math_graph_output)) dir.create(math_graph_output)
if (!dir.exists(somatic_interaction_graph_output)) dir.create(somatic_interaction_graph_output)
high_risk_color <- c(
  "Missense_Mutation" = rgb(231, 098, 084, maxColorValue = 255),
  "Frame_Shift_Del" = rgb(247, 170, 088, maxColorValue = 255),
  "Nonsense_Mutation" = rgb(255, 230, 183, maxColorValue = 255),
  "Frame_Shift_Ins" = rgb(63, 171, 71, maxColorValue = 255),
  "Splice_Site" = rgb(093, 163, 157, maxColorValue = 255),
  "In_Frame_Del" = rgb(170, 220, 224, maxColorValue = 255),
  "In_Frame_Ins" = rgb(082, 143, 173, maxColorValue = 255),
  "Translation_Start_Site" = "#0071C2",
  "Nonstop_Mutation" = rgb(13, 91, 38, maxColorValue = 255),
  "Multi_Hit" = rgb(158, 049, 080, maxColorValue = 255)
)
low_risk_color <- c(
  "Missense_Mutation" = "#0071C2",
  "Frame_Shift_Del" = rgb(082, 143, 173, maxColorValue = 255),
  "Nonsense_Mutation" = rgb(170, 220, 224, maxColorValue = 255),
  "Frame_Shift_Ins" = rgb(63, 171, 71, maxColorValue = 255),
  "Splice_Site" = rgb(093, 163, 157, maxColorValue = 255),
  "In_Frame_Del" = rgb(13, 91, 38, maxColorValue = 255),
  "In_Frame_Ins" = rgb(255, 230, 183, maxColorValue = 255),
  "Translation_Start_Site" = rgb(247, 170, 088, maxColorValue = 255),
  "Nonstop_Mutation" = rgb(231, 098, 084, maxColorValue = 255),
  "Multi_Hit" = rgb(158, 049, 080, maxColorValue = 255)
)
titv_col <- c(
  "C>T" = rgb(231, 098, 084, maxColorValue = 255),
  "T>G" = rgb(255, 230, 183, maxColorValue = 255),
  "T>C" = rgb(247, 170, 088, maxColorValue = 255),
  "T>A" = rgb(63, 171, 71, maxColorValue = 255),
  "C>G" = rgb(170, 220, 224, maxColorValue = 255),
  "C>A" = "#0071C2"
)
for (dataset in datasets) {
  mafFilePath <- dir(path = paste0(maf_dir, "/", dataset), pattern = "masked.maf.gz$", full.names = T, recursive = T)
  mafdata <- list()
  if (dataset == "CPTAC") {
    sample_list <- read.delim("data/gdc_sample_sheet.2023-09-17.tsv")
    filenames <- dir(path = paste0(maf_dir, "/", dataset), pattern = "masked.maf.gz$", full.names = F, recursive = T)
    filenames <- sapply(str_split(filenames, pattern = fixed("/")), function(x) {
      x[2]
    })
    CaseID <- sample_list$Case.ID[match(filenames, sample_list$File.Name)]
    CaseID <- sapply(str_split(CaseID, pattern = fixed(",")), function(x) {
      x[1]
    })
  }
  for (i in 1:length(mafFilePath)) {
    check_file <- read.csv(mafFilePath[i], sep = "\t", skip = 7, stringsAsFactors = F)
    if (nrow(check_file) > 0) {
      fila <- read.maf(mafFilePath[i], isTCGA = ifelse(str_detect(dataset, pattern = "TCGA"), T, F))
      if (dataset == "CPTAC") {
        fila@data$Tumor_Sample_Barcode <- CaseID[i]
        fila@clinical.data$Tumor_Sample_Barcode <- CaseID[i]
      }
      mafdata <- c(mafdata, fila)
    } else {
      print(mafFilePath[i])
    }
  }
  mutation_data <- merge_mafs(mafdata)
  clinical <- qread(file = paste0(clinical_dir, "/", dataset, ".qs"))
  clinical$CaseID <- str_replace_all(clinical$CaseID, pattern = fixed("."), replacement = "-")
  mutation_data@clinical.data <- merge(mutation_data@clinical.data, clinical, by.x = "Tumor_Sample_Barcode", by.y = "CaseID", all.x = TRUE)
  qsave(mutation_data, file = paste0(mutation_data_output, "/", dataset, "_all.qs"))
  high_risk_sample <- subset(clinical, risk_score_group == "High")
  low_risk_sample <- subset(clinical, risk_score_group == "Low")
  mutation_data_high_risk <- subsetMaf(mutation_data, tsb = high_risk_sample$CaseID)
  mutation_data_low_risk <- subsetMaf(mutation_data, tsb = low_risk_sample$CaseID)
  qsave(mutation_data_high_risk, file = paste0(mutation_data_output, "/", dataset, "_high_risk.qs"))
  qsave(mutation_data_low_risk, file = paste0(mutation_data_output, "/", dataset, "_low_risk.qs"))
  pdf(
    file = paste0(mutation_graph_output, "/", dataset, "_high_risk_summary_plot.pdf"),
    width = 6.5, height = 4, fonts = "Arial", pointsize = 10
  )
  plotmafSummary(mutation_data_high_risk, rmOutlier = TRUE, addStat = "mean", dashboard = TRUE)
  dev.off()
  pdf(
    file = paste0(mutation_graph_output, "/", dataset, "_high_risk_oncoplot.pdf"),
    width = 3.25, height = 2.2, pointsize = 7, fonts = "Arial"
  )
  oncoplot(
    top = 15, mutation_data_high_risk, draw_titv = T,
    colors = high_risk_color, titv_col = titv_col,
    bgCol = "white", borderCol = "gray80",
  )
  dev.off()
  pdf(
    file = paste0(mutation_graph_output, "/", dataset, "_low_risk_summary_plot.pdf"),
    width = 6.5, height = 4, fonts = "Arial", pointsize = 10
  )
  plotmafSummary(mutation_data_low_risk, rmOutlier = TRUE, addStat = "mean", dashboard = TRUE)
  dev.off()
  pdf(
    file = paste0(mutation_graph_output, "/", dataset, "_low_risk_oncoplot.pdf"),
    width = 3.25, height = 2.2, pointsize = 7, fonts = "Arial"
  )
  oncoplot(top = 15, mutation_data_low_risk, draw_titv = T, 
           colors = low_risk_color, 
           titv_col = titv_col, 
           bgCol = "white", borderCol = "gray80")
  dev.off()
  # Fisher test
  high_risk_genes <- getGeneSummary(mutation_data_high_risk)
  low_risk_genes <- getGeneSummary(mutation_data_low_risk)
  high_top_15 <- head(high_risk_genes, 15)$Hugo_Symbol
  low_top_15 <- head(low_risk_genes, 15)$Hugo_Symbol
  genes <- c("VHL", "PBRM1", "TTN", "SETD2", "BAP1", "MUC16", "MTOR", "KDM5C", "TP53")
  mutation_infor_high_rc <- subset(high_risk_genes, Hugo_Symbol %in% genes)
  high_rc_table <- data.frame(
    Mutated = mutation_infor_high_rc$total,
    Normal = (as.numeric(mutation_data_high_risk@summary[3, 2])
    - mutation_infor_high_rc$total)
  )
  rownames(high_rc_table) <- mutation_infor_high_rc$Hugo_Symbol
  mutation_infor_low_rc <- subset(low_risk_genes, Hugo_Symbol %in% genes)
  low_rc_table <- data.frame(
    Mutated = mutation_infor_low_rc$total,
    Normal = (as.numeric(mutation_data_low_risk@summary[3, 2]) -
      mutation_infor_low_rc$total)
  )
  row.names(low_rc_table) <- mutation_infor_low_rc$Hugo_Symbol
  odds_ratio <- c()
  lower_conf_int <- c()
  upper_conf_int <- c()
  pvalue <- c()
  High_RC_Num <- c()
  Low_RC_Num <- c()
  for (gene in genes) {
    Analysis_table <- rbind(high_rc_table[gene, ], low_rc_table[gene, ])
    rownames(Analysis_table) <- c("High_RC", "Low_RC")
    Fisher_result <- fisher.test(Analysis_table)
    odds_ratio <- c(odds_ratio, Fisher_result$estimate)
    lower_conf_int <- c(lower_conf_int, Fisher_result$conf.int[1])
    upper_conf_int <- c(upper_conf_int, Fisher_result$conf.int[2])
    pvalue <- c(pvalue, Fisher_result$p.value)
    High_RC_Num <- c(High_RC_Num, paste0(Analysis_table$Mutated[1], "/", Analysis_table$Normal[1]))
    Low_RC_Num <- c(Low_RC_Num, paste0(Analysis_table$Mutated[2], "/", Analysis_table$Normal[2]))
  }
  result_table <- data.frame(
    odds_ratio = odds_ratio, lower_conf_int = lower_conf_int,
    upper_conf_int = upper_conf_int, High_RC_Num = High_RC_Num,
    Low_RC_Num = Low_RC_Num, pvalue = pvalue
  )
  rownames(result_table) <- genes

  OR.text <- paste(round(result_table$odds_ratio, 3),
    " (", round(result_table$lower_conf_int, 3),
    "-", round(result_table$upper_conf_int, 3), ")",
    sep = ""
  )
  table.text <- cbind(
    c(NA, "Gene", rownames(result_table)),
    c(NA, "High Risk (M/WT)", result_table$High_RC_Num),
    c(NA, "Low Risk (M/WT)", result_table$Low_RC_Num),
    c(NA, "OR (95%CI)", OR.text),
    c(NA, "P Value", ifelse(result_table$pvalue < 0.001, "P < 0.001", round(result_table$pvalue, 3)))
  )
  pdf(file = paste0(mutation_graph_output, "/", dataset, "_Fisher_Test.pdf"), width = 3.25, height = 2, onefile = F)
  print(forestplot(
    labeltext = table.text,
    graph.pos = "right",
    align = "c",
    col = fpColors(box = rgb(247, 170, 088, maxColorValue = 255), lines = rgb(231, 098, 084, maxColorValue = 255), zero = "gray50"),
    mean = c(NA, NA, result_table$odds_ratio),
    lower = c(NA, NA, result_table$lower_conf_int),
    upper = c(NA, NA, result_table$upper_conf_int),
    boxsize = 0.1, lwd.ci = 1,
    ci.vertices.height = 0.08, ci.vertices = TRUE,
    zero = 1, lwd.zero = 1,
    colgap = unit(0.5, "mm"),
    lwd.xaxis = 1,
    lineheight = unit(0.4, "cm"),
    graphwidth = unit(0.15, "npc"),
    cex = 0.9, fn.ci_norm = fpDrawCircleCI,
    hrzl_lines = list(
      "2" = gpar(lwd = 2, col = "black"),
      "3" = gpar(lwd = 2, col = "black"),
      "12" = gpar(lwd = 2, col = "black")
    ),
    clip = c(0, 5),
    txt_gp = fpTxtGp(
      label = gpar(cex = 0.5, fontfamily = "Arial"),
      ticks = gpar(cex = 0.5, fontfamily = "Arial"),
      xlab = gpar(cex = 0.5, fontfamily = "Arial"),
      title = gpar(cex = 0.5, fontfamily = "Arial")
    ),
    xlab = ""
  ))
  grid.text(dataset,
    y = unit(0.95, "npc"),
    gp = gpar(fontsize = 8, fontfamily = "Arial", fontface = "bold")
  )
  dev.off()
  # somaticInteractions
  pdf(
    file = paste0(somatic_interaction_graph_output, "/", dataset, "_High_Risk_somaticInteraction.pdf"),
    width = 6, height = 6, onefile = F, fonts = "Arial"
  )
  interact_high <- somaticInteractions(
    maf = mutation_data_high_risk, top = 25,
    pvalue = c(0.05, 0.2),
    sigSymbolsSize = 1,
    fontSize = 0.5, countsFontSize = 0.5,
  )
  dev.off()
  pdf(
    file = paste0(somatic_interaction_graph_output, "/", dataset, "_Low_Risk_somaticInteraction.pdf"),
    width = 6, height = 6, onefile = F
  )
  interact_low <- somaticInteractions(
    maf = mutation_data_low_risk, top = 25, pvalue = c(0.05, 0.2),
    fontSize = 0.5, countsFontSize = 0.5
  )
  pdf(
    file = paste0(somatic_interaction_graph_output, "/", dataset, "_High_Risk_somaticInteraction_small.pdf"),
    width = 3.25, height = 3.25, onefile = F, fonts = "Arial"
  )
  interact_high <- somaticInteractions(
    maf = mutation_data_high_risk, top = 25,
    pvalue = c(0.05, 0.2),
    sigSymbolsSize = 1,
    fontSize = 0.5, countsFontSize = 0.5,
  )
  dev.off()
  pdf(
    file = paste0(somatic_interaction_graph_output, "/", dataset, "_Low_Risk_somaticInteraction_small.pdf"),
    width = 3.25, height = 3.25, onefile = F
  )
  interact_low <- somaticInteractions(
    maf = mutation_data_low_risk, top = 25, pvalue = c(0.05, 0.2),
    sigSymbolsSize = 1,
    fontSize = 0.5, countsFontSize = 0.5
  )
  dev.off()
  # tmb
  tmb_high <- tmb(maf = mutation_data_high_risk, captureSize = 50, logScale = T)
  tmb_low <- tmb(maf = mutation_data_low_risk, captureSize = 50, logScale = T)
  tmb <- c(tmb_high$total_perMB_log, tmb_low$total_perMB_log)
  group <- c(rep("High", length(tmb_high$total_perMB_log)), rep("Low", length(tmb_low$total_perMB_log)))
  tmb_table <- data.frame(TMB = tmb, Risk_Score = group)
  comparasions <- list(comparation = c("Low", "High"))
  ggplot(tmb_table, aes(x = Risk_Score, y = TMB, fill = Risk_Score)) +
    geom_violin(trim = F, color = "white") +
    geom_boxplot(width = 0.2, position = position_dodge(0.9)) +
    scale_fill_manual(values = c(
      "High" = rgb(231, 098, 084, maxColorValue = 255),
      "Low" = rgb(114, 188, 213, maxColorValue = 255)
    )) +
    labs(title = dataset, y = "TMB (log)", x = "Risk Score") +
    stat_compare_means(comparisons = comparasions, method = "t.test", hide.ns = F, label = "p.format", size = 2) +
    my_plot_theme(legend = "none")
  ggsave(filename = paste0(tmb_graph_output, "/", dataset, "_TMB.pdf"), width = 1.625, height = 2, device = "pdf")
  # MATH score
  barcode_high <- unique(mutation_data_high_risk@data$Tumor_Sample_Barcode)
  mutation_data_high_risk@data$VAF <- mutation_data_high_risk@data$t_alt_count / mutation_data_high_risk@data$t_depth
  math_high <- data.frame()
  for (i in barcode_high) {
    out.math <- inferHeterogeneity(maf = mutation_data_high_risk, tsb = i, vafCol = "VAF")
    Tumor_Sample_Barcode <- unique(out.math$clusterData$Tumor_Sample_Barcode)
    math_score <- unique(out.math$clusterData$MATH)
    out <- data.frame(Tumor_Sample_Barcode, math_score)
    math_high <- rbind(math_high, out)
  }
  barcode_low <- unique(mutation_data_low_risk@data$Tumor_Sample_Barcode)
  mutation_data_low_risk@data$VAF <- mutation_data_low_risk@data$t_alt_count / mutation_data_low_risk@data$t_depth
  math_low <- data.frame()
  for (i in barcode_low) {
    out.math <- inferHeterogeneity(maf = mutation_data_low_risk, tsb = i, vafCol = "VAF")
    Tumor_Sample_Barcode <- unique(out.math$clusterData$Tumor_Sample_Barcode)
    math_score <- unique(out.math$clusterData$MATH)
    out <- data.frame(Tumor_Sample_Barcode, math_score)
    math_low <- rbind(math_low, out)
  }
  math_all <- rbind(math_high, math_low)
  math_all$Risk_Score <- c(rep("High", nrow(math_high)), rep("Low", nrow(math_low)))
  ggplot(math_all, aes(x = Risk_Score, y = math_score, fill = Risk_Score)) +
    geom_violin(trim = F, color = "white") +
    geom_boxplot(width = 0.2, position = position_dodge(0.9)) +
    scale_fill_manual(values = c(
      "High" = rgb(231, 098, 084, maxColorValue = 255),
      "Low" = rgb(114, 188, 213, maxColorValue = 255)
    )) +
    labs(title = dataset, y = "MATH Score", x = "Risk Score") +
    stat_compare_means(comparisons = comparasions, method = "t.test", hide.ns = F, label = "p.format", size = 2) +
    my_plot_theme(legend = "none")
  ggsave(filename = paste0(math_graph_output, "/", dataset, "_MATH.pdf"), width = 1.625, height = 2, device = "pdf")
}
# MSI分析
library(cBioPortalData)
library(rapiclient)
library(AnVIL)
cbio <- cBioPortal()
studies <- getStudies(cbio)
studies$studyId
id <- "kirc_tcga_pan_can_atlas_2018"
clinical_data_cBio <- clinicalData(cbio, id)
clinical <- qread(file = paste0(clinical_dir, "/TCGA.qs"))
clinical$CaseID <- str_replace_all(clinical$CaseID, pattern = fixed("."), replacement = "-")
clinical$MSI <- as.numeric(clinical_data_cBio$MSI_SCORE_MANTIS[match(clinical$CaseID, clinical_data_cBio$patientId)])
comparisons <- list(comparison = c("Low", "High"))
ggplot(clinical, aes(x = risk_score_group, y = MSI, fill = risk_score_group)) +
  geom_violin(trim = F, color = "white") +
  geom_boxplot(width = 0.2, position = position_dodge(0.9)) +
  scale_fill_manual(values = c(
    "High" = rgb(231, 098, 084, maxColorValue = 255),
    "Low" = rgb(114, 188, 213, maxColorValue = 255)
  )) +
  labs(title = "Microsatellite Instability", y = "MSI Score", x = "Risk Score") +
  stat_compare_means(
    comparisons = comparisons,
    method = "t.test",
    hide.ns = F,
    label = "p.format",
    size = 2
  ) +
  my_plot_theme(legend = "none")
ggsave(filename = paste0(output_dir, "/MSI.pdf"), width = 1.625, height = 2, device = "pdf")
