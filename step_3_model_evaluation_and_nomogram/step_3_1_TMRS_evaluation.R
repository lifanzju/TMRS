library(qs)
library(survival)
library(survminer)
library(dplyr)
library(stringr)
library(patchwork)
library(timeROC)
library(GSVA)
library(ggstatsplot)
library(data.table)

rm(list = ls())
source("scripts/utils/plot_settings.R")
exp_dir <- "data/data_for_ml/exp"
clinical_dir <- "data/data_for_ml/clinical"
risk_scores_dir <- "results/step_2_model_construction/machine_learning/risk_score"
risk_score_genes <- read.csv(
  file = "results/step_2_model_construction/machine_learning/variables_of_final_model.csv",
  header = T, row.names = NULL
)[, 1]
result_dir <- "results/step_3_model_evaluation_and_nomogram/TMRS_performance_evaluation"
if (!dir.exists(result_dir)) dir.create(result_dir)
processed_clinical_output_dir <- "results/step_3_model_evaluation_and_nomogram/processed_clinical_data"
if (!dir.exists(processed_clinical_output_dir)) dir.create(processed_clinical_output_dir)
datasets <- str_remove_all(list.files(exp_dir), pattern = fixed(".csv"))
dataset_names <- list(
  "TCGA" = "TCGA",
  "E_MTAB_1980" = "E-MTAB-1980",
  "CPTAC" = "CPTAC",
  "TCGA_training" = "TCGA Training",
  "TCGA_validation" = "TCGA Validation"
)

# clinical data preprocessing
for (dataset in datasets) {
  clinical <- read.csv(file = paste0(clinical_dir, "/", dataset, ".csv"), row.names = NULL)
  if (dataset == "TCGA") {
    training <- read.csv(paste0(risk_scores_dir, "/", dataset, "_training_risk_scores.csv"), row.names = NULL)
    validation <- read.csv(paste0(risk_scores_dir, "/", dataset, "_validation_risk_scores.csv"), row.names = NULL)
    risk_scores <- rbind(training, validation)
  } else {
    risk_scores <- read.csv(paste0(risk_scores_dir, "/", dataset, "_risk_scores.csv"), row.names = NULL)
  }
  colnames(risk_scores)[1] <- "CaseID"
  risk_scores <- risk_scores[, c("CaseID", "risk_score")]
  clinical <- merge(clinical, risk_scores, by = "CaseID", all = T)
  clinical <- arrange(clinical, risk_score)
  clinical$Age <- as.numeric(clinical$Age)
  clinical$risk_score_group <- ifelse(clinical$risk_score >= median(clinical$risk_score), "High", "Low")
  clinical$stage_group <- ifelse(clinical$Stage %in% c("S1", "S2"), "Low",
    ifelse(clinical$Stage %in% c("S3", "S4"), "High", NA)
  )
  clinical$grade_group <- ifelse(clinical$Grade %in% c("G1", "G2"), "Low",
    ifelse(clinical$Grade %in% c("G3", "G4"), "High", NA)
  )
  clinical$risk_score_group <- factor(clinical$risk_score_group, levels = c("Low", "High"), ordered = T)
  clinical$stage_group <- factor(clinical$stage_group,
    levels = c("Low", "High"),
    labels = c("I/II", "III/IV"), ordered = T
  )
  clinical$grade_group <- factor(clinical$grade_group,
    levels = c("Low", "High"),
    labels = c("I/II", "III/IV"), ordered = T
  )
  qsave(clinical, file = paste0(processed_clinical_output_dir, "/", dataset, ".qs"))
}

# assoication with tryptophan metabolism score
cor_with_ssgsea_dir <- file.path(result_dir, "cor_with_ssgsea")
if (!dir.exists(cor_with_ssgsea_dir)) dir.create(cor_with_ssgsea_dir)
tryptophan_metabolism_genes <- read.delim(file = "data/genesets/tryptophan_metabolism.txt", header = F, row.names = NULL)[, 1]
genesets <- list(tryptophan_metabolism_genes)
for (dataset in datasets) {
  exp <- fread(file = paste0(exp_dir, "/", dataset, ".csv")) %>% data.frame(row.names = 1)
  clinical <- qread(file = paste0(processed_clinical_output_dir, "/", dataset, ".qs"))
  gsva_score <- gsva(expr = as.matrix(exp), gset.idx.list = genesets, method = "ssgsea")
  gsva_score <- gsva_score[, clinical$CaseID] %>% matrix(ncol = 1)
  colnames(gsva_score) <- "Tryptophan_Metabolism"
  rownames(gsva_score) <- clinical$CaseID
  cor_analysis_table <- cbind(gsva_score, RS = clinical$risk_score) %>% data.frame()
  cor_result <- cor.test(cor_analysis_table$RS, cor_analysis_table$Tryptophan_Metabolism)
  cor_coeff <- cor_result$estimate
  p_value <- cor_result$p.value
  p_value_to_show <- ifelse(p_value < 1e-16, 
                            "P-value < 1e-16", 
                            paste0("P-value = ",base::format.pval(p_value,digits = 3)))
  scatter_plot <- ggplot(cor_analysis_table, aes(x = RS, y = Tryptophan_Metabolism)) +
    geom_point(size = 1) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "#4DBBD5FF", size = 1) +  
    labs(
      x = "Risk Score",
      y = "Tryptophan Metabolism Score",
      title = paste0(dataset_names[[dataset]]),
      subtitle = paste("Pearson's r = ", round(cor_coeff, 2), "\n",p_value_to_show))  + 
    my_plot_theme() +
    theme(plot.subtitle = element_text(hjust = 0.5))
  ggsave(
    filename = paste0(cor_with_ssgsea_dir, "/", dataset, "_rs_ssgsea_cor.pdf"), plot = scatter_plot,
    width = 2, height = 2, device = "pdf"
  )
  qsave(scatter_plot, file = paste0(cor_with_ssgsea_dir, "/", dataset, "_rs_ssgsea_cor.qs"))
}


# KM plot
km_plot_dir <- file.path(result_dir, "km_plot")
if (!dir.exists(km_plot_dir)) dir.create(km_plot_dir)
for (dataset in datasets) {
  clinical <- qread(file = paste0(processed_clinical_output_dir, "/", dataset, ".qs"))
  k_m_model <- survfit(Surv(OS, OS.Censor) ~ risk_score_group, data = clinical)
  k_m_plot <- ggsurvplot(k_m_model,
    data = clinical,
    title = dataset_names[dataset],
    pval = TRUE,
    legend.title = "",
    legend.labs = c(
      paste0("High TMRS"),
      paste0("Low TMRS")
    ),
    palette = c("firebrick2", "dodgerblue3"),
    risk.table = F,
    cumevents = F,
    pval.size = 2.5,
    censor.size = 2.5,
    xlab = "Survival(Days)",
    ylab = "Overall Survival Probability",
    break.time.by = ifelse(dataset == "CPTAC", 500, 1000)
  )
  km_plot <- k_m_plot$plot + 
    my_plot_theme(legend = c(0.8,1) ) +
    theme(legend.key.size = unit(0.1, "cm"))
  qsave(km_plot, file = paste0(km_plot_dir, "/", dataset, "_km_plot.qs"))
  ggsave(paste0(km_plot_dir, "/", dataset, "_km_plot.pdf"), width = 2, height = 2)
}

# ROC plot
roc_curve_dir <- file.path(result_dir, "roc_curve")
if (!dir.exists(roc_curve_dir)) dir.create(roc_curve_dir)
for (dataset in datasets) {
  clinical <- qread(file = paste0(processed_clinical_output_dir, "/", dataset, ".qs"))
  clinical$OS_year <- clinical$OS / 365
  tROC <- timeROC(
    T = clinical$OS_year, delta = clinical$OS.Censor, marker = clinical$risk_score,
    cause = 1, times = c(1, 3, 5), ROC = T
  )
  roc_data <- data.frame(
    FPR = c(tROC$FP[, 1], tROC$FP[, 2], tROC$FP[, 3]),
    TPR = c(tROC$TP[, 1], tROC$TP[, 2], tROC$TP[, 3]),
    Year = factor(rep(c("1 year", "3 years", "5 years"), each = length(tROC$FP[, 1])))
  )

  colors <- c(
    "1 year" = rgb(114, 188, 213, maxColorValue = 255),
    "3 years" = rgb(255, 208, 111, maxColorValue = 255),
    "5 years" = rgb(231, 98, 84, maxColorValue = 255)
  )

  auc_labels <- c(
    paste0("1 year (AUC = ", round(tROC$AUC[1], 2), ")"),
    paste0("3 years (AUC = ", round(tROC$AUC[2], 2), ")"),
    paste0("5 years (AUC = ", round(tROC$AUC[3], 2), ")")
  )


  roc_curve <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Year)) +
    geom_line(size = 1.2) +
    scale_color_manual(values = colors, labels = auc_labels) +
    labs(
      title = dataset_names[dataset],
      x = "1 - Specificity",
      y = "Sensitivity",
      color = ""
    ) +
    my_plot_theme(legend = c(0.6, 0.3)) +
    theme(legend.key.height  = unit(0.3, "cm"))
  qsave(roc_curve, file = paste0(roc_curve_dir, "/", dataset, "_roc_curve.qs"))
  ggsave(paste0(roc_curve_dir, "/", dataset, "_roc_curve.pdf"), width = 2, height = 2)
}
