library(dplyr)
library(stringr)
library(ggplot2)
library(survival)
library(survminer)
library(ggstatsplot)
library(qs)

rm(list=ls())
source("scripts/utils/plot_settings.R")
# TIDE
tide_dir <- "data/TIDE"
clinical_dir <- "results/step_3_model_evaluation_and_nomogram/processed_clinical_data"
output_dir <- "results/step_6_therapy_response/TIDE"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

datasets <- c("TCGA","CPTAC")
for (dataset in datasets) {
  tide_score <- read.csv(file.path(tide_dir, paste0(dataset, ".csv")))
  clinical <- qread(file.path(clinical_dir, paste0(dataset, ".qs")))
  clinical <- merge(clinical,tide_score,by.x = "CaseID",by.y = "Patient")
  clinical <- arrange(clinical,risk_score)
  for (var in c("TIDE","IFNG","Dysfunction","Exclusion")) {
    cor_results <- cor.test(clinical$risk_score,clinical[[var]])
    p_value <- cor_results$p.value
    cor_coef <- cor_results$estimate
    p_value_to_show <- ifelse(p_value < 1e-16, 
                              "P-value < 1e-16", 
                              paste0("P-value = ",base::format.pval(p_value,digits = 3)))
    ggplot(clinical, aes_string(x = "risk_score", y = var)) +
      geom_point(size=1) +
      geom_smooth(method = "lm", se = T, color = "#4DBBD5FF") +
      labs(title = paste0(dataset,"-",var),
           x = "TMRS",
           y = var,
           subtitle = paste("Pearson's r = ", round(cor_coef, 2), "\n",p_value_to_show))+
      my_plot_theme() +
      theme(plot.subtitle = element_text(hjust = 0.5))
    ggsave(file.path(output_dir, paste0(dataset, "_", var, ".pdf")), width = 1.625, height = 2)
  }
}



