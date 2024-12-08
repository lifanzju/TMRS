library(stringr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(qs)
library(forestplot)

rm(list = ls())
source("scripts/utils/plot_settings.R")
output_dir <- "results/step_2_model_construction/validation_of_model_genes"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = T)
}

# expression in tumor and normal samples
exp_dir <- "data/datasets/EXP"
sample_dir <- "data/datasets/SampleType"
genes_in_model <- read.csv("results/step_2_model_construction/machine_learning/variables_of_final_model.csv")[, 1]
datasets <- c("TCGA", "CPTAC", "GSE53757", "CPTAC_proteome")
tvn_output_dir <- file.path(output_dir, "tumor_vs_normal")
if (!dir.exists(tvn_output_dir)) {
  dir.create(tvn_output_dir, recursive = T)
}
for (dataset in datasets) {
  plot_title <- dataset
  if (dataset == "GSE53757") {
    exp <- read.delim("data/GSE53757/GSE53757.txt",
      row.names = 1, header = T
    )
    sampletype <- read.delim("data/GSE53757/sample_type.txt",
      row.names = 1, header = T
    )
  } else if (dataset == "CPTAC_proteome") {
    exp <- read.delim("data/CPTAC_proteome/CPTAC_proteome.txt",
      row.names = 1, header = T
    )
    sampletype <- read.delim("data/CPTAC_proteome/CPTAC_proteome_SampleType.txt",
      row.names = 1, header = T
    )
    plot_title <- "CPTAC Proteome"
  } else {
    exp <- read.delim(paste(exp_dir, paste0(dataset, ".txt"), sep = "/"),
      row.names = 1, header = T
    )
    sampletype <- read.delim(paste(sample_dir, paste0(dataset, "_SampleType.txt"), sep = "/"),
      row.names = 1, header = T
    )
  }
  exp <- exp[, row.names(sampletype)]
  gene_exp <- exp[genes_in_model, ] %>%
    t() %>%
    as.data.frame()
  analysis_table <- cbind(Type = sampletype$Type, gene_exp)
  analysis_table <- melt(analysis_table, id.vars = "Type")
  analysis_table <- na.omit(analysis_table)
  ggplot(analysis_table, aes(x = variable, y = value, fill = Type)) +
    geom_violin(trim = F, color = "white", scale = "width", adjust = 3 / 4) +
    geom_boxplot(
      width = 0.3,
      position = position_dodge(0.9),
      linewidth = 0.3,
      outlier.size = 0.3
    ) +
    scale_fill_manual(values = c(
      rgb(114, 188, 213, maxColorValue = 255),
      rgb(231, 098, 084, maxColorValue = 255)
    )) +
    labs(x = NULL, y = ifelse(dataset == "CPTAC_proteome","Protein Abundance","Gene Expression"), 
         fill = "Sample", title = plot_title) +
    stat_compare_means(method = "t.test", hide.ns = F, label = "p.signif", size = 2) +
    my_plot_theme(x.text.angle = 45, legend = "top") +
    theme(
      legend.key.size = unit(0.3, "cm"),
      plot.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(-5, 0, -15, 0)
    )
  ggsave(file.path(tvn_output_dir, paste0(dataset, "_tumor_vs_normal.pdf")), width = 3.25, height = 2)
}

# association with clinical features
exp_dir <- "data/Datasets/proprocessed_exp"
clinical_dir <- "results/step_3_model_evaluation_and_nomogram/processed_clinical_data"
grade_output_dir <- file.path(output_dir, "clinical_features", "grade")
stage_output_dir <- file.path(output_dir, "clinical_features", "stage")
if (!dir.exists(grade_output_dir)) {
  dir.create(grade_output_dir, recursive = T)
}
if (!dir.exists(stage_output_dir)) {
  dir.create(stage_output_dir, recursive = T)
}
datasets <- c("TCGA", "CPTAC", "E_MTAB_1980", "CPTAC_proteome")
for (dataset in datasets) {
  if (dataset == "CPTAC_proteome") {
    exp <- read.delim("data/CPTAC_proteome/CPTAC_proteome.txt", row.names = 1)
    sampletype <- read.delim("data/CPTAC_proteome/CPTAC_proteome_SampleType.txt", row.names = 1)
    exp <- exp[, row.names(sampletype)[sampletype$Type == "Tumor"]]
    colnames(exp) <- str_split(colnames(exp), fixed("."), simplify = T)[, 1:2] %>% 
      apply(1, paste, collapse = ".")
    clinical <- qread(file.path(clinical_dir, "CPTAC.qs"))
    y_limits <- c(-4,3)
    plot_title <- "CPTAC Proteome"
  } else {
    exp <- read.csv(file.path(exp_dir, paste0(dataset, ".csv")), row.names = 1)
    clinical <- qread(file.path(clinical_dir, paste0(dataset, ".qs")))
    if(dataset == "E_MTAB_1980") {
      y_limits <- c(0, 20)
      plot_title <- "E-MTAB-1980"
    }else{
      y_limits <- c(-10, 15)
      plot_title <- dataset
      }
  }
  exp <- exp[, clinical$CaseID]
  gene_exp <- exp[genes_in_model, ] %>%
    t() %>%
    as.data.frame()
  for (variable in c("Grade", "Stage")) {
    analysis_table <- cbind(select(clinical, all_of(variable)), gene_exp)
    analysis_table <- melt(analysis_table, id.vars = variable, variable.name = "Gene", value.name = "Exp") %>%
      na.omit()
    if (variable == "Grade") analysis_table <- subset(analysis_table, Grade != "GX")
    analysis_table[, variable] <- as.factor(analysis_table[, variable])
    ggplot(analysis_table, aes(x = Gene, y = Exp, fill = get(variable))) +
      geom_violin(trim = F, color = "white", scale = "width", adjust = 3 / 4) +
      geom_boxplot(
        width = 0.3,
        position = position_dodge(0.9),
        linewidth = 0.3,
        outlier.size = 0.3
      ) +
      scale_fill_manual(values = c(
        rgb(055, 103, 149, maxColorValue = 255),
        rgb(114, 188, 213, maxColorValue = 255),
        rgb(255, 208, 111, maxColorValue = 255),
        rgb(231, 098, 084, maxColorValue = 255)
      )) +
      scale_y_continuous(limits = y_limits) +
      labs(x = "", y = paste0("Gene Expression"), fill = variable, title = plot_title) +
      stat_compare_means(method = "anova", hide.ns = F, label = "p.signif", size = 2) +
      my_plot_theme(legend = "top") +
      theme(
        legend.key.size = unit(0.3, "cm"),
        plot.margin = margin(0, 0, -5, 0),
        legend.box.margin = margin(-5, 0, -15, 0)
      )
    ggsave(filename = paste0(
      ifelse(variable == "Grade", grade_output_dir, stage_output_dir),
      "/", dataset, "_", variable, ".pdf"
    ), width = 6.5, height = 2)
  }
}

# multivariate cox regression
multivariate_cox_output_dir <- file.path(output_dir, "multivariate_cox")
if (!dir.exists(multivariate_cox_output_dir)) {
  dir.create(multivariate_cox_output_dir, recursive = T)
}
clinical_dir <- 'results/step_3_model_evaluation_and_nomogram/processed_clinical_data'
exp_dir <- 'data/data_for_ml/exp'
dataset_names <- list(
  "TCGA" = "TCGA",
  "CPTAC" = "CPTAC",
  "E_MTAB_1980" = "E-MTAB-1980"
)
for (dataset in c('TCGA','CPTAC','E_MTAB_1980')){
  dataset_output_dir <- file.path(multivariate_cox_output_dir, dataset)
  if (!dir.exists(dataset_output_dir)) {
    dir.create(dataset_output_dir, recursive = T)
  }
  clinical_data <- qread(file.path(clinical_dir,paste0(dataset,'.qs')))
  exp <- read.csv(file.path(exp_dir, paste0(dataset, ".csv")), row.names = 1)
  exp <- t(exp)
  exp <- exp[,genes_in_model]
  clinical_data <- merge(clinical_data,exp,by.x = "CaseID",by.y = 'row.names')
  clinical_data <- subset(clinical_data,grade_group %in% c("I/II","III/IV") & stage_group %in% c("I/II","III/IV"))
  for (gene in genes_in_model){
    coxresult <- summary(coxph(Surv(OS, OS.Censor) ~ get(gene) + Age + grade_group +
                                 stage_group, data = clinical_data))
    coef <- coxresult$coefficients[, "coef"]
    HR <- coxresult$conf.int[, "exp(coef)"]
    HR.95L <- coxresult$conf.int[, "lower .95"]
    HR.95H <- coxresult$conf.int[, "upper .95"]
    pvalue <- coxresult$coefficients[, "Pr(>|z|)"]
    result_table <- data.frame(
      coef = coef,
      HR = HR,
      HR.95L = HR.95L,
      HR.95H = HR.95H,
      pvalue = pvalue
    )
    row.names(result_table) <- c(
      gene,
      "Age (Years)",
      "Grade (G3&G4/G1&G2)",
      "Stage (III&IV/I&II)"
    )
    result_table <- as.data.frame(cbind(id = row.names(result_table), result_table))
    HR_text <- paste(
      round(result_table$HR, 3),
      "(",
      round(result_table$HR.95L, 3),
      "-",
      round(result_table$HR.95H, 3),
      ")",
      sep = ""
    )
    tabletext <- cbind(
      c(NA, "Variable", result_table$id),
      c(NA, "Coefficient", round(result_table$coef, 3)),
      c(
        NA,
        "P value",
        ifelse(
          result_table$pvalue < 0.001,
          format(result_table$pvalue, method = "e", digits = 2),
          round(result_table$pvalue, 3)
        )
      ),
      c(NA, "HR(95% CI)", HR_text)
    )
    pdf(file.path(dataset_output_dir, 
                  paste0(gene,'_',dataset,'_multi_cox_regression.pdf')),
        width = 3.25,height = 1.33,onefile = FALSE)
    forestplot(
      labeltext = tabletext,
      graph.pos = "right", # Position of the plot (p-value column)
      col = fpColors(box = "#D55E00", lines = "#CC79A7", zero = "gray50"),
      mean = c(NA, NA, result_table$HR),
      lower = c(NA, NA, result_table$HR.95L), # 95% CI lower bounds
      upper = c(NA, NA, result_table$HR.95H), # 95% CI upper bounds
      boxsize = 0.1, lwd.ci = 1, # Box size and CI line width
      ci.vertices.height = 0.08, ci.vertices = TRUE, # Confidence interval end markers
      zero = 1, lwd.zero = 1, # Reference line at HR = 1
      colgap = unit(1.5, "mm"), # Column gap for better spacing
      lwd.xaxis = 1, # X-axis line width
      lineheight = unit(0.45, "cm"), # Fixed row height
      graphwidth = unit(0.2, "npc"), # Graph width proportion
      cex = 0.9, fn.ci_norm = fpDrawCircleCI, # Circle CI points
      
      # Horizontal lines in the plot
      hrzl_lines = list(
        "2" = gpar(lwd = 2, col = "black"),
        "3" = gpar(lwd = 2, col = "black"), # Line after the header row
        "7" = gpar(lwd = 2, col = "black") # Bottom line
      ),
      
      # Font sizes and font family
      txt_gp = fpTxtGp(
        label = gpar(cex = 0.5, fontfamily = "Arial"),
        ticks = gpar(cex = 0.5, fontfamily = "Arial"),
        xlab = gpar(cex = 0.5, fontfamily = "Arial"),
        title = gpar(cex = 0.5, fontfamily = "Arial")
      ),
    ) %>% print()
    grid.text(paste0(dataset_names[[dataset]]," - ",gene),
              y = unit(0.9, "npc"),
              gp = gpar(fontsize = 8, fontfamily = "Arial", fontface = "bold")
    )
    dev.off()
  }
}



