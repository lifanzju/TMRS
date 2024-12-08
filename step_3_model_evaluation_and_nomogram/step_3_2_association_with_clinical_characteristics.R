library(qs)
library(data.table)
library(survival)
library(survminer)
library(dplyr)
library(stringr)
library(patchwork)
library(ComplexHeatmap)
library(circlize)
library(forestplot)
library(pROC)
library(gbm)

rm(list = ls())
source("scripts/utils/plot_settings.R")
exp_dir <- "data/data_for_ml/exp"
clinical_dir <- "results/step_3_model_evaluation_and_nomogram/processed_clinical_data"
risk_score_genes <- read.csv(
  file = "results/step_2_model_construction/machine_learning/variables_of_final_model.csv",
  header = T, row.names = NULL
)[, 1]
result_dir <- "results/step_3_model_evaluation_and_nomogram/association_with_clinical_characteristics"
if (!dir.exists(result_dir)) dir.create(result_dir)
datasets <- str_remove_all(list.files(exp_dir), pattern = fixed(".csv"))
dataset_names <- list(
  "TCGA" = "TCGA",
  "E_MTAB_1980" = "E-MTAB-1980",
  "CPTAC" = "CPTAC",
  "TCGA_training" = "TCGA Training",
  "TCGA_validation" = "TCGA Validation"
)

# clinical characteristics heatmap
clinical_heatmap_dir <- file.path(result_dir, "clinical_heatmap")
if (!dir.exists(clinical_heatmap_dir)) dir.create(clinical_heatmap_dir)
for (dataset in datasets) {
  exp <- fread(file = paste0(exp_dir, "/", dataset, ".csv")) %>% data.frame(row.names = 1)
  clinical <- qread(file = paste0(clinical_dir, "/", dataset, ".qs"))
  # main figure
  exp <- as.matrix(exp[risk_score_genes, clinical$CaseID])
  main_table_color <- colorRamp2(c(-2.5, -1, 0, 1, 2.5), c("dodgerblue2", "deepskyblue2", "white", "darkorange", "firebrick2"), transparency = 0.3)
  # top annotation
  risk_score_color <- colorRamp2(c(-2, -1, 2), c("yellow", "gold", "firebrick2"), transparency = 0.2)
  stage_color <- dataset_color <- c(
    "S4" = "coral1",
    "S3" = "orange",
    "S2" = "lightskyblue",
    "S1" = "deepskyblue3"
  )
  grade_color <- c(
    "G4" = "coral1",
    "G3" = "orange",
    "G2" = "lightskyblue",
    "G1" = "deepskyblue3",
    "GX" = "gray90"
  )
  sex_color <- c("Male" = "deepskyblue", "Female" = "hotpink1")
  age_color <- colorRamp2(c(20, 75, 90), c("honeydew", "seagreen3", "green4"), transparency = 0.2)
  survival_color <- ifelse(clinical$OS.Censor == 0, "dodgerblue2", "firebrick2")
  risk_score_anno <- anno_simple(clinical$risk_score, col = risk_score_color, border = T)
  stage_anno <- anno_simple(clinical$Stage, col = stage_color, border = T)
  grade_anno <- anno_simple(clinical$Grade, col = grade_color, border = T)
  age_anno <- anno_simple(clinical$Age, col = age_color, border = T)
  sex_anno <- anno_simple(clinical$Gender, col = sex_color, border = T)
  survival_anno <- anno_points(clinical$OS, gp = gpar(col = survival_color), height = unit(2, "cm"), size = unit(1, "mm"))
  Top <- columnAnnotation(
    "Overall Survival\n(Days)" = survival_anno,
    "Risk Score" = risk_score_anno,
    Age = age_anno,
    Sex = sex_anno,
    Stage = stage_anno,
    Grade = grade_anno,
    gap = unit(2, "points"),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 6, fontfamily = "Arial")
  )
  # prepare legends
  survival_legend <- Legend(
    at = c("Alive", "Dead"),
    legend_gp = gpar(col = c("Alive" = "dodgerblue2", "Dead" = "firebrick2")),
    labels_gp = gpar(fontsize = 6, fontfamily = "Arial"),  # 保持文字大小
    type = "points",
    title_gp = gpar(fontsize = 7, fontface = "plain", fontfamily = "Arial"),
    title = "Survival",
    grid_height = unit(0.1, "cm"),  
    grid_width = unit(0.1, "cm")    
  )
  risk_score_legend <- Legend(
    title = "Risk Score", col_fun = risk_score_color,
    labels_gp = gpar(fontsize = 6, fontfamily = "Arial"), 
    title_gp = gpar(fontsize = 7, fontface = "plain", fontfamily = "Arial"),
    legend_height = unit(0.4, "cm")
  )
  age_legend <- Legend(
    title = "Age", col_fun = age_color,
    labels_gp = gpar(fontsize = 6, fontfamily = "Arial"),  
    title_gp = gpar(fontsize = 7, fontface = "plain", fontfamily = "Arial"),
    legend_height = unit(0.4, "cm")  
  )
  sex_legend <- Legend(
    title = "Sex", at = c("Male", "Female"),
    legend_gp = gpar(fill = sex_color),
    labels_gp = gpar(fontsize = 6, fontfamily = "Arial"),  
    title_gp = gpar(fontsize = 7, fontface = "plain", fontfamily = "Arial"),
    grid_height = unit(0.1, "cm"), 
    grid_width = unit(0.3, "cm")    
  )
  stage_legend <- Legend(
    title = "Stage", at = names(stage_color),
    labels = c("IV", "III", "II", "I"),
    legend_gp = gpar(fill = stage_color),
    labels_gp = gpar(fontsize = 6, fontfamily = "Arial"),  
    title_gp = gpar(fontsize = 7, fontface = "plain", fontfamily = "Arial"),
    grid_height = unit(0.1, "cm"),
    grid_width = unit(0.3, "cm")
  )
  grade_legend <- Legend(
    title = "Grade", at = names(grade_color),
    legend_gp = gpar(fill = grade_color),
    labels_gp = gpar(fontsize = 6, fontfamily = "Arial"), 
    title_gp = gpar(fontsize = 7, fontface = "plain", fontfamily = "Arial"),
    grid_height = unit(0.1, "cm"),
    grid_width = unit(0.3, "cm")
  )
  exp_legend <- Legend(
    title = "Exp\n(Z-Scores)",
    col_fun = main_table_color,
    labels_gp = gpar(fontsize = 6, fontfamily = "Arial"), 
    title_gp = gpar(fontsize = 7, fontface = "plain", fontfamily = "Arial"),
    legend_height = unit(0.4, "cm")
  )
  
  legend_list <- packLegend(
    survival_legend, sex_legend, stage_legend,
    grade_legend, age_legend, risk_score_legend, exp_legend,
    gap = unit(1, "mm")  
  )
  
  ht <- Heatmap(
    exp,
    cluster_columns = FALSE, cluster_rows = FALSE, border = TRUE,
    row_names_side = "left", col = main_table_color, na_col = "white",
    show_column_names = FALSE, top_annotation = Top,
    row_names_gp = gpar(fontsize = 6, fontfamily = "Arial"),  # 保持文字大小
    show_heatmap_legend = FALSE
  )

  pdf(paste0(clinical_heatmap_dir, "/", dataset, "_clinical_heatmap.pdf"), width = 3.25, height = 6, family = "Arial")
  full_ht <- draw(ht, annotation_legend_list = legend_list, padding = unit(c(0.1, 0, 0.5, 0), "cm"))
  grid.text(dataset_names[[dataset]], 
            x = unit(0.5, "npc"), y = unit(0.98, "npc"),
            gp = gpar(fontsize = 8, fontfamily = "Arial", fontface = "bold"))
  dev.off()
  qsave(full_ht, file = paste0(clinical_heatmap_dir, "/", dataset, "_clinical_heatmap.qs"))
}

# risk score of different clinical characteristics
clinical_violin_dir <- file.path(result_dir, "clinical_violin")
stage_graph_dir <- paste0(clinical_violin_dir, "/", "stage")
grade_graph_dir <- paste0(clinical_violin_dir, "/", "grade")
if (!dir.exists(stage_graph_dir)) dir.create(stage_graph_dir, recursive = T)
if (!dir.exists(grade_graph_dir)) dir.create(grade_graph_dir, recursive = T)
for (dataset in datasets) {
  clinical <- qread(file = paste0(clinical_dir, "/", dataset, ".qs"))
  stage_analysis_table <- subset(clinical, Stage %in% c("S1", "S2", "S3", "S4"))
  grade_analysis_table <- subset(clinical, Grade %in% c("G1", "G2", "G3", "G4"))
  stage_comparasions <- list()
  for (i in 1:(length(unique(stage_analysis_table$Stage)) - 1)) {
    stage_comparasions[[i]] <- sort(unique(stage_analysis_table$Stage))[i:(i + 1)]
  }
  grade_comparasions <- list()
  for (i in 1:(length(unique(grade_analysis_table$Grade)) - 1)) {
    grade_comparasions[[i]] <- sort(unique(grade_analysis_table$Grade))[i:(i + 1)]
  }
  # stage violin plot
  stage_violin_plot <- ggplot(stage_analysis_table, aes(x = Stage, y = risk_score, fill = Stage)) +
    geom_violin(trim = F, color = "white") +
    geom_boxplot(width = 0.2, position = position_dodge(0.9)) +
    scale_fill_manual(values = c(
      rgb(055, 103, 149, maxColorValue = 255),
      rgb(114, 188, 213, maxColorValue = 255),
      rgb(255, 208, 111, maxColorValue = 255),
      rgb(239, 138, 071, maxColorValue = 255),
      rgb(231, 098, 084, maxColorValue = 255)
    )) +
    ylab("Risk Score") +
    xlab("Stage") +
    scale_x_discrete(labels = c("I", "II", "III", "IV")) +
    scale_y_continuous(limits = c(
      min(stage_analysis_table$risk_score) - 1,
      max(stage_analysis_table$risk_score) + 1.5
    )) +
    ggtitle(dataset_names[[dataset]]) +
    stat_compare_means(method = "anova", 
                       label.y = max(stage_analysis_table$risk_score) + 1.3, 
                       label.x.npc = "left",
                       size = 6, size.unit = "pt") +
    stat_compare_means(comparisons = stage_comparasions, 
                       method = "t.test", 
                       hide.ns = F, label = "p.format",
                       size = 2) +
    my_plot_theme(legend = "none")
  qsave(stage_violin_plot, file = paste0(stage_graph_dir, "/", dataset, "_stage_violin.qs"))
  ggsave(filename = paste0(stage_graph_dir, "/", dataset, "_stage_violin.pdf"), device = "pdf", width = 1.625, height = 2)
  # grade violin plot
  grade_violin_plot <- ggplot(grade_analysis_table, aes(x = Grade, y = risk_score, fill = Grade)) +
    geom_violin(trim = F, color = "white") +
    geom_boxplot(width = 0.2, position = position_dodge(0.9)) +
    scale_fill_manual(values = c(
      rgb(055, 103, 149, maxColorValue = 255),
      rgb(114, 188, 213, maxColorValue = 255),
      rgb(255, 208, 111, maxColorValue = 255),
      rgb(239, 138, 071, maxColorValue = 255),
      rgb(231, 098, 084, maxColorValue = 255)
    )) +
    ylab("Risk Score") +
    xlab("Grade") +
    scale_x_discrete(labels = c("I", "II", "III", "IV", "GX")) +
    scale_y_continuous(limits = c(
      min(grade_analysis_table$risk_score) - 1,
      max(grade_analysis_table$risk_score) + 1.5
    )) +
    ggtitle(dataset_names[[dataset]]) +
    stat_compare_means(method = "anova", 
                       label.y = max(grade_analysis_table$risk_score) + 1.3, 
                       label.x.npc = "left",
                       size = 6, size.unit = "pt") +
    stat_compare_means(comparisons = grade_comparasions, 
                       method = "t.test", 
                       hide.ns = F, label = "p.format",
                       size = 2) +
    my_plot_theme(legend = "none")
  qsave(grade_violin_plot, file = paste0(grade_graph_dir, "/", dataset, "_grade_violin.qs"))
  ggsave(filename = paste0(grade_graph_dir, "/", dataset, "_grade_violin.pdf"), device = "pdf", width = 1.625, height = 2)
}

# multivariate cox regression
multivariate_cox_dir <- file.path(result_dir, "multivariate_cox_regression")
if (!dir.exists(multivariate_cox_dir)) dir.create(multivariate_cox_dir)
for (dataset in datasets) {
  clinical <- qread(file = paste0(clinical_dir, "/", dataset, ".qs"))
  clinical <- na.omit(clinical, cols = c("risk_score", "Age", "grade_group", "stage_group"))
  # fit the Cox model
  coxresult <- summary(coxph(Surv(OS, OS.Censor) ~ risk_score + Age + grade_group + stage_group, data = clinical))

  # extract coefficients, hazard ratios, and p-values
  coef <- coxresult$coefficients[, "coef"]
  HR <- coxresult$conf.int[, "exp(coef)"]
  HR.95L <- coxresult$conf.int[, "lower .95"]
  HR.95H <- coxresult$conf.int[, "upper .95"]
  pvalue <- coxresult$coefficients[, "Pr(>|z|)"]
  result_table <- data.frame(
    Variable = c("Risk Score", "Age(Years)", "Grade(III&IV/I&II)", "Stage(III&IV/I&II)"),
    coef = round(coef, 2),
    HR = round(HR, 2),
    HR.95L = round(HR.95L, 2),
    HR.95H = round(HR.95H, ),
    pvalue = ifelse(pvalue < 0.001, "<0.001", round(pvalue, 3))
  )

  # prepare table text for the forestplot function
  HR_text <- paste0(result_table$HR, " (", result_table$HR.95L, "-", result_table$HR.95H, ")")
  tabletext <- cbind(
    c(NA, "Variable", result_table$Variable),
    c(NA, "HR (95% CI)", HR_text),
    c(NA, "P Value", as.character(result_table$pvalue))
  )

  # generate the forest plot
  pdf(
    file = paste0(multivariate_cox_dir, "/", dataset, "_Forest_Plot_with_Table.pdf"),
    width = 3.25, height = 2, onefile = F, family = "Arial"
  )
  # draw the forest plot
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
    colgap = unit(3.5, "mm"), # Column gap for better spacing
    lwd.xaxis = 1, # X-axis line width
    lineheight = unit(0.75, "cm"), # Fixed row height
    graphwidth = unit(0.4, "npc"), # Graph width proportion
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
  grid.text(dataset_names[[dataset]],
    y = unit(0.95, "npc"),
    gp = gpar(fontsize = 8, fontfamily = "Arial", fontface = "bold")
  )
  dev.off()
}

# ability to distinguish normal and tumor samples (diagnostic value)
diag_value_dir <- file.path(result_dir, "diagnostic_value")
if (!dir.exists(diag_value_dir)) dir.create(diag_value_dir)
all_sample_exp_dir <- "data/Datasets/EXP"
sample_type_dir <- "data/Datasets/SampleType"
TMRS_model <- qread("results/step_2_model_construction/machine_learning/final_model.qs")
for (dataset in c("TCGA","CPTAC")) {
  exp <- fread(file = paste0(all_sample_exp_dir, "/", dataset, ".txt")) %>%
    data.frame(row.names = 1) %>%
    t()
  sample_type <- read.delim(file = paste0(sample_type_dir, "/", dataset, "_SampleType.txt"), row.names = 1)
  exp <- exp[rownames(sample_type), ]
  exp <- as.data.frame(scale(exp))
  risk_score <- as.numeric(predict(TMRS_model, newdata = exp, type = "link"))
  analysis_table <- data.frame(sample_type = sample_type$Type, risk_score = risk_score)
  diag_roc <- roc(analysis_table$sample_type, analysis_table$risk_score)
  diag_auc <- auc(diag_roc)
  diag_roc_data <- data.frame(
    FPR = 1 - diag_roc$specificities,
    TPR = diag_roc$sensitivities
  )
  best_threshold <- coords(diag_roc, "best", ret = c("threshold", "specificity", "sensitivity"))
  best_point <- data.frame(
    FPR = 1 - best_threshold[1, 2],
    TPR = best_threshold[1, 3]
  )
  best_label <- paste0(
    "Best Threshold: ", round(best_threshold["threshold"], 2), "\n",
    "Specificity: ", round(best_threshold["specificity"], 2), "\n",
    "Sensitivity: ", round(best_threshold["sensitivity"], 2)
  )
  diag_roc_curve <- ggplot(diag_roc_data, aes(x = FPR, y = TPR)) +
    geom_line(size = 1.2, color = rgb(114, 188, 213, maxColorValue = 255)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
    geom_point(data = best_point, aes(x = FPR, y = TPR), color = "black", size = 2) +
    annotate("text",
      x = best_point$FPR + 0.05, y = best_point$TPR - 0.1,
      label = best_label, size = 2, hjust = 0, color = "black"
    ) +
    annotate("text", x = 0.8, y = 0.2, label = paste("AUC =", round(diag_auc, 3)), size = 2, color = "black") +
    labs(
      title = dataset_names[[dataset]],
      x = "1 - Specificity",
      y = "Sensitivity",
      color = ""
    ) +
    my_plot_theme(legend = c(0.8, 0.2))
  qsave(diag_roc_curve, file = paste0(diag_value_dir, "/", dataset, "_roc_curve.qs"))
  ggsave(paste0(diag_value_dir, "/", dataset, "_roc_curve.pdf"),
    plot = diag_roc_curve, width = 1.625, height = 2
  )

  diag_violin_plot <- ggplot(analysis_table, aes(x = sample_type, y = risk_score, fill = sample_type)) +
    geom_violin(trim = F, color = "white") +
    geom_boxplot(width = 0.1, position = position_dodge(0.9)) +
    scale_fill_manual(values = c(
      rgb(114, 188, 213, maxColorValue = 255),
      rgb(239, 138, 071, maxColorValue = 255)
    )) +
    labs(x = NULL, y = "Risk Score", title = dataset_names[[dataset]]) +
    stat_compare_means(method = "t.test", label.x.npc = 0.3, label.y = max(analysis_table$risk_score) + 0.5, size = 2) +
    my_plot_theme(legend = "none")
  ggsave(filename = paste0(diag_value_dir, "/", dataset, "_violin.pdf"), device = "pdf",  width = 1.625, height = 2)
  qsave(diag_violin_plot, file = paste0(diag_value_dir, "/", dataset, "_violin.qs"))
}

# subgroup analysis
subgroup_dir <- file.path(result_dir, "subgroup_analysis")
stage_subgroup_dir <- file.path(subgroup_dir, "stage")
grade_subgroup_dir <- file.path(subgroup_dir, "grade")
if (!dir.exists(stage_subgroup_dir)) dir.create(stage_subgroup_dir, recursive = T)
if (!dir.exists(grade_subgroup_dir)) dir.create(grade_subgroup_dir, recursive = T)
# stage subgroup analysis
for (dataset in datasets) {
  clinical <- qread(file = paste0(clinical_dir, "/", dataset, ".qs"))
  for (stage in c("S1","S2","S3","S4")) {
    analysis_table <- subset(clinical,Stage==stage)
    result<- summary(coxph(Surv(OS,OS.Censor)~risk_score,data = analysis_table))
    result_table <- data.frame(result$concordance)%>%t()%>%data.frame()
    rownames(result_table) <- NULL
    result_table$dataset <- dataset_names[[dataset]] %>% str_wrap(width = 4)
    result_table$Stage <- stage
    result_table$sample_number <- nrow(analysis_table)
    if(stage=="S1"){combined_result <- result_table}else{combined_result <- rbind(combined_result,result_table)}
  }
  if(dataset==datasets[1]){
    combined_dataset_result <- combined_result
  }else{
    combined_dataset_result <- rbind(combined_dataset_result,combined_result)
  }
}
write.csv(combined_dataset_result,file = paste0(stage_subgroup_dir,"/stage_subgroup.csv"))
bar_col <- rev(c(rgb(231,098,084,maxColorValue = 255,alpha = 220),
                 rgb(247,170,088,maxColorValue = 255,alpha = 220),
                 rgb(170,220,224,maxColorValue = 255,alpha = 220),
                 rgb(030,070,110,maxColorValue = 255,alpha = 220)))
combined_dataset_result$Stage <- factor(combined_dataset_result$Stage,levels = c("S1","S2","S3","S4"),labels = c("I","II","III","IV"),ordered = T)
satge_subgroup_plot <- ggplot(combined_dataset_result, aes(x=dataset, y=C, fill=Stage)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.6) +
  geom_errorbar(aes(ymin=C-se.C., ymax=C+se.C.), width=.2, position=position_dodge(.6)) +
  scale_fill_manual(values = bar_col)+
  scale_y_continuous(limits = c(0,1.05),breaks = c(0,0.2,0.4,0.6,0.8,1),expand = c(0.02,0))+
  labs(title="Subgroup Analysis - Stage",y="C-Index",x=NULL,fill="Stage") +
  my_plot_theme(x.text.angle = 0,legend = "top") + 
  theme(legend.key.size = unit(0.3, "cm"),
        legend.box.margin = margin(0,0,-15,0))
qsave(satge_subgroup_plot,file = paste0(stage_subgroup_dir,"/stage_subgroup.qs"))
ggsave(filename = paste0(stage_subgroup_dir,"/stage_subgroup.pdf"),plot = satge_subgroup_plot,width = 3.25,height = 2)

# grade subgroup analysis
for (dataset in datasets) {
  clinical <- qread(file = paste0(clinical_dir, "/", dataset, ".qs"))
  for (grade_group in c("I/II","III/IV")) {
    analysis_table <- clinical[clinical$grade_group == grade_group,] %>% na.omit()
    result<- summary(coxph(Surv(OS,OS.Censor)~risk_score,data = analysis_table))
    result_table <- data.frame(result$concordance)%>%t()%>%data.frame()
    rownames(result_table) <- NULL
    result_table$dataset <- dataset_names[[dataset]] %>% str_wrap(width = 4)
    result_table$grade <- grade_group
    result_table$sample_number <- nrow(analysis_table)
    if(grade_group=="I/II"){combined_result <- result_table}else{combined_result <- rbind(combined_result,result_table)}
  }
  if(dataset==datasets[1]){
    combined_dataset_result <- combined_result
  }else{
    combined_dataset_result <- rbind(combined_dataset_result,combined_result)
  }
}
write.csv(combined_dataset_result,file = paste0(grade_subgroup_dir,"/grade_subgroup.csv"))
bar_col <- c(  rgb(030,070,110,maxColorValue = 255,alpha = 220),
               rgb(231,098,084,maxColorValue = 255,alpha = 220))
grade_subgroup_plot <- ggplot(combined_dataset_result, aes(x=dataset, y=C, fill=grade)) +
  geom_bar(stat="identity", position=position_dodge(), width=0.6) +
  geom_errorbar(aes(ymin=C-se.C., ymax=C+se.C.), width=.2, position=position_dodge(.6)) +
  scale_fill_manual(values = bar_col)+
  scale_y_continuous(limits = c(0,1.05),breaks = c(0,0.2,0.4,0.6,0.8,1),expand = c(0.02,0))+
  labs(title="Subgroup Analysis - Grade",y="C-Index",x=NULL,fill = "Grade") +
  my_plot_theme(x.text.angle = 0,legend = "top") + 
  theme(legend.key.size = unit(0.3, "cm"),
        legend.box.margin = margin(0,0,-15,0))
qsave(grade_subgroup_plot,file = paste0(grade_subgroup_dir,"/grade_subgroup.qs"))
ggsave(filename = paste0(grade_subgroup_dir,"/grade_subgroup.pdf"),plot = grade_subgroup_plot,width = 3.25,height = 2)



