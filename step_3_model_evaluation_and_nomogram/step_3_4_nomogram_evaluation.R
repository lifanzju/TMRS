library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(stringr)
library(data.table)
library(rms)
library(pec)
library(dcurves)
library(qs)
library(reshape2)

rm(list = ls())
source("scripts/utils/plot_settings.R")
exp_dir <- "data/data_for_ml/exp"
clinical_dir <- "results/step_3_model_evaluation_and_nomogram/processed_clinical_data"
result_dir <- "results/step_3_model_evaluation_and_nomogram/nomogram_evaluation"
if (!dir.exists(result_dir)) dir.create(result_dir)
datasets <- str_remove_all(list.files(exp_dir), pattern = fixed(".csv"))
dataset_names <- list(
  "TCGA" = "TCGA",
  "E_MTAB_1980" = "E-MTAB-1980",
  "CPTAC" = "CPTAC",
  "TCGA_training" = "TCGA Training",
  "TCGA_validation" = "TCGA Validation"
)
os_nomogram <- qread("results/step_3_model_evaluation_and_nomogram/nomogram_construction/os_nomogram_cox.qs")
pfs_nomogram <- qread("results/step_3_model_evaluation_and_nomogram/nomogram_construction/pfs_nomogram_cox.qs")

# compare Harrell's c-index
harrell_c_index_dir <- file.path(result_dir, "harrell_c_index")
os_harrell_c_index_dir <- file.path(harrell_c_index_dir, "os")
pfs_harrell_c_index_dir <- file.path(harrell_c_index_dir, "pfs")
if (!dir.exists(os_harrell_c_index_dir)) dir.create(os_harrell_c_index_dir, recursive = T)
if (!dir.exists(pfs_harrell_c_index_dir)) dir.create(pfs_harrell_c_index_dir, recursive = T)
for (dataset in datasets) {
  clinical <- qread(file.path(clinical_dir, paste0(dataset, ".qs")))
  colnames(clinical)[colnames(clinical) == "risk_score"] <- "TMRS"
  os_nomogram_score <- predict(os_nomogram, type = "lp", newdata = clinical)
  clinical$os_nomogram_score <- os_nomogram_score
  os_cox_list <- list(
    "Nomogram" = coxph(Surv(OS, OS.Censor) ~ os_nomogram_score, data = clinical, x = T, model = T),
    "TMRS" = coxph(Surv(OS, OS.Censor) ~ TMRS, data = clinical, x = T, model = T),
    "Age+Stage+Grade" = coxph(Surv(OS, OS.Censor) ~ Age + stage_group + grade_group, data = clinical, x = T, model = T),
    "Stage" = coxph(Surv(OS, OS.Censor) ~ stage_group, data = clinical, x = T, model = T),
    "Grade" = coxph(Surv(OS, OS.Censor) ~ grade_group, data = clinical, x = T, model = T),
    "Age" = coxph(Surv(OS, OS.Censor) ~ Age, data = clinical, x = T, model = T)
  )
  os_c_index <- c()
  se <- c()
  for (i in 1:length(os_cox_list)) {
    summary <- summary(os_cox_list[[i]])
    os_c_index <- c(os_c_index, summary$concordance[1])
    se <- c(se, summary$concordance[2])
  }
  os_c_index_table <- data.frame(c_index = os_c_index, se = se, dataset = dataset_names[[dataset]])
  rownames(os_c_index_table) <- names(os_cox_list)
  os_c_index_table <- arrange(os_c_index_table, desc(c_index))
  write.csv(os_c_index_table, file = file.path(os_harrell_c_index_dir, paste0(dataset, "_os_c_index.csv")))
  if (str_detect(dataset, "TCGA")) {
    pfs_nomogram_score <- predict(pfs_nomogram, type = "lp", newdata = clinical)
    clinical$pfs_nomogram_score <- pfs_nomogram_score
    pfs_cox_list <- list(
      "Nomogram" = coxph(Surv(PFS, PFS.Censor) ~ pfs_nomogram_score, data = clinical, x = T, model = T),
      "TMRS" = coxph(Surv(PFS, PFS.Censor) ~ TMRS, data = clinical, x = T, model = T),
      "Age+Stage+Grade" = coxph(Surv(PFS, PFS.Censor) ~ Age + stage_group + grade_group,
        data = clinical, x = T, model = T
      ),
      "Stage" = coxph(Surv(PFS, PFS.Censor) ~ stage_group, data = clinical, x = T, model = T),
      "Grade" = coxph(Surv(PFS, PFS.Censor) ~ grade_group, data = clinical, x = T, model = T),
      "Age" = coxph(Surv(PFS, PFS.Censor) ~ Age, data = clinical, x = T, model = T)
    )
    pfs_c_index <- c()
    se <- c()
    for (i in 1:length(pfs_cox_list)) {
      summary <- summary(pfs_cox_list[[i]])
      pfs_c_index <- c(pfs_c_index, summary$concordance[1])
      se <- c(se, summary$concordance[2])
    }
    pfs_c_index_table <- data.frame(c_index = pfs_c_index, se = se, dataset = dataset_names[[dataset]])
    rownames(pfs_c_index_table) <- names(pfs_cox_list)
    pfs_c_index_table <- arrange(pfs_c_index_table, desc(c_index))
    write.csv(pfs_c_index_table, file = file.path(pfs_harrell_c_index_dir, paste0(dataset, "_pfs_c_index.csv")))
  }
}

bar_col <- c(
  "Nomogram" = rgb(231, 098, 084, maxColorValue = 255, alpha = 220),
  "TMRS" = rgb(247, 170, 088, maxColorValue = 255, alpha = 220),
  "Age+Stage+Grade" = rgb(255, 230, 183, maxColorValue = 255, alpha = 220),
  "Stage" = rgb(170, 220, 224, maxColorValue = 255, alpha = 220),
  "Grade" = rgb(082, 143, 173, maxColorValue = 255, alpha = 220),
  "Age" = rgb(030, 070, 110, maxColorValue = 255, alpha = 220)
)
# os plot
os_all_harrll_c_index <- rbindlist(lapply(list.files(os_harrell_c_index_dir, full.names = T,pattern = "\\.csv$"), fread))
colnames(os_all_harrll_c_index)[1] <- "Variables"
os_all_harrll_c_index$Variables <- factor(os_all_harrll_c_index$Variables,
  levels = c("Nomogram", "TMRS", "Age+Stage+Grade", "Stage", "Grade", "Age")
)
os_all_harrll_c_index$dataset <- factor(os_all_harrll_c_index$dataset,
  levels = c("TCGA", "E_MTAB_1980", "CPTAC", "TCGA Training", "TCGA Validation")
)
os_harrll_c_index_plot <- ggplot(os_all_harrll_c_index, aes(x = dataset, y = c_index, fill = Variables)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.6) +
  geom_errorbar(aes(ymin = c_index - se, ymax = c_index + se), width = 0.3, position = position_dodge(.6)) +
  scale_fill_manual(values = bar_col) +
  scale_x_discrete(labels = str_wrap(levels(os_all_harrll_c_index$dataset),width = 4))+
  scale_y_continuous(limits = c(0, 0.9), breaks = c(0, 0.2, 0.4, 0.6, 0.8), expand = c(0.02, 0)) +
  labs(title = "Harrell's C Index (OS)", y = "C-Index", x = NULL, fill = "") +
  my_plot_theme(x.text.angle = 0,legend = "top") +
  theme(legend.key.height = unit(0.2, "cm"),
        legend.box.margin = margin(1, 0, -15, 0))
ggsave(
  file = file.path(os_harrell_c_index_dir, "os_harrell_c_index.pdf"),
  plot = os_harrll_c_index_plot, width = 3.25, height = 2
)
qsave(os_all_harrll_c_index, file = file.path(os_harrell_c_index_dir, "os_harrell_c_index.qs"))

# pfs plot
pfs_all_harrll_c_index <- rbindlist(lapply(list.files(pfs_harrell_c_index_dir, full.names = T,pattern = "\\.csv$"), fread))
colnames(pfs_all_harrll_c_index)[1] <- "Variables"
pfs_all_harrll_c_index$Variables <- factor(pfs_all_harrll_c_index$Variables,
  levels = c("Nomogram", "TMRS", "Age+Stage+Grade", "Stage", "Grade", "Age")
)
pfs_all_harrll_c_index$dataset <- factor(pfs_all_harrll_c_index$dataset,
  levels = c("TCGA", "TCGA Training", "TCGA Validation")
)
pfs_harrll_c_index_plot <- ggplot(pfs_all_harrll_c_index, aes(x = dataset, y = c_index, fill = Variables)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.6) +
  geom_errorbar(aes(ymin = c_index - se, ymax = c_index + se), width = 0.3, position = position_dodge(.6)) +
  scale_fill_manual(values = bar_col) +
  scale_x_discrete(labels = str_wrap(levels(pfs_all_harrll_c_index$dataset),width = 4))+
  scale_y_continuous(limits = c(0, 0.9), breaks = c(0, 0.2, 0.4, 0.6, 0.8), expand = c(0.02, 0)) +
  labs(title = "Harrell's C Index (PFS)", y = "C-Index", x = NULL, fill = "") +
  my_plot_theme(x.text.angle = 0,legend = "top")+
  theme(legend.key.height = unit(0.2, "cm"),
        legend.box.margin = margin(1, 0, -15, 0))
ggsave(
  file = file.path(pfs_harrell_c_index_dir, "pfs_harrell_c_index.pdf"),
  plot = pfs_harrll_c_index_plot, width = 3.25, height = 2
)
qsave(pfs_all_harrll_c_index, file = file.path(pfs_harrell_c_index_dir, "pfs_harrell_c_index.qs"))

# compare time-dependent c-index
time_dependent_c_index_dir <- file.path(result_dir, "time_dependent_c_index")
os_time_dependent_c_index_dir <- file.path(time_dependent_c_index_dir, "os")
pfs_time_dependent_c_index_dir <- file.path(time_dependent_c_index_dir, "pfs")
if (!dir.exists(os_time_dependent_c_index_dir)) dir.create(os_time_dependent_c_index_dir, recursive = T)
if (!dir.exists(pfs_time_dependent_c_index_dir)) dir.create(pfs_time_dependent_c_index_dir, recursive = T)
for (dataset in datasets) {
  clinical <- qread(file.path(clinical_dir, paste0(dataset, ".qs")))
  colnames(clinical)[colnames(clinical) == "risk_score"] <- "TMRS"
  os_nomogram_score <- predict(os_nomogram, type = "lp", newdata = clinical)
  clinical$os_nomogram_score <- os_nomogram_score
  os_cox_list <- list(
    "Nomogram" = coxph(Surv(OS, OS.Censor) ~ os_nomogram_score, data = clinical, x = T, model = T),
    "TMRS" = coxph(Surv(OS, OS.Censor) ~ TMRS, data = clinical, x = T, model = T),
    "Age+Stage+Grade" = coxph(Surv(OS, OS.Censor) ~ Age + stage_group + grade_group, data = clinical, x = T, model = T),
    "Stage" = coxph(Surv(OS, OS.Censor) ~ stage_group, data = clinical, x = T, model = T),
    "Grade" = coxph(Surv(OS, OS.Censor) ~ grade_group, data = clinical, x = T, model = T),
    "Age" = coxph(Surv(OS, OS.Censor) ~ Age, data = clinical, x = T, model = T)
  )
  if (dataset == "CPTAC") {
    eval_times <- seq(0, 2000, 100)
  } else {
    eval_times <- seq(0, 4000, 100)
  }
  os_time_c_index <- cindex(os_cox_list,
    data = clinical,
    formula = Surv(OS, OS.Censor) ~ 1,
    eval.times = eval_times
  )
  os_time_c_index_data <- os_time_c_index$AppCindex %>% data.frame()
  os_time_c_index_data$time <- os_time_c_index$time
  os_time_c_index_data <- melt(os_time_c_index_data, id.vars = "time", variable.name = "model", value.name = "c_index")
  os_time_c_index_data <- subset(os_time_c_index_data, c_index != "NaN")
  os_time_c_index_data$model <- str_replace_all(os_time_c_index_data$model,
    pattern = "Age.Stage.Grade",
    replacement = "Age+Stage+Grade"
  ) %>%
    factor(levels = c("Nomogram", "TMRS", "Age+Stage+Grade", "Stage", "Grade", "Age"))

  os_time_c_plot <- ggplot(os_time_c_index_data, aes(x = time, y = c_index, color = model)) +
    geom_line(size = 1) +
    labs(
      x = "Time (days)",
      y = "Time-dependent C-Index",
      color = "",
      title = paste0(dataset_names[[dataset]], " (OS)")
    ) +
    scale_y_continuous(limits = c(0.4, 1.15), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    scale_color_manual(values = bar_col) +
    my_plot_theme(legend = c(0.6, 0.9)) +
    theme(
      legend.key.height = unit(0.1, "cm"),
    )

  ggsave(file.path(os_time_dependent_c_index_dir, paste0(dataset, "_os_time_dependent_c_index.pdf")),
    plot = os_time_c_plot, width = 1.625, height = 2
  )
  qsave(os_time_c_plot, file = file.path(os_time_dependent_c_index_dir, paste0(dataset, "_os_time_dependent_c_index.qs")))
  if (str_detect(dataset, pattern = "TCGA")) {
    pfs_nomogram_score <- predict(pfs_nomogram, type = "lp", newdata = clinical)
    clinical$pfs_nomogram_score <- pfs_nomogram_score
    clinical <- subset(clinical, !is.na(PFS) & !is.na(PFS.Censor))
    pfs_cox_list <- list(
      "Nomogram" = coxph(Surv(PFS, PFS.Censor) ~ pfs_nomogram_score, data = clinical, x = T, model = T),
      "TMRS" = coxph(Surv(PFS, PFS.Censor) ~ TMRS, data = clinical, x = T, model = T),
      "Age+Stage+Grade" = coxph(Surv(PFS, PFS.Censor) ~ Age + stage_group + grade_group, data = clinical, x = T, model = T),
      "Stage" = coxph(Surv(PFS, PFS.Censor) ~ stage_group, data = clinical, x = T, model = T),
      "Grade" = coxph(Surv(PFS, PFS.Censor) ~ grade_group, data = clinical, x = T, model = T),
      "Age" = coxph(Surv(PFS, PFS.Censor) ~ Age, data = clinical, x = T, model = T)
    )
    pfs_time_c_index <- cindex(pfs_cox_list,
      data = clinical,
      formula = Surv(PFS, PFS.Censor) ~ 1,
      eval.times = eval_times
    )
    pfs_time_c_index_data <- pfs_time_c_index$AppCindex %>% data.frame()
    pfs_time_c_index_data$time <- pfs_time_c_index$time
    pfs_time_c_index_data <- melt(pfs_time_c_index_data, id.vars = "time", variable.name = "model", value.name = "c_index")
    pfs_time_c_index_data <- subset(pfs_time_c_index_data, c_index != "NaN")
    pfs_time_c_index_data$model <- str_replace_all(pfs_time_c_index_data$model,
      pattern = "Age.Stage.Grade",
      replacement = "Age+Stage+Grade"
    ) %>%
      factor(levels = c("Nomogram", "TMRS", "Age+Stage+Grade", "Stage", "Grade", "Age"))
    pfs_time_c_plot <- ggplot(pfs_time_c_index_data, aes(x = time, y = c_index, color = model)) +
      geom_line(size = 1) +
      labs(
        x = "Time (days)",
        y = "Time-dependent Concordance Index",
        color = "",
        title = paste0(dataset_names[[dataset]], " (PFS)")
      ) +
      scale_y_continuous(limits = c(0.4, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
      scale_color_manual(values = bar_col) +
      my_plot_theme(legend = c(0.6, 0.9)) +
      theme(
        legend.key.height = unit(0.1, "cm"),
      )
    ggsave(file.path(pfs_time_dependent_c_index_dir, paste0(dataset, "_pfs_time_dependent_c_index.pdf")),
      plot = pfs_time_c_plot, width = 1.625, height = 2
    )
    qsave(pfs_time_c_plot, file = file.path(pfs_time_dependent_c_index_dir, paste0(dataset, "_pfs_time_dependent_c_index.qs")))
  }
}


# calibration plot
calibration_dir <- file.path(result_dir, "calibration")
os_calibration_dir <- file.path(calibration_dir, "os")
pfs_calibration_dir <- file.path(calibration_dir, "pfs")
if (!dir.exists(os_calibration_dir)) dir.create(os_calibration_dir, recursive = T)
if (!dir.exists(pfs_calibration_dir)) dir.create(pfs_calibration_dir, recursive = T)
for (dataset in datasets) {
  clinical <- qread(file.path(clinical_dir, paste0(dataset, ".qs")))
  colnames(clinical)[colnames(clinical) == "risk_score"] <- "TMRS"
  os_nomogram_score <- predict(os_nomogram, type = "lp", newdata = clinical)
  clinical$os_nomogram_score <- os_nomogram_score
  ddDD <- datadist(clinical)
  options(datadist = "ddDD")
  os_calibration_results <- data.frame()
  for (year in c(1, 3, 5)) {
    os_cox <- cph(Surv(OS, OS.Censor) ~ os_nomogram_score,
      data = clinical, x = TRUE, y = TRUE,
      surv = TRUE, time.inc = 365 * year
    )
    cal_os <- rms::calibrate(os_cox,
      cmethod = "KM",
      method = "boot",
      u = 365 * year,
      m = round(nrow(clinical) / 3.5),
      B = 1000
    )
    cal_os_data <- data.frame(
      predicted_values = cal_os[, "mean.predicted"],
      actual_values = cal_os[, "KM"],
      std_err = cal_os[, "std.err"],
      year = paste0(year, ifelse(year == 1, " year", " years"))
    )
    os_calibration_results <- rbind(os_calibration_results, cal_os_data)
  }
  os_calibration_results$year <- factor(os_calibration_results$year)
  os_calibration_plot <- ggplot(os_calibration_results, aes(x = predicted_values, y = actual_values, color = year)) +
    geom_point() +
    geom_line() +
    geom_errorbar(
      aes(
        ymin = actual_values - std_err,
        ymax = actual_values + std_err
      ),
      width = 0.02
    ) +
    scale_color_manual(values = c(
      "1 year" = rgb(170, 220, 224, maxColorValue = 255, alpha = 220),
      "3 years" = rgb(247, 170, 088, maxColorValue = 255, alpha = 220),
      "5 years" = rgb(231, 098, 084, maxColorValue = 255, alpha = 220)
    )) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    labs(
      x = "Predicted Survival",
      y = "Observed Survival",
      title = paste0(dataset_names[[dataset]], " (OS)"),
      color = ""
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    my_plot_theme(legend = c(0.75, 0.25))+
    theme(legend.key.height = unit(0.3, "cm"))
  ggsave(file.path(os_calibration_dir, paste0(dataset, "_os_calibration.pdf")),
    width = 1.625,
    height = 2, plot = os_calibration_plot
  )
  qsave(os_calibration_plot, file = file.path(os_calibration_dir, paste0(dataset, "_os_calibration.qs")))
  if (str_detect(dataset, "TCGA")) {
    pfs_nomogram_score <- predict(pfs_nomogram, type = "lp", newdata = clinical)
    clinical$pfs_nomogram_score <- pfs_nomogram_score
    clinical <- subset(clinical, !is.na(PFS) & !is.na(PFS.Censor))
    ddDD <- datadist(clinical)
    options(datadist = "ddDD")
    pfs_calibration_results <- data.frame()
    for (year in c(1, 3, 5)) {
      pfs_cox <- cph(Surv(PFS, PFS.Censor) ~ pfs_nomogram_score,
        data = clinical, x = TRUE, y = TRUE,
        surv = TRUE, time.inc = 365 * year
      )
      cal_pfs <- rms::calibrate(pfs_cox,
        cmethod = "KM",
        method = "boot",
        u = 365 * year,
        m = round(nrow(clinical) / 3.5),
        B = 1000
      )
      cal_pfs_data <- data.frame(
        predicted_values = cal_pfs[, "mean.predicted"],
        actual_values = cal_pfs[, "KM"],
        std_err = cal_pfs[, "std.err"],
        year = paste0(year, ifelse(year == 1, " year", " years"))
      )
      pfs_calibration_results <- rbind(pfs_calibration_results, cal_pfs_data)
    }
    pfs_calibration_results$year <- factor(pfs_calibration_results$year)
    pfs_calibration_plot <- ggplot(pfs_calibration_results, aes(x = predicted_values, y = actual_values, color = year)) +
      geom_point() +
      geom_line() +
      geom_errorbar(
        aes(
          ymin = actual_values - std_err,
          ymax = actual_values + std_err
        ),
        width = 0.02
      ) +
      scale_color_manual(values = c(
        "1 year" = rgb(170, 220, 224, maxColorValue = 255, alpha = 220),
        "3 years" = rgb(247, 170, 088, maxColorValue = 255, alpha = 220),
        "5 years" = rgb(231, 098, 084, maxColorValue = 255, alpha = 220)
      )) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
      labs(
        x = "Predicted Survival",
        y = "Observed Survival",
        title = paste0(dataset_names[[dataset]], " (PFS)"),
        color = ""
      ) +
      scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
      scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
      my_plot_theme(legend = c(0.75, 0.25))+
      theme(legend.key.height = unit(0.3, "cm"))
    ggsave(file.path(pfs_calibration_dir, paste0(dataset, "_pfs_calibration.pdf")),
      width = 1.625,
      height = 2, plot = pfs_calibration_plot
    )
    qsave(pfs_calibration_plot, file = file.path(pfs_calibration_dir, paste0(dataset, "_pfs_calibration.qs")))
  }
}

# decision curve analysis
dca_dir <- file.path(result_dir, "dca")
os_dca_dir <- file.path(dca_dir, "os")
pfs_dca_dir <- file.path(dca_dir, "pfs")
if (!dir.exists(os_dca_dir)) dir.create(os_dca_dir, recursive = T)
if (!dir.exists(pfs_dca_dir)) dir.create(pfs_dca_dir, recursive = T)
for (dataset in datasets) {
  clinical <- qread(file.path(clinical_dir, paste0(dataset, ".qs")))
  clinical <- na.omit(clinical)
  colnames(clinical)[colnames(clinical) == "risk_score"] <- "TMRS"
  os_nomogram_score <- predict(os_nomogram, type = "lp", newdata = clinical)
  clinical$os_nomogram_score <- os_nomogram_score
  os_cox_list <- list(
    "Nomogram" = coxph(Surv(OS, OS.Censor) ~ os_nomogram_score, data = clinical, x = T, model = T),
    "TMRS" = coxph(Surv(OS, OS.Censor) ~ TMRS, data = clinical, x = T, model = T),
    "Age+Stage+Grade" = coxph(Surv(OS, OS.Censor) ~ Age + stage_group + grade_group, data = clinical, x = T, model = T),
    "Stage" = coxph(Surv(OS, OS.Censor) ~ stage_group, data = clinical, x = T, model = T),
    "Grade" = coxph(Surv(OS, OS.Censor) ~ grade_group, data = clinical, x = T, model = T),
    "Age" = coxph(Surv(OS, OS.Censor) ~ Age, data = clinical, x = T, model = T)
  )
  # calculate probability of survival at 1, 3, 5 years
  for (year in c(1, 3, 5)) {
    formula_terms <- c()
    for (i in 1:length(os_cox_list)) {
      model_name <- names(os_cox_list)[i]
      os_prob_column <- paste0(model_name, "_prob_", year, "_year") %>% str_replace_all(pattern = fixed("+"), replacement = "_")
      clinical[, os_prob_column] <- c(1 - predictSurvProb(
        os_cox_list[[i]],
        newdata = clinical,
        times = year * 365
      )[, 1])
      formula_terms <- c(formula_terms, paste0(os_prob_column))
    }
    formula_str <- paste("Surv(OS, OS.Censor) ~", paste(formula_terms, collapse = " + "))
    formula_obj <- as.formula(formula_str)
    dca_result <- dca(formula_obj,
      data = clinical,
      label = setNames(as.list(names(os_cox_list)), formula_terms),
      time = year * 365
    )
    dca_plot <- plot(dca_result, smooth = T) +
      scale_color_manual(values = c(bar_col, "Treat All" = "gray", "Treat None" = "black")) +
      ggtitle(str_wrap(
        paste0(dataset_names[[dataset]], 
               " (OS - ", year, " ", 
               ifelse(year == 1, "year)", "years)")
               ),
        width = 20)) +
      my_plot_theme(legend = c(0.65, 0.85)) +
      theme(legend.key.height = unit(0.15, "cm"),
            legend.key.width = unit(0.2, "cm"))
    ggsave(file.path(os_dca_dir, paste0(dataset, "_os_dca_", year, "_year.pdf")),
      width = 1.625, height = 2, plot = dca_plot
    )
    qsave(dca_plot, file = file.path(os_dca_dir, paste0(dataset, "_os_dca_", year, "_year.qs")))
  }
  if (str_detect(dataset, "TCGA")) {
    pfs_nomogram_score <- predict(pfs_nomogram, type = "lp", newdata = clinical)
    clinical$pfs_nomogram_score <- pfs_nomogram_score
    pfs_cox_list <- list(
      "Nomogram" = coxph(Surv(PFS, PFS.Censor) ~ pfs_nomogram_score, data = clinical, x = T, model = T),
      "TMRS" = coxph(Surv(PFS, PFS.Censor) ~ TMRS, data = clinical, x = T, model = T),
      "Age+Stage+Grade" = coxph(Surv(PFS, PFS.Censor) ~ Age + stage_group + grade_group, data = clinical, x = T, model = T),
      "Stage" = coxph(Surv(PFS, PFS.Censor) ~ stage_group, data = clinical, x = T, model = T),
      "Grade" = coxph(Surv(PFS, PFS.Censor) ~ grade_group, data = clinical, x = T, model = T),
      "Age" = coxph(Surv(PFS, PFS.Censor) ~ Age, data = clinical, x = T, model = T)
    )
    # calculate probability of survival at 1, 3, 5 years
    for (year in c(1, 3, 5)) {
      formula_terms <- c()
      for (i in 1:length(pfs_cox_list)) {
        model_name <- names(pfs_cox_list)[i]
        pfs_prob_column <- paste0(model_name, "_prob_", year, "_year") %>% str_replace_all(pattern = fixed("+"), replacement = "_")
        clinical[, pfs_prob_column] <- c(1 - predictSurvProb(
          pfs_cox_list[[i]],
          newdata = clinical,
          times = year * 365
        )[, 1])
        formula_terms <- c(formula_terms, paste0(pfs_prob_column))
      }
      formula_str <- paste("Surv(PFS, PFS.Censor) ~", paste(formula_terms, collapse = " + "))
      formula_obj <- as.formula(formula_str)
      dca_result <- dca(formula_obj,
        data = clinical,
        label = setNames(as.list(names(pfs_cox_list)), formula_terms),
        time = year * 365
      )
      dca_plot <- plot(dca_result, smooth = T) +
        scale_color_manual(values = c(bar_col, "Treat All" = "gray", "Treat None" = "black")) +
        ggtitle(str_wrap(
          paste0(dataset_names[[dataset]], 
                 " (PFS - ", year, " ", 
                 ifelse(year == 1, "year", "years)")
                 ),
          width = 20)) +
        my_plot_theme(legend = c(0.65, 0.85)) +
        theme(legend.key.height = unit(0.15, "cm"),
              legend.key.width = unit(0.2, "cm"))
      ggsave(file.path(pfs_dca_dir, paste0(dataset, "_pfs_dca_", year, "_year.pdf")),
        width = 1.625, height = 2, plot = dca_plot
      )
      qsave(dca_plot, file = file.path(pfs_dca_dir, paste0(dataset, "_pfs_dca_", year, "_year.qs")))
    }
  }
}

