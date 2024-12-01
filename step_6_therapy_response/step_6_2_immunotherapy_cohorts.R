library(dplyr)
library(stringr)
library(ggplot2)
library(gbm)
library(survival)
library(survminer)
library(data.table)
library(qs)

rm(list=ls())
source("scripts/utils/plot_settings.R")

input_dir <- "data/immunotherapy_cohort"
output_dir <- "results/step_6_therapy_response/immunotherapy_cohorts"
tmrs_model <- qread("results/step_2_model_construction/machine_learning/final_model.qs")
datasets <- list.files(input_dir)

for (dataset in datasets) {
  exp <- fread(file.path(input_dir, dataset, "exp.csv")) %>% 
    data.frame(row.names = 1) %>% 
    t() %>% 
    scale() %>% 
    data.frame()
  clinical <- fread(file.path(input_dir, dataset, "clinical.csv")) %>%
    data.frame(row.names = 1)
  clinical <- clinical[rownames(exp),]
  tmrs_score <- predict(tmrs_model, newdata = exp, type = "link")
  clinical$tmrs_score <- tmrs_score
  clinical$dataset <- dataset
  if (dataset == "GSE91061") {
    clinical$OS <- clinical$OS/30
  }
  if (dataset != "JAVELIN_Renal_101") {
    cutoff <- surv_cutpoint(clinical, time = "OS", event = "censor",
                            variables = "tmrs_score")$cutpoint$cutpoint
    clinical$tmrs_score_group <- ifelse(clinical$tmrs_score > cutoff ,"High","Low")
    K_M_model <- survfit(Surv(OS,censor)~tmrs_score_group,data = clinical)
    K_M_Plot <- ggsurvplot(K_M_model, 
                           data = clinical,
                           title = dataset,
                           pval = TRUE,
                           palette = c("firebrick2","dodgerblue3"), 
                           legend.title = "",
                           legend.labs = c(
                             paste0("High"),
                             paste0("Low")
                           ),
                           risk.table = FALSE,
                           cumevents = FALSE,  
                           tables.y.text = FALSE,
                           conf.int=F,
                           pval.size = 2.5,
                           censor.size = 2.5,
                           xlab = "Survival (Months)", 
                           ylab="OS") 
    K_M_Plot$plot<- K_M_Plot$plot+ 
      my_plot_theme(legend = "top")+
      theme(legend.box.margin = margin(-5,0,-15,0))
    ggsave(file = file.path(output_dir, paste0(dataset, "_KM_plot.pdf")), K_M_Plot$plot, width = 1.625, height = 2)
    if (dataset == "checkmate025") {
      pfs_cutoff <- surv_cutpoint(clinical, time = "PFS", event = "PFS_CNSR",
                                  variables = "tmrs_score")$cutpoint$cutpoint
      clinical$tmrs_score_pfs_group <- ifelse(clinical$tmrs_score > pfs_cutoff ,"High","Low")
      K_M_model_pfs <- survfit(Surv(PFS,PFS_CNSR)~tmrs_score_pfs_group,data = clinical)
      K_M_Plot_pfs <- ggsurvplot(K_M_model_pfs, 
                                 data = clinical,
                                 title = dataset,
                                 pval = TRUE,
                                 palette = c("firebrick2","dodgerblue3"), 
                                 legend.title = "",
                                 legend.labs = c(
                                   paste0("High"),
                                   paste0("Low")
                                 ),
                                 risk.table = FALSE,
                                 cumevents = FALSE,  
                                 tables.y.text = FALSE,
                                 conf.int=F,
                                 pval.size = 2.5,
                                 censor.size = 2.5,
                                 xlab = "Survival (Months)", 
                                 ylab="PFS")
      K_M_Plot_pfs$plot<- K_M_Plot_pfs$plot+ 
        my_plot_theme(legend = "top")+
        theme(legend.box.margin = margin(-5,0,-15,0))
      ggsave(file = file.path(output_dir, paste0(dataset, "_KM_plot_pfs.pdf")), K_M_Plot_pfs$plot, width = 1.625, height = 2)
    }
  }else{
    cutoff <- surv_cutpoint(clinical, time = "PFS_P", event = "censor",
                            variables = "tmrs_score")$cutpoint$cutpoint
    clinical$tmrs_score_group <- ifelse(clinical$tmrs_score > cutoff ,"High","Low")
    K_M_model <- survfit(Surv(PFS_P,censor)~tmrs_score_group,data = clinical)
    K_M_Plot <- ggsurvplot(K_M_model, 
                           data = clinical,
                           title = dataset,
                           pval = TRUE,
                           palette = c("firebrick2","dodgerblue3"), 
                           legend.title = "",
                           legend.labs = c(
                             paste0("High"),
                             paste0("Low")
                           ),
                           risk.table = FALSE,
                           cumevents = FALSE,  
                           tables.y.text = FALSE,
                           conf.int=F,
                           pval.size = 2.5,
                           censor.size = 2.5,
                           xlab = "Survival (Months)", 
                           ylab="PFS") 
    K_M_Plot$plot<- K_M_Plot$plot+ 
      my_plot_theme(legend = "top")+
      theme(legend.box.margin = margin(-5,0,-15,0))
    ggsave(file = file.path(output_dir, paste0(dataset, "_KM_plot.pdf")), K_M_Plot$plot, width = 1.625, height = 2)
  }
}
