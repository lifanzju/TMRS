library(survival)
library(dplyr)
library(stringr)
library(parallel)
library(tibble)
library(BART)
library(qs)
library(combinat)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(ggplotify)
rm(list = ls())

# data preprocessing
unicox_intersection <- read.csv("results/step_2_model_construction/uni_cox/sig_uni_cox_intersection.csv")[, 1]
exp_dir <- "data/Datasets/proprocessed_exp"
clinical_dir <- "data/Datasets/Clinical"
results_dir <- "results/step_2_model_construction/machine_learning"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = T)
processed_data_for_ml_dir <- "data/data_for_ml"
if (!dir.exists(processed_data_for_ml_dir)) dir.create(processed_data_for_ml_dir, recursive = T)

datasets <- str_remove_all(list.files(exp_dir), pattern = fixed(".csv"))
all_datasets <- list()
for (dataset in datasets) {
  all_datasets[[dataset]] <- list(
    EXP = t(read.csv(file = paste0(exp_dir, "/", dataset, ".csv"), header = T, row.names = 1)),
    Clinical = read.delim(file = paste0(clinical_dir, "/", dataset, "_Clinical.txt"), header = T, row.names = NULL)
  )
}
preparation_step_1 <- function(x) {
  NA_ID <- x$Clinical$CaseID[is.na(x$Clinical$OS) | is.na(x$Clinical$OS.Censor)]
  Negative_ID <- x$Clinical$CaseID[x$Clinical$OS <= 0]
  Bad_ID <- c(NA_ID, Negative_ID)
  x$Clinical <- x$Clinical[!x$Clinical$CaseID %in% Bad_ID, ]
  x$EXP <- x$EXP[!row.names(x$EXP) %in% Bad_ID, ]
  x$EXP <- scale(x$EXP)
  return(x)
}
preparation_step_2 <- function(x, genelist) {
  x$EXP <- x$EXP[match(x$Clinical$CaseID, row.names(x$EXP)), ]
  x$EXP <- x$EXP[, genelist]
  x$Survival <- data.matrix(Surv(time = as.double(x$Clinical$OS), event = as.double(x$Clinical$OS.Censor)))
  x$Analysis_table <- as.data.frame(cbind(x$Survival, x$EXP))
  rownames(x$Survival) <- x$Clinical$CaseID
  return(x)
}
all_datasets <- mclapply(all_datasets, preparation_step_1, mc.cores = 3)

ml_exp_dir <- paste0(processed_data_for_ml_dir, "/exp")
ml_clinical_dir <- paste0(processed_data_for_ml_dir, "/clinical")
if (!dir.exists(ml_exp_dir)) dir.create(ml_exp_dir)
if (!dir.exists(ml_clinical_dir)) dir.create(ml_clinical_dir)

for (dataset in datasets) {
  exp <- all_datasets[[dataset]]$EXP
  clinical <- all_datasets[[dataset]]$Clinical
  exp <- t(exp[match(clinical$CaseID, row.names(exp)), ])
  write.csv(exp, file = paste0(ml_exp_dir, "/", dataset, ".csv"))
  write.csv(clinical, file = paste0(ml_clinical_dir, "/", dataset, ".csv"),row.names = F)
}
# split TCGA dataset to training dataset and internal validation dataset
set.seed(527)
training_samples <- sample(rownames(all_datasets$TCGA$EXP), size = round(0.7 * nrow(all_datasets$TCGA$EXP)))
TCGA_training_exp <- all_datasets$TCGA$EXP[rownames(all_datasets$TCGA$EXP) %in% training_samples, ] 
TCGA_training_clinical <- subset(all_datasets$TCGA$Clinical, CaseID %in% training_samples)
TCGA_validation_exp <- all_datasets$TCGA$EXP[!rownames(all_datasets$TCGA$EXP) %in% training_samples, ] 
TCGA_validation_clinical <- subset(all_datasets$TCGA$Clinical, CaseID %in% row.names(TCGA_validation_exp))
write.csv(t(TCGA_training_exp), file = paste0(ml_exp_dir, "/TCGA_training.csv"))
write.csv(TCGA_training_clinical, file = paste0(ml_clinical_dir, "/TCGA_training.csv"),row.names = F)
write.csv(t(TCGA_validation_exp), file = paste0(ml_exp_dir, "/TCGA_validation.csv"))
write.csv(TCGA_validation_clinical, file = paste0(ml_clinical_dir, "/TCGA_validation.csv"),row.names = F)
all_datasets[["TCGA_validation"]] <- list(EXP = TCGA_validation_exp, Clinical = TCGA_validation_clinical)
all_datasets[["TCGA_training"]] <- list(EXP = TCGA_training_exp, Clinical = TCGA_training_clinical)
all_datasets$TCGA <- NULL
all_datasets <- mclapply(all_datasets, preparation_step_2, unicox_intersection, mc.cores = 4)
qsave(all_datasets, file = paste0(processed_data_for_ml_dir, "/all_datasets.qs"))

# machine learning
# source the algorithm functions
source("scripts/utils/ml_algorithm.R")

# combine algorithms
screenable <- c("RSF", "Coxboost", "superPC", "Lasso", "Elastic_Net", "Stepwise_cox")
unscreenable <- c("SurvivalSVM", "GBM", "plsRcox", "Ridge")
All_methods <- c(screenable, unscreenable)
combinations_of_screenable <- combn(screenable, 2)
Permutations_of_scrrenable <- apply(combinations_of_screenable, 2, FUN = permn)
Permutations_list <- list()
for (i in c(1:length(All_methods))) {
  Permutations_list[[i]] <- All_methods[i]
}
for (i in c(1:length(Permutations_of_scrrenable))) {
  Permutations_list[[length(All_methods) + 2 * i - 1]] <- Permutations_of_scrrenable[[i]][[1]]
  Permutations_list[[length(All_methods) + 2 * i]] <- Permutations_of_scrrenable[[i]][[2]]
}
Permutations_list <- list()
combinations <- expand.grid(i = 1:length(screenable), a = 1:length(unscreenable))
for (row in seq_len(nrow(combinations))) {
  i <- combinations$i[row]
  a <- combinations$a[row]
  Permutations_list[[row]] <- c(screenable[i], unscreenable[a])
}
# run ml 
training_dataset <- all_datasets$TCGA_training
seed <- 527
args <- list(TrainingDataset=training_dataset,ValidationDatasets=all_datasets,seed=seed)

cores <- makeCluster(6)
All_functions <- paste0(All_methods,"_function")
clusterExport(cores, All_functions)#######将定义好的函数读入到每一个核心中去
clusterEvalQ(cores, {
  for (pkg in c("glmnet","survival","dplyr","stringr","randomForestSRC","plsRcox","superpc","gbm","CoxBoost","survivalsvm","tibble","BART","combinat")) {
    library(pkg, character.only = TRUE)
  }
})
All_Results <- parLapply(cl=cores,Permutations_list,fun  = ml_run,args=args)
ML_names <- sapply(Permutations_list, function(x){if (length(x)==1) {name <- x}else{name <- paste0(x[1]," + ",x[2])}})
names(All_Results) <- ML_names
qsave(All_Results,file = paste0(results_dir,"/all_ml_results.qs"))
stopCluster(cores)
All_ML_C_Indexs <- data.frame()
for (ML in names(All_Results)) {
  if(!(str_detect(ML,pattern = fixed("Stepwise_cox"))|str_detect(ML,pattern = fixed("Elastic_Net")))){
    C_index <- t(All_Results[[ML]]$C_Index)
    rownames(C_index) <- ML
    All_ML_C_Indexs <- rbind(All_ML_C_Indexs,C_index)
  }else{
    functions_used <- str_split(ML,pattern = fixed(" + "),simplify = T)[1,]
    if(length(functions_used)==1){
      ML_result <- All_Results[[ML]]
      C_indexs_table <- t(sapply(ML_result, function(x){x$C_Index}))
      original_names <- rownames(C_indexs_table)
      newrownames <- paste0(ML," [",original_names,"]")
      rownames(C_indexs_table) <- newrownames
      All_ML_C_Indexs <- rbind(All_ML_C_Indexs,C_indexs_table)
    }else{
      if(functions_used[1]%in%c("Elastic_Net","Stepwise_cox")&functions_used[2]%in%c("Elastic_Net","Stepwise_cox")){
        ML_result <- All_Results[[ML]]
        for (index in names(ML_result)) {
          ML_result_index <- ML_result[[index]]
          C_indexs_table <- t(sapply(ML_result_index, function(x){x$C_Index}))
          original_names <- rownames(C_indexs_table)
          newrownames <- paste0(functions_used[1]," [",index,"] ","+ ",functions_used[2]," [",original_names,"]")
          rownames(C_indexs_table) <- newrownames
          All_ML_C_Indexs <- rbind(All_ML_C_Indexs,C_indexs_table)
        }
      }else{
        ML_result <- All_Results[[ML]]
        C_indexs_table <- t(sapply(ML_result, function(x){x$C_Index}))
        original_names <- rownames(C_indexs_table)
        if(functions_used[1]%in%c("Elastic_Net","Stepwise_cox")){
          newrownames <- paste0(functions_used[1]," [",original_names,"]"," + ",functions_used[2])
          rownames(C_indexs_table) <- newrownames
          All_ML_C_Indexs <- rbind(All_ML_C_Indexs,C_indexs_table)
        }else{
          newrownames <- paste0(ML," [",original_names,"]")
          rownames(C_indexs_table) <- newrownames
          All_ML_C_Indexs <- rbind(All_ML_C_Indexs,C_indexs_table)
        }
      }
    }
  }
}
All_ML_C_Indexs$All_Mean <- apply(All_ML_C_Indexs,MARGIN = 1,mean)
All_ML_C_Indexs$Validation_Mean <- apply(All_ML_C_Indexs[,paste0(c("CPTAC","E_MTAB_1980","TCGA_validation"),".concordance")],MARGIN = 1,mean)
All_ML_C_Indexs <- arrange(All_ML_C_Indexs,desc(All_Mean))
write.csv(All_ML_C_Indexs,file = paste0(results_dir,"/All_ML_C_Indexs.csv"))
# get the info of selected model
Screen_Model <- step(coxph(Surv(time,status)~.,data = training_dataset$Analysis_table),direction = "backward")
qsave(Screen_Model,file = paste0(results_dir,"/Screen_Model.qs"))
Screen_result <- Stepwise_cox_function(TrainingDataset = training_dataset,ValidationDatasets = all_datasets,seed = seed)
Screen_result <- Screen_result$backward
Variables_for_model_building <- Screen_result$Selected_Variable
Training_dataset_for_model_building <- list()
Training_dataset_for_model_building[["EXP"]] <- training_dataset$EXP[,Variables_for_model_building]
Training_dataset_for_model_building[["Analysis_table"]] <- cbind(training_dataset$Analysis_table[,1:2],training_dataset$Analysis_table[,Variables_for_model_building])
set.seed(seed)
fit0 <- gbm(formula = Surv(time,status)~.,data = Training_dataset_for_model_building$Analysis_table,distribution = 'coxph',
            n.trees = 10000,
            interaction.depth = 3,
            n.minobsinnode = 10,
            shrinkage = 0.001,
            cv.folds = 10,n.cores = 6)
best <- which.min(fit0$cv.error)
set.seed(seed)
final_model <- gbm(formula = Surv(time,status)~.,data = Training_dataset_for_model_building$Analysis_table,distribution = 'coxph',
                   n.trees = best,
                   interaction.depth = 3,
                   n.minobsinnode = 10,
                   shrinkage = 0.001,
                   cv.folds = 10,n.cores = 6)
Risk_Score_List <- lapply(all_datasets, function(x){
  risk_score <- as.numeric(predict(final_model,newdata=x$Analysis_table,type = "link",n.trees = best))
  names(risk_score) <- row.names(x$Analysis_table)
  rs <- cbind(x$Survival,risk_score)
  rs <- as.data.frame(rs)
  return(rs)
})
C_indexs <- sapply(Risk_Score_List, function(x){
  cox_regression <- coxph(Surv(time,status)~risk_score,data = x)
  c_index <- cox_regression$concordance[6]
  return(c_index)
})
for (dataset_name in names(Risk_Score_List)) {
  write.csv(Risk_Score_List[[dataset_name]],file = paste0(results_dir,"/",dataset_name,"_risk_scores.csv"))
}
qsave(final_model,file = paste0(results_dir,"/final_model.qs"))

# relative importance
source("scripts/utils/plot_settings.R")
Relative_influence <- summary(final_model)
Relative_influence$var <- factor(Relative_influence$var,levels = rev(Relative_influence$var),ordered = T)
write.csv(Relative_influence$var,file = paste0(results_dir,"/variables_of_final_model.csv"),col.names = F,row.names = F)
col.low <- rgb(253,217,133,maxColorValue=255,alpha = 220)
col.high <- rgb(231,098,084,maxColorValue = 255,alpha = 220)
relative_influence_plot <- ggplot(Relative_influence,aes(x=rel.inf, y= var, fill=rel.inf)) + geom_col(width = 0.7)+ 
  scale_fill_gradient(low=col.low,high=col.high)+theme_bw()+
  xlab("Relative Importance")+ylab("")+
  scale_x_continuous(breaks = c(0,5,10,15),limits = c(0,15))+
  my_plot_theme(legend = "none")
qsave(relative_influence_plot,file = paste0(results_dir,"/Relative_influence.qs"))
ggsave(filename = paste0(results_dir,"/relative_importance_of_variables.pdf"),
       plot = relative_influence_plot,width = 2,height = 2)

# Partial Dependence Plots
variable_names <- summary(final_model)$var
Partial_Dependence <- data.frame()
for (i in variable_names) {
  PP <- plot(final_model,i.var = i,return.grid = T,continuous.resolution = 100,type = "link")
  colnames(PP)[1] <- "x"
  PP$Gene <- i
  Partial_Dependence<- rbind(Partial_Dependence,PP)
}

start_label_data <- Partial_Dependence %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarize(
    x_start = dplyr::first(x),
    y_start = dplyr::first(y),
    .groups = "drop"
  )

partial_dependence_plot <- ggplot(Partial_Dependence,aes(x=x,y=y,color=Gene))+
  geom_line()+theme_bw()+labs(x="Gene Expression",y="Hazard Ratio (log)",title = "Partial Dependence Plots")+
  scale_color_manual(values = c(pal_npg(alpha = 0.8)(7),pal_jama(alpha = 0.8)(7)))+
  geom_text(data = start_label_data, 
            aes(x = x_start, y = y_start, label = Gene),check_overlap = F, size = 2,
            hjust = 0.1, vjust = -0.3)+
  my_plot_theme()+
  theme(legend.position = "none")
qsave(partial_dependence_plot,file = paste0(results_dir,"/partial_dependence_plot.qs"))
ggsave(filename = paste0(results_dir,"/Partial_Dependence_Plots.pdf"),width = 2,height = 2)

# c-index heatmap

Selected_ML_C_Indexs <- All_ML_C_Indexs[c(1:50),]
Selected_C_indexs_main_table <- Selected_ML_C_Indexs[,c(4,3,2,1)]
dataset_name<- str_remove_all(colnames(Selected_C_indexs_main_table),pattern = fixed(".concordance"))
colnames(Selected_C_indexs_main_table) <- dataset_name
shorten_row_names <- function(df) {
  row_names <- rownames(df)
  replacements <- list(
    "Stepwise_cox" = "SC",
    "Elastic_Net" = "EN",
    "Lasso" = "Lasso",
    "Coxboost" = "CB",
    "superPC" = "sPC",
    "SurvivalSVM" = "SSVM",
    "RSF" = "RSF",
    "\\[backward\\]" = "(bwd)",
    "\\[forward\\]" = "(fwd)",
    "\\[both\\]" = "(both)",
    "alpha=" = "α="
  )
  for (pattern in names(replacements)) {
    replacement <- replacements[[pattern]]
    row_names <- gsub(pattern, replacement, row_names)
  }
  rownames(df) <- gsub(" ", "", row_names)
  return(df)
}
Selected_C_indexs_main_table <- shorten_row_names(Selected_C_indexs_main_table)




Selected_C_indexs_main_table <- as.matrix(Selected_C_indexs_main_table)
main_table_color <- colorRamp2(c(0.63,0.7,0.8,0.9),c("dodgerblue2","white","gold","firebrick2"))
# right annotation
Mean_of_All_Datasets <- Selected_ML_C_Indexs$All_Mean
names(Mean_of_All_Datasets) <- rownames(Selected_ML_C_Indexs)
Mean_of_Validation_Datasets <- Selected_ML_C_Indexs$Validation_Mean
names(Mean_of_Validation_Datasets) <- rownames(Selected_ML_C_Indexs)
Mean_of_All_Datasets_scaled <- (Mean_of_All_Datasets-0.7)/(0.78-0.7)
nrow <- length(Mean_of_All_Datasets_scaled)
Right <- rowAnnotation(
  gap_anno = anno_empty(width = unit(1.5, "mm"), border = FALSE),
  Mean_of_Validation_Datasets = anno_barplot(
    Mean_of_Validation_Datasets,
    baseline = 0.7,
    ylim = c(0.7, 0.78),
    axis_param = list(at = c(0.7, 0.78), side = "bottom", labels_rot = "0",gp = gpar(fontsize = 6, fontfamily = "Arial")),
    numbers_gp = gpar(fontsize = 6, fontfamily = "Arial"),
    gp = gpar(col = "dodgerblue4", fill = "dodgerblue2", fontfamily = "Arial"),
    width = unit(0.25, "inches")
  ),
  annotation_name_side = "top",
  annotation_label = c(gap_anno = NULL, Mean_of_Validation_Datasets = "Average\nC-Index\n  of Validation\nDatasets"),
  annotation_name_rot = c(0),
  annotation_name_gp = gpar(fontsize = 6, col = "dodgerblue3", fontface = "bold", fontfamily = "Arial")
)
# top annotation
dataset_type <- c("Training","Validation","Validation","Validation")
dataset_color <- c(
  "TCGA_training" = rgb(231,098,084,maxColorValue = 255,alpha = 220),
  "TCGA_validation" = rgb(255,208,111,maxColorValue = 255,alpha = 220),
  "E_MTAB_1980" = rgb(021,151,165,maxColorValue = 255,alpha = 220),
  "CPTAC"=rgb(055,103,149,maxColorValue = 255,alpha = 220)
)
Top <- columnAnnotation(
  datasetname=anno_simple(dataset_name,border = T,col = dataset_color,gp = gpar(col = "black")),
  show_annotation_name = F)
Top_legend <- Legend(
  at = dataset_name,
  labels = c("TCGA Training", "TCGA Validation", "E-MTAB-1980", "CPTAC"),
  legend_gp = gpar(fill = dataset_color, fontfamily = "Arial"),
  labels_gp = gpar(fontsize = 6, fontfamily = "Arial"),
  title_gp = gpar(fontsize = 7, fontface = "plain", fontfamily = "Arial"),
  title = "Datasets",
  grid_height = unit(0.3, "cm")
)
ht <- Heatmap(
  Selected_C_indexs_main_table,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  col = main_table_color,
  na_col = "white",
  rect_gp = gpar(col = "black"),
  column_names_rot = 0,
  show_column_names = FALSE,
  right_annotation = Right,
  top_annotation = Top,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 6, fontfamily = "Arial"),
  name = "C-Index",
  show_heatmap_legend = FALSE,
  column_split = c(1, rep(2, ncol(Selected_C_indexs_main_table) - 1)),
  column_title = NULL,
  column_gap = unit(2, "mm"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.2f", Selected_C_indexs_main_table[i, j]), x, y, 
              gp = gpar(fontsize = 6, fontfamily = "Arial"))
  }
)
lgd_Heatmap <- Legend(
  title = "C-Index",
  col_fun = main_table_color,
  at = c(0.6, 0.7, 0.8, 0.9),
  labels = c("0.6", "0.7", "0.8", "0.9"),
  labels_gp = gpar(fontsize = 6, fontfamily = "Arial"),
  title_gp = gpar(fontsize = 7, fontface = "plain", fontfamily = "Arial"),
  legend_height = unit(0.5, "cm")
)

c_index_heatmap <- grid.grabExpr({
  draw(ht,padding = unit(c(1, 0, 0, 5), "mm"),
       annotation_legend_list = packLegend(Top_legend,lgd_Heatmap,gap = unit(0.1, "cm"),direction = "horizontal"),
       annotation_legend_side = "bottom")
  decorate_annotation("Mean_of_Validation_Datasets", {
    for (i in 1:(nrow - 1)) {
      grid.segments(
        x0 = unit(Mean_of_All_Datasets_scaled[i], "npc"),
        y0 = unit((nrow - i + 1) / nrow - 0.5 / nrow, "npc"),
        x1 = unit(Mean_of_All_Datasets_scaled[i + 1], "npc"),
        y1 = unit((nrow - i) / nrow - 0.5 / nrow, "npc"),
        gp = gpar(col = "firebrick2", lwd = 1.3)
      )
    }
    grid.text(
      "Average\nC-Index of\nAll Datasets",
      x = unit(0.5, "npc"),
      y = unit(-0.065, "npc"),
      just = c("center", "bottom"),
      gp = gpar(fontsize = 6, col = "firebrick2", fontface = "bold",fontfamily = "Arial")
    )
  })
}) %>% 
  as.ggplot()
qsave(c_index_heatmap,file = paste0(results_dir,"/c_index_heatmap.qs"))
ggsave(filename = paste0(results_dir,"/C_Index_Heatmap.pdf"),width = 2.4,height = 8)




