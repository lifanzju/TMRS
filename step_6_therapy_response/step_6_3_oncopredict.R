library(oncoPredict)
library(data.table)
library(gtools)
library(reshape2)
library(ggpubr)
library(dplyr)
library(stringr)
library(pls)
library(cowplot)
library(scales)
library(ggforce)

rm(list = ls());gc()
options(stringsAsFactors = F)
output_dir <- "results/step_6_therapy_response/oncopredict"
if(!dir.exists(output_dir)) dir.create(output_dir,recursive = T)


exp_dir <- "data/Datasets/proprocessed_exp"
datasets <- list.files(exp_dir)%>%str_remove_all(fixed(".csv"))
GDSC_Expr <- readRDS(file = "data/oncopredict/GDSC1_Expr (RMA Normalized and Log Transformed).rds")
GDSC_Res <- readRDS(file = "data/oncopredict/GDSC1_Res.rds")%>%data.frame()
drugs <- c("Pazopanib_199","Axitinib_1021","Sorafenib_30","Cabozantinib_249","Sunitinib_5","Temsirolimus_1016") 
GDSC_Res <- select(GDSC_Res,all_of(drugs))%>%
  as.matrix()%>%
  exp()
datasets <- c("TCGA","CPTAC","E_MTAB_1980")

work_dir <- getwd()
for (dataset in datasets) {
  result_dir <- paste0(output_dir,"/",dataset)
  if(!dir.exists(result_dir)) dir.create(result_dir)
  data <- fread(paste0(exp_dir,"/",dataset,".csv"),header = T)%>%
    data.frame(row.names = 1)%>%
    as.matrix()
  setwd(result_dir)
  calcPhenotype(
    trainingExprData = GDSC_Expr,
    trainingPtype = GDSC_Res,
    testExprData = data,
    powerTransformPhenotype = T,
    printOutput = TRUE,
    rsq = F,
    cc=F,
    pcr=T,
    batchCorrect = ifelse(dataset%in%c("CPTAC","TCGA"),"standardize","eb"),
    minNumSamples = 10,
    removeLowVaringGenesFrom = "rawData"
  )
  setwd(work_dir)
}



clinical_dir <- "results/step_3_model_evaluation_and_nomogram/processed_clinical_data"
graph_dir <- file.path(output_dir,"graph")
if(!dir.exists(graph_dir)) dir.create(graph_dir)
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

all_result <- data.frame()
for (dataset in datasets){
  prediction <- read.csv(file = file.path(output_dir,dataset,"calcPhenotype_Output","DrugPredictions.csv"),row.names = 1) %>% 
    log() %>% 
    select(drugs) %>% 
    t()
  clinical <- qread(file.path(clinical_dir,paste0(dataset,".qs")))
  prediction <- prediction[,clinical$CaseID]
  risk_score <- select(clinical,"risk_score") %>% t()
  analysis_table <- rbind(risk_score,prediction)
  correlation <- cor_fun(analysis_table,method = "pearson")
  result <- cbind(Dataset=dataset,Drug=str_split(row.names(correlation),pattern = fixed("_"),simplify = T)[,1],correlation)
  all_result <- rbind(all_result,result)
}
all_result$logP <- -log10(all_result$pvalue)
all_result$sig <- ifelse(all_result$pvalue<0.05,"sig","not")
write.csv(all_result,file = paste0(output_dir,"/","cor_result_all_datasets.csv"))

all_result$Dataset <- str_replace_all(all_result$Dataset,pattern = fixed("_"),replacement = "-")
all_result$Drug <- factor(all_result$Drug,levels = c("Cabozantinib","Sorafenib","Temsirolimus","Sunitinib","Axitinib","Pazopanib"))
ggplot(all_result,aes(x=Drug, y= Dataset,fill=cor,size=logP,color = sig)) + 
  geom_point(shape =21)+
  scale_fill_gradientn(name = "Cor",
                       colours = c("#1663A9",rgb(033,158,188,maxColorValue = 255),
                                   "white",rgb(255,183,005,maxColorValue = 255),"darkorange"),
                       values = rescale(c(min(all_result$cor),-0.25,0,0.1,max(all_result$cor))),
                       breaks = c(-0.4,-0.15,0,0.1))+
  scale_size_continuous(name = "p-value\n(-log10)",
                        limit = c(0,max(all_result$logP)),
                        breaks = c(0,5,10),
                        range = c(2,5))+
  labs(x= NULL,
       y= NULL,
       title = "Correlation with predicted IC50(log)")+
  scale_color_manual(values = c("sig" = "black","not"="gray"))+
  my_plot_theme(x.text.angle = 45)+
  theme(panel.grid.major = element_line(colour = "gray"),
        legend.key.height = unit(0.2,"cm"),
        legend.spacing.y = unit(-0.4,"cm"))
ggsave(filename = paste0(graph_dir,"/all_cor.pdf"),width = 3.25,height = 2,device = "pdf")
