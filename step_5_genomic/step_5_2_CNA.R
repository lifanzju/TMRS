library(jsonlite)
library(stringr)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggpubr)
library(qs)
library(parallel)
library(data.table)
rm(list = ls());gc()
source("scripts/utils/plot_settings.R")
output_dir <- "results/step_5_genomic/CNV"
if (!dir.exists(output_dir)) {
  dir.create(output_dir,recursive = T)
}
# TCGA
CNV_Input_dir <- "data/genomic/CNV/TCGA/Gene Level Copy Number"
clinical <- qread("results/step_3_model_evaluation_and_nomogram/processed_clinical_data/TCGA.qs")
datasets <- list.files(CNV_Input_dir)
metadata <- jsonlite::fromJSON("data/genomic/CNV/TCGA/metadata.json")  
SampleID <- metadata$submitter_id%>%str_split(pattern = fixed("_"),simplify = T) 
SampleID <- SampleID[,1]%>%str_replace_all(pattern = "-",replacement = ".")
SampleID_FileName <- data.frame(SampleID=SampleID,FileName=metadata$file_name)
CNV_files <- list.files(CNV_Input_dir,pattern = '*.tsv$',recursive = TRUE) 
CNV_files_names <- str_split(CNV_files,pattern  = fixed("/"),simplify = T)[,2]
process_file <- function(i) {
  data <- fread (paste0(CNV_Input_dir, "/", CNV_files[i]))
  data <- dplyr::select(data, gene_id,copy_number)
  colnames(data)[2] <- substr(SampleID_FileName$SampleID[which(SampleID_FileName$FileName == CNV_files_names[i])], start = 1, stop = 12)
  return(data)
}

file_list <- lapply(
  X = 1:length(CNV_files), 
  FUN = process_file
)
mat_data <- Reduce(function(x, y) {
  merged <- merge(x, y, by = "gene_id", all = TRUE)
  return(merged)
}, 
file_list)

mat_data <- data.frame(mat_data,row.names = 1)
Annotation_data <- read.delim(paste0(CNV_Input_dir,"/",CNV_files[1]),fill = T,row.names = 1,header = T)[,c(1:4)]
Annotation_data <- subset(Annotation_data,!chromosome%in%c("chrX","chrY"))
clinical_of_CNV_data <- subset(clinical,CaseID%in%colnames(mat_data))
clinical_of_CNV_data <- arrange(clinical_of_CNV_data,risk_score)
mat_data <- mat_data[rownames(Annotation_data),clinical_of_CNV_data$CaseID]

discrete_mat_data <- apply(mat_data,MARGIN = 2,FUN = function(col){ifelse(col >= 4, "High Amplification",
                                                                          ifelse(col == 3, "Amplification",
                                                                                 ifelse(col==2, "No Change",
                                                                                        ifelse(col==1, "Deletion",
                                                                                               "High Deletion"))))})
combined_data_dir <- "data/genomic/CNV/combined_data"
if (!dir.exists(combined_data_dir)) {
  dir.create(combined_data_dir,recursive = T)
}
write.csv(mat_data,file = file.path(combined_data_dir, "TCGA.csv"))
write.csv(discrete_mat_data,file = file.path(combined_data_dir, "TCGA_CNV_discrete.csv"))
# draw heatmap
main_table_color <- c("Amplification"=rgb(247,170,088,maxColorValue = 255),
                      "High Amplification"=rgb(183,034,048,maxColorValue = 255,alpha = 200),
                      "Deletion"=rgb(049,124,183,maxColorValue = 255),
                      "High Deletion"=rgb(057,081,162,maxColorValue = 255),
                      "No Change"="gray95",
                      "NA"="gray")
Chr_group <- select(Annotation_data, chromosome)
Chr_group <- factor(Chr_group$chromosome,levels = paste0("chr",c(1:22)),ordered = T)
names(Chr_group) <- rownames(Annotation_data)
RS_group <- select(clinical_of_CNV_data,risk_score_group)
RS_group <- factor(RS_group$risk_score_group,levels = c("Low","High"),ordered = T)
names(RS_group) <- clinical_of_CNV_data$CaseID
# top annotation
Alterted_gene_percentages <- apply(discrete_mat_data, MARGIN =2, FUN = function(x){
  No_Change_num <- length(x[x=="No Change"])
  Prop <- (length(x)-No_Change_num)/length(x)*100
  return(Prop)
})
Risk_Score <- select(clinical_of_CNV_data,risk_score)
rownames(Risk_Score) <- clinical_of_CNV_data$CaseID
Altered_Ratio_Color <- colorRamp2(c(0,7,15,45,100),
                                  c(rgb(057,081,162,maxColorValue = 255),
                                    rgb(049,124,183,maxColorValue = 255),
                                    "white",
                                    rgb(247,170,088,maxColorValue = 255),
                                    rgb(183,034,048,maxColorValue = 255)),
                                  transparency = 0.2)
Risk_Score_Color <- colorRamp2(c(-1,-0.5,0,0.5,1),
                               c(rgb(057,081,162,maxColorValue = 255),
                                 rgb(049,124,183,maxColorValue = 255),
                                 "white",
                                 rgb(247,170,088,maxColorValue = 255),
                                 rgb(183,034,048,maxColorValue = 255)),
                               transparency = 0.2)
RS_group_anno <- anno_simple(clinical_of_CNV_data["risk_score_group"],
                             col = c("Low"=rgb(049,124,183,maxColorValue = 255),
                                     "High"=rgb(183,034,048,maxColorValue = 255)),
                             border = T,
                             height = unit(0.2,"cm"))
Risk_Score_anno <- anno_simple(Risk_Score,col = Risk_Score_Color,border = T,height = unit(0.2,"cm"))
Altered_Ratio_anno <- anno_simple(Alterted_gene_percentages,border = T,col = Altered_Ratio_Color,height = unit(0.2,"cm"))
Top_Anno <- columnAnnotation(" "=RS_group_anno,
                             "TMRS"=Risk_Score_anno,
                             "%Genes\nAffected"=Altered_Ratio_anno,
                             annotation_name_side="left",
                             gap=unit(1,"mm"),
                             annotation_name_gp=gpar(fontsize=6))
# prepare legends
main_table_legend <- Legend(at = c("Amplification","High Amplification","Deletion","High Deletion","No Change","NA"),
                            labels = c("Amp\n(CN=3)","High Amp\n(CN>=4)","Del\n(CN=1)","High Del\n(CN=0)","NC\n(CN=2)","NA"),
                            legend_gp = gpar(fill =main_table_color),
                            labels_gp = gpar(fontsize = 6, fontfamily = "Arial"), 
                            title_gp = gpar(fontsize = 7, fontface = "plain", fontfamily = "Arial"),
                            grid_height = unit(0.1, "cm"),
                            grid_width = unit(0.3, "cm"),
                            title = "CNA")
Risk_Score_legend <-  Legend(title = "TMRS", col_fun =  Risk_Score_Color,
                             labels_gp = gpar(fontsize = 6, fontfamily = "Arial"), 
                             title_gp = gpar(fontsize = 7, fontface = "plain", fontfamily = "Arial"),
                             legend_height = unit(0.4, "cm")
                             )
proportion_affected_genes_legend <- Legend(title = "Genes\nAffected (%)", 
                                           col_fun =  Altered_Ratio_Color,
                                           labels_gp = gpar(fontsize = 6, fontfamily = "Arial"),
                                           title_gp = gpar(fontsize = 7, fontface = "plain", fontfamily = "Arial"),
                                           legend_height = unit(0.4, "cm")
                                           )
legend_list <- list(main_table_legend,Risk_Score_legend,proportion_affected_genes_legend)
heatmap_graph <- Heatmap(discrete_mat_data,cluster_rows = F,cluster_columns = F,show_heatmap_legend = F,
                         show_column_names = F,show_row_names = F,border = T,
                         col = main_table_color,row_split = Chr_group,row_title_rot = 0,
                         row_title_gp = gpar(fontsize = 6, fontfamily = "Arial"),
                         column_split = RS_group,column_gap = unit(1,"mm"),
                         column_title = c("Low Risk","High Risk"),
                         column_title_gp = gpar(fontsize = 6, fontfamily = "Arial"),
                         top_annotation = Top_Anno,
                         use_raster = T,raster_quality = 10)
pdf(file = file.path(output_dir,"tcga_heatmap.pdf"),width = 3.25,height = 4)
draw(heatmap_graph,annotation_legend_list=legend_list)
grid.text("TCGA", 
          x = unit(0.45, "npc"), y = unit(0.98, "npc"),
          gp = gpar(fontsize = 8, fontfamily = "Arial", fontface = "bold"))
dev.off()
Alterted_gene_percentages_compare <- data.frame(Percentages=Alterted_gene_percentages,RS_group=clinical_of_CNV_data$risk_score_group)
wilcox.test(Percentages~RS_group,data = Alterted_gene_percentages_compare)
comparisons <- list(comparison=c("Low","High"))
ggplot(Alterted_gene_percentages_compare, aes(x=RS_group, y=Percentages, fill=RS_group)) +
  geom_violin(trim = F,color="white")+
  geom_boxplot(width=0.1,position = position_dodge(0.9))+
  scale_fill_manual(values = c("High"=rgb(231,098,084,maxColorValue = 255),
                               "Low"=rgb(114,188,213,maxColorValue = 255)))+
  labs(title= "TCGA - CNA",y="Proportion of Genes Affected (%)",x="Risk Score") +
  stat_compare_means(comparisons = comparisons,
                     method ="wilcox.test",
                     hide.ns = F,
                     label = "p.format",
                     size = 2)+
  my_plot_theme(legend = "none")
ggsave(filename = file.path(output_dir,"TCGA_CNA_vlnplot.pdf"),width = 1.625,height = 2,device = "pdf")

# CPTAC
rm(list = ls());gc()
source("scripts/utils/plot_settings.R")
output_dir <- "results/step_5_genomic/CNV"
CNV_Input_dir <- "data/genomic/CNV/CPTAC/Gene Level Copy Number"
clinical <- qread("results/step_3_model_evaluation_and_nomogram/processed_clinical_data/CPTAC.qs")
SampleSheet <- read.delim(file = "data/genomic/CNV/CPTAC/CPTAC_Sample_Sheet.tsv")
SampleID_used_in_RS_cal <- read.delim(file = "data/genomic/CNV/CPTAC/CPTAC_SampleID.txt")
combined_data_dir <- "data/genomic/CNV/combined_data"
SampleID_used_in_RS_cal <- SampleID_used_in_RS_cal$SampleID
clinical$CaseID <- str_replace_all(clinical$CaseID,pattern = "-",replacement = ".")
SampleID <- SampleSheet$Sample.ID%>%str_replace_all(pattern = fixed("-"),replacement = ".")
CaseID <- SampleSheet$Case.ID%>%str_split(pattern = fixed(","),simplify = T)
CaseID <- CaseID[,1]%>%str_replace_all(pattern = fixed("-"),replacement = fixed("."))
SampleID_CaseID_FileName <- data.frame(SampleID=SampleID,CaseID=CaseID,FileName=SampleSheet$File.Name)
SampleID_CaseID_FileName <- subset(SampleID_CaseID_FileName,CaseID%in%clinical$CaseID)
rownames(SampleID_CaseID_FileName) <- NULL
processed_CaseID_FileName <- data.frame()
for (i in 1:nrow(SampleID_CaseID_FileName)) {
  info <- SampleID_CaseID_FileName[i,] 
  SampleID_list <- str_split(info$SampleID,pattern = fixed(","),simplify = T)%>%str_trim()
  logi <- SampleID_list%in%SampleID_used_in_RS_cal
  if(T%in%logi) processed_CaseID_FileName <- rbind(processed_CaseID_FileName,info)
}
CNV_files <- list.files(CNV_Input_dir,pattern = '*.tsv$',recursive = TRUE)  #Counts文件夹名
CNV_files_split <- str_split(CNV_files,pattern = fixed("/"),simplify = T)%>%data.frame()
CNV_files_split <- subset(CNV_files_split,X2%in%processed_CaseID_FileName$FileName)
CNV_files_processed <- paste0(CNV_files_split$X1,"/",CNV_files_split$X2)
CNV_files_names <- CNV_files_split$X2


process_file <- function(i) {
  data <- fread (paste0(CNV_Input_dir, "/", CNV_files_processed[i]))
  data <- select(data, gene_id,copy_number)
  colnames(data)[2] <- processed_CaseID_FileName$CaseID[which(processed_CaseID_FileName$FileName==CNV_files_names[i])]
  return(data)
}

file_list <- lapply(
  X = 1:length(CNV_files_processed), 
  FUN = process_file
)
mat_data <- Reduce(function(x, y) {
  merged <- merge(x, y, by = "gene_id", all = TRUE)
  return(merged)
}, 
file_list)
mat_data <- data.frame(mat_data,row.names = 1)
Annotation_data <- read.delim(paste0(CNV_Input_dir,"/",CNV_files[1]),fill = T,row.names = 1,header = T)[,c(1:4)]
Annotation_data <- subset(Annotation_data,!chromosome%in%c("chrX","chrY"))
clinical_of_CNV_data <- subset(clinical,CaseID%in%colnames(mat_data))
clinical_of_CNV_data <- arrange(clinical_of_CNV_data,risk_score)
mat_data <- mat_data[rownames(Annotation_data),clinical_of_CNV_data$CaseID]

discrete_mat_data <- apply(mat_data,MARGIN = 2,FUN = function(col){ifelse(col >= 4, "High Amplification",
                                                                          ifelse(col == 3, "Amplification",
                                                                                 ifelse(col==2, "No Change",
                                                                                        ifelse(col==1, "Deletion",
                                                                                               "High Deletion"))))})
write.csv(mat_data,file = file.path(combined_data_dir, "CPTAC.csv"))
write.csv(discrete_mat_data,file = file.path(combined_data_dir, "CPTAC_CNV_discrete.csv"))
# draw heatmap
main_table_color <- c("Amplification"=rgb(247,170,088,maxColorValue = 255),
                      "High Amplification"=rgb(183,034,048,maxColorValue = 255,alpha = 200),
                      "Deletion"=rgb(049,124,183,maxColorValue = 255),
                      "High Deletion"=rgb(057,081,162,maxColorValue = 255),
                      "No Change"="gray95",
                      "NA"="gray")
Chr_group <- select(Annotation_data, chromosome)
Chr_group <- factor(Chr_group$chromosome,levels = paste0("chr",c(1:22)),ordered = T)
names(Chr_group) <- rownames(Annotation_data)
RS_group <- select(clinical_of_CNV_data,risk_score_group)
RS_group <- factor(RS_group$risk_score_group,levels = c("Low","High"),ordered = T)
names(RS_group) <- clinical_of_CNV_data$CaseID
# top annotation
Alterted_gene_percentages <- apply(discrete_mat_data, MARGIN =2, FUN = function(x){
  No_Change_num <- length(x[x=="No Change"])
  Prop <- (length(x)-No_Change_num)/length(x)*100
  return(Prop)
})
Risk_Score <- select(clinical_of_CNV_data,risk_score)
rownames(Risk_Score) <- clinical_of_CNV_data$CaseID
Altered_Ratio_Color <- colorRamp2(c(0,7,15,45,100),
                                  c(rgb(057,081,162,maxColorValue = 255),
                                    rgb(049,124,183,maxColorValue = 255),
                                    "white",
                                    rgb(247,170,088,maxColorValue = 255),
                                    rgb(183,034,048,maxColorValue = 255)),
                                  transparency = 0.2)
Risk_Score_Color <- colorRamp2(c(-1,-0.5,0,0.5,1),
                               c(rgb(057,081,162,maxColorValue = 255),
                                 rgb(049,124,183,maxColorValue = 255),
                                 "white",
                                 rgb(247,170,088,maxColorValue = 255),
                                 rgb(183,034,048,maxColorValue = 255)),
                               transparency = 0.2)
RS_group_anno <- anno_simple(clinical_of_CNV_data["risk_score_group"],
                             col = c("Low"=rgb(049,124,183,maxColorValue = 255),
                                     "High"=rgb(183,034,048,maxColorValue = 255)),
                             border = T,
                             height = unit(0.2,"cm"))
Risk_Score_anno <- anno_simple(Risk_Score,col = Risk_Score_Color,border = T,height = unit(0.2,"cm"))
Altered_Ratio_anno <- anno_simple(Alterted_gene_percentages,border = T,col = Altered_Ratio_Color,height = unit(0.2,"cm"))
Top_Anno <- columnAnnotation(" "=RS_group_anno,
                             "TMRS"=Risk_Score_anno,
                             "%Genes\nAffected"=Altered_Ratio_anno,
                             annotation_name_side="left",
                             gap=unit(1,"mm"),
                             annotation_name_gp=gpar(fontsize=6))
# prepare legends
main_table_legend <- Legend(at = c("Amplification","High Amplification","Deletion","High Deletion","No Change","NA"),
                            labels = c("Amp\n(CN=3)","High Amp\n(CN>=4)","Del\n(CN=1)","High Del\n(CN=0)","NC\n(CN=2)","NA"),
                            legend_gp = gpar(fill =main_table_color),
                            labels_gp = gpar(fontsize = 6, fontfamily = "Arial"), 
                            title_gp = gpar(fontsize = 7, fontface = "plain", fontfamily = "Arial"),
                            grid_height = unit(0.1, "cm"),
                            grid_width = unit(0.3, "cm"),
                            title = "CNA")
Risk_Score_legend <-  Legend(title = "TMRS", col_fun =  Risk_Score_Color,
                             labels_gp = gpar(fontsize = 6, fontfamily = "Arial"), 
                             title_gp = gpar(fontsize = 7, fontface = "plain", fontfamily = "Arial"),
                             legend_height = unit(0.4, "cm")
)
proportion_affected_genes_legend <- Legend(title = "Genes\nAffected (%)", 
                                           col_fun =  Altered_Ratio_Color,
                                           labels_gp = gpar(fontsize = 6, fontfamily = "Arial"),
                                           title_gp = gpar(fontsize = 7, fontface = "plain", fontfamily = "Arial"),
                                           legend_height = unit(0.4, "cm")
)
legend_list <- list(main_table_legend,Risk_Score_legend,proportion_affected_genes_legend)
heatmap_graph <- Heatmap(discrete_mat_data,cluster_rows = F,cluster_columns = F,show_heatmap_legend = F,
                         show_column_names = F,show_row_names = F,border = T,
                         col = main_table_color,row_split = Chr_group,row_title_rot = 0,
                         row_title_gp = gpar(fontsize = 6, fontfamily = "Arial"),
                         column_split = RS_group,column_gap = unit(1,"mm"),
                         column_title = c("Low Risk","High Risk"),
                         column_title_gp = gpar(fontsize = 6, fontfamily = "Arial"),
                         top_annotation = Top_Anno,
                         use_raster = T,raster_quality = 10)
pdf(file = file.path(output_dir,"cptac_heatmap.pdf"),width = 3.25,height = 4)
draw(heatmap_graph,annotation_legend_list=legend_list)
grid.text("CPTAC", 
          x = unit(0.45, "npc"), y = unit(0.98, "npc"),
          gp = gpar(fontsize = 8, fontfamily = "Arial", fontface = "bold"))
dev.off()
Alterted_gene_percentages_compare <- data.frame(Percentages=Alterted_gene_percentages,RS_group=clinical_of_CNV_data$risk_score_group)
wilcox.test(Percentages~RS_group,data = Alterted_gene_percentages_compare)
comparisons <- list(comparison=c("Low","High"))
ggplot(Alterted_gene_percentages_compare, aes(x=RS_group, y=Percentages, fill=RS_group)) +
  geom_violin(trim = F,color="white")+
  geom_boxplot(width=0.1,position = position_dodge(0.9))+
  scale_fill_manual(values = c("High"=rgb(231,098,084,maxColorValue = 255),
                               "Low"=rgb(114,188,213,maxColorValue = 255)))+
  labs(title= "CPTAC - CNA",y="Proportion of Genes Affected (%)",x="Risk Score") +
  stat_compare_means(comparisons = comparisons,
                     method ="wilcox.test",
                     hide.ns = F,
                     label = "p.format",
                     size = 2)+
  my_plot_theme(legend = "none")
ggsave(filename = file.path(output_dir,"CPTAC_CNA_vlnplot.pdf"),width = 1.625,height = 2,device = "pdf")


