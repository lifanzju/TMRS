library(qs)
library(data.table)
library(survival)
library(survminer)
library(dplyr)
library(stringr)
library(rms)

rm(list = ls())
source("scripts/utils/plot_settings.R")
exp_dir <- "data/data_for_ml/exp"
clinical_dir <- "results/step_3_model_evaluation_and_nomogram/processed_clinical_data"
result_dir <- "results/step_3_model_evaluation_and_nomogram/nomogram_construction"
if (!dir.exists(result_dir)) dir.create(result_dir)
datasets <- str_remove_all(list.files(exp_dir), pattern = fixed(".csv"))
dataset_names <- list(
  "TCGA" = "TCGA",
  "E_MTAB_1980" = "E-MTAB-1980",
  "CPTAC" = "CPTAC",
  "TCGA_training" = "TCGA Training",
  "TCGA_validation" = "TCGA Validation"
)

tcga_clinical <- qread(file.path(clinical_dir, "TCGA.qs"))
colnames(tcga_clinical)[10] <- "TMRS"

# os nomogram
ddDD=datadist(tcga_clinical)
options(datadist="ddDD")
os_coxresult <- cph(Surv(OS,OS.Censor)~TMRS+Age+Stage,surv=T,x=T, y=T,data = tcga_clinical)
os_nomogram_cox <- coxph(Surv(OS,OS.Censor)~TMRS+Age+Stage,data = tcga_clinical)
os_surv <- Survival(os_coxresult)
sur_1_year<-function(x) os_surv(1*365*1,lp=x)  
sur_3_year<-function(x) os_surv(1*365*3,lp=x)
sur_5_year<-function(x) os_surv(1*365*5,lp=x)
os_nomogram <- nomogram(os_coxresult,fun=list(sur_1_year,sur_3_year,sur_5_year),
                        lp= T,
                        funlabel=c('Pr(OS < 1-year)','Pr(OS < 3-year)','Pr(OS < 5-year)'),
                        maxscale=100,
                        fun.at=c('0.9','0.7','0.5','0.3','0.1'))
pdf(file = paste0(result_dir, "/os_nomogram.pdf"), width = 3.25, height = 2, fonts = "Arial",pointsize = 7)
par(mar = c(0.5, 0, 0.5, 0))
plot(os_nomogram, gp = gpar(fontfamily = "Arial"))
dev.off()
qsave(os_nomogram_cox, file = paste0(result_dir, "/os_nomogram_cox.qs"))
# pfs nomogram
label(tcga_clinical$grade_group) <- "Grade Group"
options(contrasts = c("contr.treatment", "contr.treatment"))
pfs_coxresult <- cph(Surv(PFS,PFS.Censor)~TMRS+Gender+Stage+grade_group,surv=T,x=T, y=T,data = tcga_clinical)
pfs_nomogram <- coxph(Surv(PFS,PFS.Censor)~TMRS+Gender+Stage+grade_group,data = tcga_clinical)
pfs_surv <- Survival(pfs_coxresult)
sur_1_year<-function(x) pfs_surv(1*365*1,lp=x)  
sur_3_year<-function(x) pfs_surv(1*365*3,lp=x)
sur_5_year<-function(x) pfs_surv(1*365*5,lp=x)

nom_sur_pfs <- nomogram(pfs_coxresult,fun=list(sur_1_year,sur_3_year,sur_5_year),
                        lp= T,
                        funlabel=c('Pr(PFS < 1-year)','Pr(PFS < 3-year)','Pr(PFS < 5-year)'),
                        maxscale=100,
                        fun.at=c('0.9','0.7','0.5','0.3','0.1'))
pdf(file = paste0(result_dir, "/pfs_nomogram.pdf"), width = 3.25, height = 2,pointsize = 7, fonts = "Arial")
par(mar = c(0.5, 0, 0.5, 0))
plot(nom_sur_pfs,gp = gpar(fontfamily = "Arial"))
dev.off()
qsave(pfs_coxresult, file = paste0(result_dir,"/pfs_nomogram_cox.qs"))
