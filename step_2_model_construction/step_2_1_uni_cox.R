library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(stringr)
library(venn)
library(qs)
library(data.table)

exp_dir <- "data/Datasets/EXP"
clinical_dir <- "data/Datasets/Clinical"
sample_type_dir <- "data/Datasets/SampleType"
wgcna_hdwgcna_intersection <- read.csv(file = "results/step_1_hdwgcna_wgcna/intersection/wgcna_hdwgcna_intersection.csv")[, 1]
results_dir <- "results/step_2_model_construction/uni_cox"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

datasets <- str_remove_all(list.files(exp_dir), pattern = fixed(".txt"))
datasets_including_normal <- c("TCGA", "CPTAC")

# preprocess expression data
preprocessed_exp_dir <- "data/Datasets/proprocessed_exp"
if (!dir.exists(preprocessed_exp_dir)) dir.create(preprocessed_exp_dir, recursive = TRUE)
for (dataset in datasets) {
  exp <- fread(file = paste0(exp_dir, "/", dataset, ".txt"), header = T) %>% 
    data.frame(row.names = 1)
  exp <- t(exp) %>% data.frame()
  exp$SampleID <- rownames(exp)
  if (dataset %in% datasets_including_normal) {
    sample_type <- read.delim(file = paste0(sample_type_dir, "/", dataset, "_SampleType.txt"), header = T, row.names = NULL)
    tumor_samples <- subset(sample_type, Type == "Tumor")[, 1]
    tumor_exp <- subset(exp, SampleID %in% tumor_samples)
  } else {
    tumor_exp <- exp
  }
  if (dataset == "CPTAC") {
    case_id <- strsplit(tumor_exp$SampleID, split = ".", fixed = T)
    case_id <- sapply(case_id, function(x) {
      paste(x[1], x[2], sep = ".")
    })
    tumor_exp$CaseID <- case_id
    tumor_exp <- tumor_exp[, -(ncol(tumor_exp) - 1)]
  }
  if (dataset == "TCGA") {
    case_id <- strsplit(tumor_exp$SampleID, split = ".", fixed = T)
    case_id <- sapply(case_id, function(x) {
      paste(x[1], x[2], x[3], sep = ".")
    })
    tumor_exp$CaseID <- case_id
    tumor_exp <- tumor_exp[, -(ncol(tumor_exp) - 1)]
  }
  if (dataset == "E_MTAB_1980") {
    case_id <- tumor_exp$SampleID
    tumor_exp$CaseID <- case_id
    tumor_exp <- tumor_exp[, -(ncol(tumor_exp) - 1)]
  }
  tumor_exp <- distinct(tumor_exp, CaseID, .keep_all = T)
  rownames(tumor_exp) <- tumor_exp$CaseID
  tumor_exp <- tumor_exp[, colnames(tumor_exp) != "CaseID"]
  write.csv(t(tumor_exp), file = paste0(preprocessed_exp_dir, "/", dataset, ".csv"))
}

# univariate cox regression
all_cox_results_dir <- paste0(results_dir, "/all")
sig_cox_results_dir <- paste0(results_dir, "/sig")
if (!dir.exists(all_cox_results_dir)) dir.create(all_cox_results_dir, recursive = TRUE)
if (!dir.exists(sig_cox_results_dir)) dir.create(sig_cox_results_dir, recursive = TRUE)

for (dataset in datasets) {
  exp <- read.csv(file = paste0(preprocessed_exp_dir, "/", dataset, ".csv"), header = T, row.names = 1) %>%
    t() %>%
    data.frame()
  clinical <- read.delim(file = paste0(Clinical.dir, "/", dataset, "_Clinical.txt"), header = T, row.names = 1)
  coef <- c()
  HR <- c()
  HR.95L <- c()
  HR.95H <- c()
  pvalue <- c()
  for (gene in wgcna_hdwgcna_intersection) {
    gene_EXP <- subset(exp, select = gene)
    Analysis_table <- clinical
    Analysis_table$gene <- gene_EXP[match(rownames(Analysis_table), rownames(gene_EXP)), gene]
    coxresult <- summary(coxph(Surv(OS, OS.Censor) ~ gene, data = Analysis_table))
    coef <- c(coef, coxresult$coefficients["gene", "coef"])
    HR <- c(HR, coxresult$conf.int["gene", "exp(coef)"])
    HR.95L <- c(HR.95L, coxresult$conf.int["gene", "lower .95"])
    HR.95H <- c(HR.95H, coxresult$conf.int["gene", "upper .95"])
    pvalue <- c(pvalue, coxresult$coefficients["gene", "Pr(>|z|)"])
  }
  univariate_cox <- data.frame(coef = coef, HR = HR, HR.95L = HR.95L, HR.95H = HR.95H, pvalue = pvalue)
  row.names(univariate_cox) <- wgcna_hdwgcna_intersection
  univariate_cox_sig <- subset(univariate_cox, pvalue < 0.05)
  write.csv(univariate_cox, file = paste0(all_cox_results_dir, "/", dataset, "_cox_regression_all.csv"))
  write.csv(univariate_cox_sig, file = paste0(sig_cox_results_dir, "/", dataset, "_cox_regression_sig.csv"))
}

# take the intersection of significant genes
intersection_list <- list()
for (dataset in datasets) {
  intersection_list[[dataset]] <- read.csv(file = paste0(sig_cox_results_dir, "/", dataset, "_cox_regression_sig.csv"), header = T, row.names = NULL)[, 1]
}
names(intersection_list)[names(intersection_list) == "E_MTAB_1980"] <- "E-MTAB-1980"
intersects <- function(x) {
  Reduce(intersect, x)
}
intersection <- intersects(intersection_list)
write.csv(intersection, file = paste0(results_dir, "/sig_uni_cox_intersection.csv"), col.names = F, row.names = F)
venn_plot <- venn(intersection_list,
  zcolor = c(
    rgb(114, 188, 213, maxColorValue = 256),
    rgb(253, 217, 133, maxColorValue = 256),
    rgb(231, 098, 084, maxColorValue = 256)
  ),
  opacity = 0.6,
  ilabels = "counts",
  box = F,
  ilcs = 0.5,
  sncs = 0.5,
  ggplot = T
) +
  coord_fixed(ratio = 1)
qsave(venn_plot, file = paste0(results_dir, "/sig_uni_cox_venn_plot.qs"))
ggsave(file = paste0(results_dir, "/sig_uni_cox_venn_plot.tiff"), plot = venn_plot, width = 3.25, height = 2.67)




