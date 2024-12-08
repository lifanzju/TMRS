library(venn)
library(clusterProfiler)
library(org.Hs.eg.db)
library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(qs)
library(stringr)
rm(list = ls())
source("scripts/utils/plot_settings.R")

results_dir <- "results/step_1_hdwgcna_wgcna/intersection"
if (!dir.exists(results_dir)) dir.create(results_dir)

# take intersection of WGCNA and hdWGCNA
WGCNA <- read.csv(
  file = "results/step_1_hdwgcna_wgcna/WGCNA/GeneOfEachModules/Selected/Selected_Gene_Information_of_MMpurple.csv",
  header = T, row.names = NULL
)[, 1]
hdWGCNA <- read.csv(file = "results/step_1_hdwgcna_wgcna/hdWGCNA/selected_gene_information.csv", header = T, row.names = NULL)[, 1]
wgcna_hdwgcna_intersection <- intersect(WGCNA, hdWGCNA)
vennlist <- list(TCGA = WGCNA, GSE156632 = hdWGCNA)
wgcna_hdwgcna_venn_plot <- venn(vennlist,
  zcolor = c(
    rgb(114, 188, 213, maxColorValue = 256),
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
qsave(wgcna_hdwgcna_venn_plot,file.path(results_dir, "venn_plot.qs"))
ggsave(
  file.path(results_dir, "/venn_plot.tiff"),
  wgcna_hdwgcna_venn_plot,
  width = 3.25,
  height = 2.67
)
write.csv(wgcna_hdwgcna_intersection, file = file.path(results_dir, "wgcna_hdwgcna_intersection.csv"), row.names = F)

# function analysis of intersection
GO_BP_genesets <- read.gmt("data/genesets/other_genesets/c5.go.bp.v2023.1.Hs.symbols.txt")
GO_BP_result <- enricher(gene = wgcna_hdwgcna_intersection,TERM2GENE = GO_BP_genesets)@result
GO_BP_result$logP <- -log10(GO_BP_result$p.adjust)
GO_BP_result <- arrange(GO_BP_result,desc(logP))
GO_BP_result <- GO_BP_result[c(1:10),]
GeneRatio <- str_split(GO_BP_result$GeneRatio,pattern = fixed("/"))%>%
        sapply(FUN = function(x){as.numeric(x[1])/as.numeric(x[2])})
BgRatio <- str_split(GO_BP_result$BgRatio,pattern = fixed("/"))%>%
        sapply(FUN = function(x){as.numeric(x[1])/as.numeric(x[2])})
Proper_name <- function(data,geneset){
        data <- str_remove_all(data,pattern = paste0(geneset,"_"))
        data <- str_replace_all(data,pattern = "_",replacement = " ")
        data <- str_to_title(data)
        data <- str_replace_all(data,pattern = "To",replacement = "to")
        data <- str_replace_all(data,pattern = "Atp", replacement = "ATP")
        data <- str_replace_all(data,pattern = "Tca", replacement = "TCA")
        data <- str_replace_all(data,pattern = "Of", replacement = "of")
        data <- str_replace_all(data,pattern = "On", replacement = "on")
        data <- str_replace_all(data,pattern = "Or", replacement = "or")
        data <- str_replace_all(data,pattern = "Nad P H", replacement = "NADPH")
        data <- str_replace_all(data,pattern = "Nad", replacement = "NAD")
        data <- str_replace_all(data,pattern = "NADh", replacement = "NADH")
        data <- str_replace_all(data,pattern = "As", replacement = "as")
        data <- str_replace_all(data,pattern = "Il", replacement = "IL")
        data <- str_replace_all(data,pattern = "Nfkb", replacement = "NFKB")
        data <- str_replace_all(data,pattern = "E2f", replacement = "E2F")
        data <- str_replace_all(data,pattern = "G2m", replacement = "G2M")
        data <- str_replace_all(data,pattern = "Tgf", replacement = "TGF")
        data <- str_replace_all(data,pattern = "Kras", replacement = "KRAS")
        data <- str_replace_all(data,pattern = "Uv", replacement = "UV")
        data <- str_replace_all(data,pattern = "Jak", replacement = "JAK")
        data <- str_replace_all(data,pattern = "Stat", replacement = "STAT")
        data <- str_replace_all(data,pattern = "Mtorc", replacement = "mTORC")
        data <- str_replace_all(data,pattern = "Pi3k", replacement = "PI3K")
        data <- str_replace_all(data,pattern = "Mtor", replacement = "mTOR")
        data <- str_replace_all(data,pattern = "Dna", replacement = "DNA")
        data <- str_replace_all(data,pattern = "Via", replacement = "via")
        data <- str_replace_all(data,pattern = "Tnfa", replacement = "TNFA")
        data <- str_replace_all(data,pattern = "Rrna", replacement = "rRNA")
        data <- str_replace_all(data,pattern = "Mrna", replacement = "mRNA")
        data <- str_replace_all(data,pattern = "Utr", replacement = "UTR")
        data <- str_replace_all(data,pattern = "Ecm", replacement = "ECM")
        data <- str_replace_all(data,pattern = "assembly", replacement = "Assembly")
        data <- str_replace_all(data,pattern = "And", replacement = "and")
        data <- str_replace_all(data,pattern = "organic", replacement = "Organic")
}
GO_BP_result$RichFactor <- GeneRatio/BgRatio
GO_BP_result$ID <- Proper_name(GO_BP_result$ID,geneset="GOBP")
GO_BP_result$ID <- str_wrap(GO_BP_result$ID, width = 15)
GO_BP_result$ID <- factor(GO_BP_result$ID,levels = rev(GO_BP_result$ID),ordered = T)

GO_BP_plot <- ggplot(GO_BP_result,aes(x=logP, y= ID, size= Count, col=RichFactor,fill=RichFactor)) + geom_point(shape=21)+ 
        theme_bw()+ 
        scale_color_gradient(low = rgb(255,208,111,maxColorValue = 255,alpha=200), high = "firebrick3")+
        scale_fill_gradient(low = rgb(255,208,111,maxColorValue = 255,alpha = 200), high = "firebrick3")+
        scale_size_continuous(breaks = c(15,20,30), range = c(2,6))+
        scale_x_continuous(limits = c(0,max(GO_BP_result$logP)+2))+
        labs(y=NULL,x="-log10 P.adj",title ="GO - BP",fill="Rich\nFactor",color="Rich\nFactor")+
        my_plot_theme()+
        theme(legend.position = c(0.8,0.35),
              legend.key.height = unit(0.25, "cm")) 
qsave(GO_BP_plot,file.path(results_dir, "GO_BP_plot.qs"))
ggsave(
  file.path(results_dir, "/GO_BP_plot.tiff"),
  GO_BP_plot,
  width = 3.25,
  height = 2.67*2
)

