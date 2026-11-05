# FigS1AB

# Figure S1. LKB1 loss impairs the majority of IFN𝛾-responsive genes in human and mouse lung adenocarcinoma
# (A) RNA-seq analysis of sgNT LKR10 cells treated with vehicle or IFNg (16h, 25ng/mL). Heatmap of hierarchically clustered IFNg-dependent,
# differentially expressed genes “DEG” (FC ³ 2.0, FDR ≤ 0.05). (B) Volcano plot depicting all DE genes between sgNT IFN𝛾 versus sgLkb1 IFN𝛾. Red
# dots represent genes significantly induced by loss of LKB1. Blue dots represent genes significantly downregulated by loss of LKB1.

# R barrel 4.3.2
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(here)

# Shared helper functions
source(here("R", "plot_helpers.R"))

# ---- Figure S1A ----
set.seed(12345678)

mouse2 <- readRDS("run2_QC.RDS")
# ---- Load IFNg-dependent gene list (n = 286; Table S1) ----
ifng_dependent <- read.delim(here("data", "tableS1.txt"),
                             header = TRUE)$gene_name

subdat <- mouse2[[2]][,mouse2[[2]]$KO %in% c("NT1")]
subdf <- as.data.frame(colData(subdat)[,  c("KO", "treatment")])
names(subdf) <- c("KO","Treatment")

subdat <- assay(subdat)[rownames(assay(subdat)) %in% ifng_dependent, ]
subdf$Treatment <- factor(subdf$Treatment, levels=c("v", "F"))

subdat <- subdat[,rownames(subdf[order(subdf$KO, subdf$Treatment),])]
zero_var_rows <- apply(subdat, 1, var) == 0
subdat <- subdat[!zero_var_rows, ]

ann_colors <- list(
  KO = c("NT1" = "grey"),
  Treatment = c( "v" = "turquoise", "F" = "gold2")
)

p1 <- pheatmap::pheatmap(subdat,
                         cluster_rows = T,
                         show_rownames = F,
                         show_colnames = F,
                         cluster_cols = F,
                         annotation_col = subdf,
                         annotation_colors = ann_colors,
                         main = "IFNg response genes in NT Sample",
                         scale = "row",
                         clustering_distance_rows = "euclidean",
                         clustering_method = "ward.D2",
                         treeheight_row = 25,
                         fontsize = 8,
                         color = colorRampPalette(c("#3F51A3", "white", "#EE2B2B"))(50))
save_pheatmap_jpeg(p1, "FigS1A.jpeg",width=6,height=10)

# ---- Figure S1B ----
set.seed(12345678)

deresult <- readRDS("run2_results.RDS")

geneList=c("Tmem173",  "Irf7",  "Stat1", "Stat2", "Cd274") # "Irf1","Irf9" not signif
name="sgNT IFNg VS sgLKB1 IFNg"
gV <- gVolcano4review(deresult[[6]],geneList = geneList,name=name)
ggsave("FigS1B.pdf", gV, width = 6, height = 8, units = "in")
