# Figure S6G

library(ggplot2)
library(pheatmap)
library(here)

# Shared helper functions
source(here("R", "plot_helpers.R"))
set.seed(12345678)

deg <- read.delim("deg from table S2.txt")$gene_name

mouse2 <- readRDS(here("results", "QC", "run2_QC.RDS"))
subdat <- mouse2[[2]]
subdat <- mouse2[[2]][,mouse2[[2]]$treatment =="v"]
subdf <- as.data.frame(colData(subdat)[,  c("guide_KO", "treatment")])
names(subdf) <- c("KO","Treatment")

subdat <- assay(subdat)[rownames(assay(subdat)) %in% deg, ]
dim(subdat)
subdf$KO <- factor(subdf$KO, levels=c("sg1.NT1", "sg1.LKB1", "sg1.Nr4a1", "sg2.Nr4a1"))

subdat <- subdat[,rownames(subdf[order(subdf$Treatment, subdf$KO),])]
colnames(subdat)
"Tmem173" %in% rownames(subdat)

zero_var_rows <- apply(subdat, 1, var) == 0
sum(zero_var_rows)  # Count zero variance rows
subdat <- subdat[!zero_var_rows, ]
colnames(subdat)

ann_colors <- list(
  KO = c("sg1.NT1" = "grey", "sg1.LKB1" = "#FF004D", "sg1.Nr4a1"="#CC00FF", "sg2.Nr4a1"="#00FF66"),
  Treatment = c("v" = "turquoise","F" = "gold2")
)

p1 <- pheatmap::pheatmap(subdat,
                         cluster_rows = T,
                         show_rownames = T,
                         show_colnames = F,
                         cluster_cols = F,
                         annotation_col = subdf,
                         annotation_colors = ann_colors,
                         scale = "row",
                         clustering_distance_rows = "euclidean",
                         clustering_method = "ward.D2",
                         treeheight_row = 25,
                         fontsize = 8,
                         color = colorRampPalette(c("blue", "white", "red"))(50))
save_pheatmap_jpeg(p1, "FigS6.jpeg",width=6,height=29)
