# ==========================================================
# Fig1A.R
# Figure 1A: RNA-seq of CRISPR-Cas9 sgNT and sgLkb1 LKR10 cells treated
# with vehicle or IFNg (16 h, 25 ng/mL). Heatmap of hierarchically
# clustered IFNg-dependent differentially expressed genes
# (FC >= 2.0, FDR <= 0.05; n = 286) in sgNT and sgLkb1 cells.
#
# Requires: run scripts/01_RNAseq_preprocess.R first
# Input:    results/QC/run2_QC.RDS   (run 2, from 01_RNAseq_preprocess.R)
#           data/tableS1.txt         (IFNg-dependent gene list, Table S1)
# Output:   figures/Figure_1A.jpeg
# R version used: 4.3.2
# ==========================================================

library(DESeq2)
library(pheatmap)
library(ggplot2)
library(here)

# Shared helper functions
source(here("R", "plot_helpers.R"))

set.seed(12345678)

dir.create(here("figures"), showWarnings = FALSE)

# ---- Load IFNg-dependent gene list (n = 286; Table S1) ----
ifng_dependent <- read.delim(here("data", "tableS1.txt"),
                             header = TRUE)$gene_name

# ---- Load QC object ----
# [[1]] = DESeq2 object, [[2]] = transformed expression (used for heatmaps)
mouse2 <- readRDS(here("results", "QC", "run2_QC.RDS"))
subdat <- mouse2[[2]][, mouse2[[2]]$KO %in% c("NT1", "LKB1")]

subdf <- as.data.frame(colData(subdat)[, c("guide_KO", "treatment")])
names(subdf) <- c("KO", "Treatment")
subdf$Treatment <- factor(subdf$Treatment, levels = c("v", "F"))   # vehicle, IFNg
subdf$KO        <- factor(subdf$KO, levels = c("sg1.NT1", "sg1.LKB1"))
subdat <- assay(subdat)[rownames(assay(subdat)) %in% ifng_dependent, ]
subdat <- subdat[, rownames(subdf[order(subdf$Treatment, subdf$KO), ])]

zero_var_rows <- apply(subdat, 1, var) == 0
message("Zero-variance genes removed: ", sum(zero_var_rows))
subdat <- subdat[!zero_var_rows, ]
message("Genes plotted: ", nrow(subdat))

ann_colors <- list(
  KO        = c("sg1.NT1" = "grey", "sg1.LKB1" = "darkred"),
  Treatment = c("v" = "turquoise", "F" = "gold2")
)

p1 <- pheatmap(subdat,
               cluster_rows             = TRUE,
               show_rownames            = FALSE,
               show_colnames            = FALSE,
               cluster_cols             = FALSE,
               annotation_col           = subdf,
               annotation_colors        = ann_colors,
               main                     = "IFNg response genes in sgLKB1 and NT Sample\n(n=286)",
               scale                    = "row",
               clustering_distance_rows = "euclidean",
               clustering_method        = "ward.D2",
               treeheight_row           = 25,
               fontsize                 = 8,
               color = colorRampPalette(c("#3F51A3", "white", "#EE2B2B"))(50))

save_pheatmap_jpeg(p1, here("figures", "Figure_1A.jpeg"), width = 6, height = 12)
