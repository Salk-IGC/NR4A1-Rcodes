# ==========================================================
# fig2D.R
# Figure 2. The LKB1-CRTC-CREB signaling axis dominantly controls
# interferon response.
# (D) RNA-seq analysis of CRISPR-Cas9-mediated sgNT, sgLkb1,
# sgLkb1+sgCrtc2, and sgLkb1+sgCreb1 LKR10 cells treated with vehicle or
# IFNg (16 h, 25 ng/mL). Heatmap of hierarchically clustered
# LKB1-dependent IFNg response genes (n = 166) in sgNT, sgLkb1, sgCrtc2,
# and sgCreb1 LKR10 cells.
#
# Requires: run scripts/01_RNAseq_preprocess.R first
# Input:    results/QC/run1_QC.RDS   (run 1, from 01_RNAseq_preprocess.R)
#           data/tableS2.txt         (LKB1-dependent gene list, Table S2)
# Output:   figures/Figure_2D.jpeg
# Samples:  sgNT, sgLKB1, sg2.CREB1, sgCRTC2 (sg1.CREB1 excluded);
#           vehicle, then IFNg
# ==========================================================

library(DESeq2)
library(pheatmap)
library(ggplot2)
library(here)

# Shared helper functions
source(here("R", "plot_helpers.R"))

set.seed(12345678)

dir.create(here("figures"), showWarnings = FALSE)

# ---- Load LKB1-dependent gene list (n = 166, Table S2) ----
lkb1_dependent <- read.delim(here("data", "tableS2.txt"),
                             header = TRUE)$gene_name

# ---- Load QC object (run 1) ----
# [[1]] = DESeq2 object, [[2]] = transformed expression (used for heatmaps)
mouse1 <- readRDS(here("results", "QC", "run1_QC.RDS"))
subdat <- mouse1[[2]]

# ---- Sample annotation: drop sg1.CREB1, order vehicle then IFNg ----
subdf <- as.data.frame(colData(subdat)[, c("guide_KO", "treatment")])
names(subdf) <- c("KO", "Treatment")
subdf <- subdf[subdf$KO != "sg1.CREB1", ]
subdf$Treatment <- factor(subdf$Treatment, levels = c("v", "F"))   # vehicle, IFNg
subdf$KO <- factor(subdf$KO, levels = c("sg1.NT", "sg1.LKB1", "sg2.CREB1",
                                        "sg1.CRTC2", "sg2.CRTC2"))
subdf <- subdf[order(subdf$Treatment, subdf$KO), ]
subdat <- assay(subdat)[rownames(assay(subdat)) %in% lkb1_dependent, ]
subdat <- subdat[, -c(1:6)]   # drops the first 6 columns (sg1.CREB1 samples)

zero_var_rows <- apply(subdat, 1, var) == 0
message("Zero-variance genes removed: ", sum(zero_var_rows))
subdat <- subdat[!zero_var_rows, ]
message("Genes plotted: ", nrow(subdat))
subdat <- subdat[, rownames(subdf)]

ann_colors <- list(
  Treatment = c("v" = "turquoise", "F" = "gold2"),
  KO = c("sg1.NT"    = "grey",
         "sg1.LKB1"  = "darkred",
         "sg2.CREB1" = "coral",
         "sg1.CRTC2" = "#00B3FF",
         "sg2.CRTC2" = "#00B3FF")
)

p1 <- pheatmap(subdat,
               cluster_rows             = TRUE,
               show_rownames            = FALSE,
               show_colnames            = FALSE,
               cluster_cols             = FALSE,
               annotation_col           = subdf,
               annotation_colors        = ann_colors,
               main                     = "LKB1-dependent genes\n(n=166)",
               scale                    = "row",
               clustering_distance_rows = "euclidean",
               clustering_method        = "ward.D2",
               treeheight_row           = 25,
               fontsize                 = 8,
               color = colorRampPalette(c("#3F51A3", "white", "#EE2B2B"))(50))

save_pheatmap_jpeg(p1, here("figures", "Figure_2D.jpeg"), width = 6, height = 12)
