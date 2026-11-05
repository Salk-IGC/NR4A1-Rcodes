###FigS5EG

# Figure S5. NR4A1 is a CRTC2-CREB target gene upregulated in LKB1-deficient lung cancer cells
# (E) Volcano plot depicting all DEGs between sgLkb1 versus sgLkb1+sgCrtc2
# stimulated with IFNg. Red dots represent genes significantly induced by loss of Crtc2. Blue dots represent genes significantly downregulated
# by loss of Crtc2. Labeled genes are transcription factors elevated in sgLkb1 cells.
# (G) Heatmap of all expressed Nuclear Receptor family members. Note
# only NR4A1 and NR4A2 are up in LKB1-mutant cells which is reversed by Crtc2 or Creb1 co-deletion. **P<0.01; **** P< 0.0001; ns, not
# significant; unpaired t test.

# R barrel 4.3.2
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(here)

# Shared helper functions
source(here("R", "plot_helpers.R"))

# ---- Figure S5E ----
set.seed(12345678)

deresult <- readRDS("run1_results.RDS")

geneList3=c("Nr4a1", "Mafb", "Irx3", "Irf8", "Foxf1", "Msantd3")
name="sgLKB1 IFNg VS sgCRTC2 IFNg"
gV <- gVolcano4review(deresult[[8]],geneList = geneList3,name=name)
ggsave("FigS5E.pdf", gV, width = 6, height = 8, units = "in")

# ---- Figure S5G ----
set.seed(12345678)

nuc.recep.genes <- c("THRA","THRB","RARA","RARB","RARG","PPARA","PPARD","PPARG","NR1D1","NR1D2",
                     "RORA","RORB","RORC","NR1H3","NR1H2","NR1H4","VDR","NR1I2","NR1I3","HNF4A",
                     "HNF4G","RXRA","RXRB","NR2C1","NR2C2","NR2F1","NR2F2","NR2F6","ESR1","ESR2",
                     "ESRRA","ESRRB","ESRRG","NR3C1","NR3C2","AR","NR4A1","NR4A2","NR4A3","NR5A1",
                     "NR6A1","NR0B2")
nr_genes <- paste0(toupper(substr(nuc.recep.genes, 1, 1)),
                   tolower(substring(nuc.recep.genes, 2)))

mouse1 <- readRDS("run1_QC.RDS")
dds <- mouse1[[2]]

# ---- vehicle samples only ----
subdat <- dds[, dds$treatment == "v"]
subdf  <- as.data.frame(colData(subdat)[, c("guide_KO", "treatment")])
names(subdf) <- c("KO", "treatment")
subdf <- subdf[subdf$KO != "sg1.CREB1", ]

ko_levels <- c("sg1.NT", "sg1.LKB1", "sg1.CRTC2", "sg2.CRTC2", "sg2.CREB1")
ko_labels <- c("sgNT", "sgLkb1", "sgLkb1+sgCrtc2", "sgLkb1+sgCrtc2", "sgLkb1+sgCreb1")
subdf$KO <- factor(subdf$KO, levels = ko_levels)
subdf$KO <- factor(ko_labels[as.integer(subdf$KO)],
                   levels = c("sgNT", "sgLkb1", "sgLkb1+sgCrtc2", "sgLkb1+sgCreb1"))
subdf$treatment <- factor("veh")
subdf <- subdf[order(subdf$KO), c("treatment","KO")]

# ---- matrix: nuclear receptors only ----
mat <- assay(subdat)
mat <- mat[rownames(mat) %in% nr_genes, rownames(subdf)]
mat <- mat[apply(mat, 1, var) > 0, ]
dim(mat)

ann_colors <- list(
  treatment = c(veh = "#B5B82A"),
  KO = c("sgNT"            = "#40C4FF",
         "sgLkb1"          = "#4CD694",
         "sgLkb1+sgCrtc2"  = "#EE82EE",
         "sgLkb1+sgCreb1"  = "#F28B82")
)

z <- t(scale(t(mat)))
lim <- ceiling(max(abs(z), na.rm = TRUE))
brks <- seq(-lim, lim, length.out = 101)

p1 <- pheatmap::pheatmap(mat,
                         scale = "row",
                         breaks = brks,
                         color = colorRampPalette(c("blue", "white", "red"))(100),
                         cluster_rows = TRUE,
                         clustering_distance_rows = "euclidean",
                         clustering_method = "ward.D2",
                         treeheight_row = 50,
                         cluster_cols = FALSE,
                         annotation_col = subdf,
                         annotation_colors = ann_colors,
                         show_colnames = FALSE,
                         show_rownames = TRUE,
                         border_color = "grey60",
                         fontsize = 10,
                         main = "Heatmap of Nuclear Receptors")

save_pheatmap_jpeg(p1, "FigS5G.jpeg", width = 10, height = 7)
