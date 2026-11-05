# Fig3A,C,G

# Figure 3. NR4A1 (NUR77) is a CRTC2 - CREB target gene that suppresses IFN-dependent transcription.
# (A) Heatmap displaying the top 25 CRTC2-dependent, differentially expressed (DE) genes (FC ³ 2, FDR ≤ 0.05) in sgNT, sgLkb1, and sgLkb1+sgCrtc2 LKR10 cells.
# (C) Heatmap depicting mRNA expression of LKB1-dependent IFNg-response genes in sgNT, sgLkb1, and sgLkb1+sgNr4a1 LKR10 cells.

# R barrel 4.3.2
library(DESeq2)
library(pheatmap)
library(ggplot2)
library(cowplot)
library(dplyr)
library(here)

# Shared helper functions
source(here("R", "plot_helpers.R"))

# ---- Figure 3A ----
# ---- Load QC object (run 1) ----
# [[1]] = DESeq2 object, [[2]] = transformed expression (used for heatmaps)
mouse1 <- readRDS(here("results", "QC", "run1_QC.RDS"))
subdat <- mouse1[[2]]

degsub <- c("Anxa10", "Slc17a6", "Tmem252", "Tmem179", "Gzmd", "Tff1", "Fabp6", "Cd93",
         "Ndnf", "Ppbp", "Gkn1", "Krt20", "Has2", "Prl2c3", "Galnt9", "Prl2c4", "Clca3b",
         "Ramp3", "Prl2c2", "Mmp10", "Plat", "Gm19461", "Nr4a1", "Ivl", "Fam25c")

# ---- Sample annotation: drop sg1.CREB1, order vehicle then IFNg ----
subdf <- as.data.frame(colData(subdat)[, c("guide_KO", "treatment")])
names(subdf) <- c("KO", "Treatment")
subdf <- subdf[!subdf$KO %in% c("sg1.CREB1", "sg2.CREB1"), ]
subdf$Treatment <- factor(subdf$Treatment, levels = c("v", "F"))   # vehicle, IFNg
subdf$KO <- factor(subdf$KO, levels = c("sg1.NT", "sg1.LKB1",
                                        "sg1.CRTC2", "sg2.CRTC2"))
subdf <- subdf[order(subdf$Treatment, subdf$KO), ]
subdat <- assay(subdat)[rownames(assay(subdat)) %in% degsub, ]
subdat <- subdat[, -c(1:6, 13:18)]   # drops the first 6 columns ( sg1.CREB1 samples)

zero_var_rows <- apply(subdat, 1, var) == 0
message("Zero-variance genes removed: ", sum(zero_var_rows))
subdat <- subdat[!zero_var_rows, ]
message("Genes plotted: ", nrow(subdat))

subdat <- subdat[, rownames(subdf)]
subdat <- subdat[degsub,]

lin <- counts(mouse1[[1]], normalized = TRUE)[rownames(subdat), colnames(subdat)]
mm <- t(apply(lin, 1, function(x) (x - min(x)) / (max(x) - min(x))))
subdf2 <- subdf[, c("Treatment", "KO")]
names(subdf2) <- c("Treatment","Cell Line")

ann_colors2 <- list(
  Treatment   = c("v" = "turquoise", "F" = "gold2"),
  `Cell Line` = c("sg1.NT"    = "grey",
                  "sg1.LKB1"  = "darkred",
                  "sg2.CREB1" = "coral",
                  "sg1.CRTC2" = "#00B3FF",
                  "sg2.CRTC2" = "#00B3FF")
)

p1 <- pheatmap::pheatmap(mm,
                         scale             = "none",
                         breaks            = seq(0, 1, length.out = 101),
                         color             = colorRampPalette(c("#3F51A3", "white", "#EE2B2B"))(100),
                         legend_breaks     = c(0, 1),
                         legend_labels     = c("row min", "row max"),
                         cluster_rows      = FALSE,
                         cluster_cols      = FALSE,
                         show_rownames     = TRUE,
                         show_colnames     = FALSE,
                         border_color      = NA,
                         main              = NA,
                         annotation_col    = subdf2,
                         annotation_colors = ann_colors2,
                         annotation_legend = FALSE,
                         fontsize          = 10)

save_pheatmap_jpeg(p1, "Fig3A.jpeg", width = 5, height = 4)

# ---- Figure 3C ----
mouse2 <- readRDS(here("results", "QC", "run2_QC.RDS"))
# ---- Load LKB1-dependent gene list (n = 166, Table S2) ----
lkb1_dependent <- read.delim(here("data", "tableS2.txt"),
                             header = TRUE)$gene_name

subdat <- mouse2[[2]]
subdf <- as.data.frame(colData(subdat)[,  c("guide_KO", "treatment")])
names(subdf) <- c("KO","Treatment")
subdat <- assay(subdat)[rownames(assay(subdat)) %in% lkb1_dependent, ]
subdf$Treatment <- factor(subdf$Treatment, levels=c("v", "F"))
subdf$KO <- factor(subdf$KO, levels=c("sg1.NT1", "sg1.LKB1", "sg1.Nr4a1", "sg2.Nr4a1"))
subdat <- subdat[,rownames(subdf[order(subdf$Treatment, subdf$KO),])]
zero_var_rows <- apply(subdat, 1, var) == 0
subdat <- subdat[!zero_var_rows, ]

ann <- data.frame(
  `Cell Line` = factor(subdf$KO,
                       levels = c("sg1.NT1", "sg1.LKB1", "sg1.Nr4a1", "sg2.Nr4a1")),
  Treatment   = factor(ifelse(subdf$Treatment == "v", "Vehicle", "IFNγ"),
                       levels = c("Vehicle", "IFNγ")),
  check.names = FALSE,
  row.names   = rownames(subdf)
)

ord    <- order(ann$Treatment, ann$`Cell Line`)
subdat <- subdat[, rownames(ann)[ord]]
ann    <- ann[ord, , drop = FALSE]
trt_cols <- c("Vehicle" = "#8ED6D3", "IFNγ" = "#FFEE58")
cl_cols  <- c("sg1.NT1" = "#BDBDBD", "sg1.LKB1" = "#B5403F",
              "sg1.Nr4a1" = "#8A6BB1", "sg2.Nr4a1" = "#3DB54A")
ann_colors <- list(Treatment = trt_cols, `Cell Line` = cl_cols)
hm_cols   <- colorRampPalette(c("#3F51A3", "white", "#EE2B2B"))(100)
hm_breaks <- seq(-2.5, 2.5, length.out = 101)

p1 <- pheatmap::pheatmap(
  subdat,
  scale = "row",
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  treeheight_row = 0,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = ann,
  annotation_colors = ann_colors,
  annotation_legend = T,
  annotation_names_col = TRUE,
  border_color = NA,
  main = NA,
  color = hm_cols,
  breaks = hm_breaks,
  legend = FALSE,
  fontsize = 8,
  silent = TRUE
)

d <- data.frame(x = 1, y = 1, v = 0)
leg_plot <- ggplot(d, aes(x, y, fill = v)) + geom_tile() +
  scale_fill_gradientn(
    "Expression",
    colors = hm_cols, limits = c(-2.5, 2.5), breaks = -2:2,
    guide = guide_colorbar(title.position = "top", title.hjust = 0,
                           barwidth = 1, barheight = 8,
                           frame.colour = NA, ticks = FALSE)
  ) +
  theme(legend.title = element_text(size = 10),
        legend.text  = element_text(size = 10))
leg_cb <- cowplot::get_legend(leg_plot)

final <- ggdraw() +
  draw_grob(p1$gtable, x = 0.08, y = 0, width = 0.78, height = 1) +
  draw_plot(leg_cb, x = 0.88, y = 0.5, width = 0.12, height = 0.4) +
  draw_label("LKB1-dependent\nIFNγ Response Genes (n=166)",
             x = 0.03, y = 0.4, angle = 90, size = 10) +
  draw_line(x = c(0.065, 0.065), y = c(0.03, 0.8), linewidth = 0.6)
ggsave("Fig3C.jpeg", final, width = 6, height = 5.5, units = "in", dpi = 300)

# ---- Figure 3G ----
# (G) Heatmap depicting mRNA expression of LKB1-dependent IFNg-response genes in EV and NR4A1-WT overexpression (o.e.) LKR10 cells

set.seed(12345678)

lkb1_dependent <- read.delim(here("data", "tableS2.txt"), header=T)
lkb1_dependent <- lkb1_dependent$gene_name

# Generate dds from the GSE351227 raw counts
count_file <- "rsem.merged.gene_counts.tsv"
count.table <- read.delim(count_file, header = TRUE, row.names = NULL, check.names = FALSE)

is_dup <- duplicated(count.table$gene_name) |
  duplicated(count.table$gene_name, fromLast = TRUE)
rownames(count.table) <- ifelse(is_dup, count.table$gene_id, count.table$gene_name)

counts <- count.table[, 3:ncol(count.table)]

parts    <- strsplit(colnames(counts), "_")
metadata <- sapply(parts, `[`, 1)
vars     <- do.call(rbind, strsplit(metadata, ".", fixed = TRUE))

coldata <- data.frame(
  Sample     = colnames(counts),
  Variable1  = vars[, 1],
  Variable2  = vars[, 2],
  Metadata   = metadata,
  Replicate  = sapply(parts, `[`, 2),
  row.names  = colnames(counts)
)

dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData   = coldata,
  design    = ~Metadata
)
dds <- DESeq(dds)

assay(dds, "vst") <- assay(varianceStabilizingTransformation(dds))

# generate figure
subdat <- dds
subdf <- as.data.frame(colData(subdat)[,  c("Variable1", "Variable2")])
names(subdf) <- c("condition","Treatment")
subdat <- assay(subdat, "vst")[rownames(assay(subdat, "vst")) %in% lkb1_dependent, ]
dim(subdat)
subdf$Treatment <- factor(subdf$Treatment, levels=c("veh", "IFN"))
subdf$condition <- factor(subdf$condition, levels=c("EV", "WT" ))
subdat <- subdat[,rownames(subdf[order(subdf$Treatment, subdf$condition),])]
colnames(subdat)
zero_var_rows <- apply(subdat, 1, var) == 0
sum(zero_var_rows)  # Count zero variance rows
subdat <- subdat[!zero_var_rows, ]
colnames(subdat)
dim(subdat)

ann <- data.frame(
  `Cell Line` = factor(subdf$condition, levels=c("EV", "WT")),
  Treatment   = factor(ifelse(subdf$Treatment == "veh", "Vehicle", "IFNγ"),
                       levels = c("Vehicle", "IFNγ")),
  check.names = FALSE,
  row.names   = rownames(subdf)
)

ord    <- order(ann$Treatment, ann$`Cell Line`)
subdat <- subdat[, rownames(ann)[ord]]
ann    <- ann[ord, , drop = FALSE]

trt_cols <- c("Vehicle" = "#8ED6D3", "IFNγ" = "#FFEE58")
cl_cols  <- c("EV" = "#BDBDBD", "WT" = "green4")
ann_colors <- list(Treatment = trt_cols, `Cell Line` = cl_cols)
hm_cols   <- colorRampPalette(c("#3F51A3", "white", "#EE2B2B"))(100)
hm_breaks <- seq(-2.5, 2.5, length.out = 101)

# Heatmap
p1 <- pheatmap::pheatmap(
  subdat,
  scale = "row",
  cluster_rows = TRUE,
  clustering_distance_rows = "euclidean",
  clustering_method = "ward.D2",
  treeheight_row = 0,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = ann,
  annotation_colors = ann_colors,
  annotation_legend = T,
  annotation_names_col = TRUE,
  border_color = NA,
  main = NA,
  color = hm_cols,
  breaks = hm_breaks,
  legend = F,                 # color bar is drawn separately below
  fontsize = 8,
  silent = TRUE
)

d <- data.frame(x = 1, y = 1, v = 0)
leg_plot <- ggplot(d, aes(x, y, fill = v)) + geom_tile() +
  scale_fill_gradientn(
    "Expression",
    colors = hm_cols, limits = c(-2.5, 2.5), breaks = -2:2,
    guide = guide_colorbar(title.position = "top", title.hjust = 0,
                           barwidth = 1, barheight = 8,
                           frame.colour = NA, ticks = FALSE)
  ) +
  theme(legend.title = element_text(size = 10),
        legend.text  = element_text(size = 10))
leg_cb <- cowplot::get_legend(leg_plot)

final <- ggdraw() +
  draw_grob(p1$gtable, x = 0.08, y = 0, width = 0.78, height = 1) +
  draw_plot(leg_cb, x = 0.88, y = 0.5, width = 0.12, height = 0.4) +
  draw_label("LKB1-dependent\nIFNγ Response Genes",
             x = 0.03, y = 0.4, angle = 90, size = 10) +
  draw_line(x = c(0.065, 0.065), y = c(0.03, 0.8), linewidth = 0.6)
ggsave("Fig3G.pdf", final, width = 6, height = 5.5, units = "in")
