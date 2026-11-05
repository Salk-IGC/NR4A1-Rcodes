# FIG4: A,B,C,D
# Figure 4. NR4A1 (NUR77) is a CRTC2 - CREB target gene that suppresses IFN-dependent transcription.
# (A) NR4A1 CUT&RUN signal coverage histograms of all NR4A1-bound sites in IFN𝛾-treated sgLkb1 cells across sgNT, sgLkb1, and
# sgLkb1+sgNr4a1 cells treated with IFN𝛾 or vehicle. (B) NR4A1 CUT&RUN signal coverage histograms of promoter-associated Nr4a1-bound sites in
# IFN𝛾-treated sgLkb1 cells across sgNT, sgLkb1, and sgLkb1+sgNr4a1 cells treated with IFN𝛾 or vehicle. (C) NR4A1 CUT&RUN signal coverage
# heatmaps and histograms of IFN𝛾-dependent and –independent promoter-associated peaks in sgNT, sgLkb1, and sgLkb1+sgNr4a1 cells treated
# with IFN𝛾 or vehicle. (D) Heatmap of gene expression for NR4A1-bound DEGs at promoter regions from IFN𝛾-treated sgLkb1 cells across sgNT,
# sgLkb1, and sgLkb1+sgNr4a1 cells treated with IFN𝛾 or vehicle. Black arrows indicate Stat1, Stat2, and IRF family genes

# CUT&RUN data from GSE302652

library(EnrichedHeatmap)
library(rtracklayer)
library(circlize)
library(data.table)
library(reshape2)
library(ggplot2)
library(grid)
library(gtable)
library(RColorBrewer)
library(DESeq2)
library(pheatmap)
library(here)

# Shared helper functions
source(here("R", "cutandrun_helpers.R"))

tss = fread("mm10.tss", header = FALSE, data.table = FALSE)

peakfnm <- c("peaks_sgLkb1_IFNg_anno.txt", # peaks from LKR10 cells, sgLkb1, IFNg treated, Nur77, replicate 1&2
             "peaks_sgLkb1_sg1Nr4a1_IFNg_anno.txt", # peaks from LKR10 cells, sgLkb1+sg1.Nr4a1, IFNg treated, Nur77, replicate 1&2
             "peaks_sgLkb1_sg2Nr4a1_IFNg_anno.txt") # LKR10 cells, sgLkb1+sg2.Nr4a1, IFNg treated, Nur77, replicate 1&2

for (i in 1:length(peakfnm)){
  peakanno = fread(peakfnm[i], header = TRUE, data.table = FALSE)
  peakanno <- unique(peakanno$`Nearest PromoterID`)
  df <- tss[tss$V1 %in% peakanno,]
  head(df)
  write.table(df[, -1], paste0("hist/TSS_", strsplit(peakfnm[i], "./peaks_|_anno.txt")[[1]][2], "_anno.txt"),
              sep="\t", col.names = F, row.names = T, quote=F)
}

# ---- Figure 4A ----
# generate the histogram for peaks using HOMER annotatePeaks.pl, extract raw counts
dat <- read.delim("peaks_sgLkb1_IFNg_anno_allSample_hist.txt", header=T)
names(dat)
head(dat)
dat <- dat[, c(1,2,5,8,11,14,17,20,23)]
names(dat) <- c("Distance from peak center",
                "sgNT Veh", "sgLKB1 Veh", "sgAB Veh", "sgAD Veh",
                "sgNT IFNg",
                "sgLKB1 IFNg",
                "sgAB IFNg",
                "sgAD IFNg")
dat_re <- melt(dat, id="Distance from peak center")
names(dat_re) <- c("Distance", "Condition", "Coverage")
head(dat_re)

dat_re$Condition <- factor(dat_re$Condition, levels=c( "sgNT Veh","sgNT IFNg",
                                                       "sgLKB1 Veh","sgLKB1 IFNg",
                                                      "sgAB Veh","sgAB IFNg",
                                                      "sgAD Veh", "sgAD IFNg"))

gg <- ggplot(dat_re , aes(x=Distance,y=Coverage,  color=Condition)) +
  geom_line(alpha = 1, size = 1.2)+
  labs(title="All peaks in sgLKB1 IFNg",
       x="Distance from peak center",
       y = "Coverage (per bp per peak)") +

  scale_color_manual(values = c("sgNT Veh"="black","sgNT IFNg"= "#D9D9D9",
                                "sgLKB1 Veh"="red","sgLKB1 IFNg"="#FF9999",
                                "sgAB Veh"="#CC00FF", "sgAB IFNg"="#F0B3FF",
                                "sgAD Veh"="#00C49A", "sgAD IFNg"= "#00FF66"))+

  scale_x_continuous(   breaks = c(-3000, 0, 3000),
                        labels = c("-3k", "0", "3k")) +
  geom_vline(xintercept = 0, linetype="dashed", color = "grey", size=0.5) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    legend.title = element_blank(),
    legend.position = c(0.75, 0.5),
    legend.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )
ggsave("Fig4A.jpeg", width=6, height=4,  units = "in", dpi=600)

# ---- Figure 4B ----
dat <- read.delim("peaks_sgLkb1_sg1Nr4a1_IFNg_anno_allSample_hist.txt", header=T)
names(dat)
head(dat)
dat <- dat[, c(1,2,5,8,11,14,17,20,23)]
names(dat) <- c("Distance from TSS",
                "sgNT Veh", "sgLKB1 Veh", "sgAB Veh", "sgAD Veh",
                "sgNT IFNg",
                "sgLKB1 IFNg",
                "sgAB IFNg",
                "sgAD IFNg")
dat_re <- melt(dat, id="Distance from TSS")
names(dat_re) <- c("Distance", "Condition", "Coverage")
head(dat_re)

dat_re$Condition <- factor(dat_re$Condition, levels=c( "sgNT Veh","sgNT IFNg",
                                                       "sgLKB1 Veh","sgLKB1 IFNg",
                                                       "sgAB Veh","sgAB IFNg",
                                                       "sgAD Veh", "sgAD IFNg"))

gg <- ggplot(dat_re , aes(x=Distance,y=Coverage,  color=Condition)) +
  geom_line(alpha=1, size=1.2)+
  labs(title="sgLKB1 IFNg nearby genes",
       x="Distance from TSS",
       y = "Coverage (per bp per peak)") +

  scale_color_manual(values = c("sgNT Veh"="black","sgNT IFNg"= "#D9D9D9",
                                "sgLKB1 Veh"="red","sgLKB1 IFNg"="#FF9999",
                                "sgAB Veh"="#CC00FF", "sgAB IFNg"="#F0B3FF",
                                "sgAD Veh"="#00C49A", "sgAD IFNg"= "#00FF66"))+

  scale_x_continuous(   breaks = c(-3000, 0, 3000),
                        labels = c("-3k", "0", "3k")) +
  geom_vline(xintercept = 0, linetype="dashed", color = "grey", size=0.5) +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    legend.title = element_blank(),
    legend.position = c(0.75, 0.5),
    legend.background = element_rect(fill = "white", color = NA),
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )
ggsave("Fig4B.jpeg", width=6, height=4,  units = "in", dpi=600)

# ---- Figure 4C ----
# bigwigs from GSE302652
bigwigfnm <- c(
  "control_vehicle.bw", # LKR10 cells, non-targeting control, vehicle treated, Nur77
  "sgLkb1_vehicle.bw", # LKR10 cells, sgLkb1, vehicle treated, Nur77
  "sgLkb1_sg1Nr4a1_vehicle.bw", # LKR10 cells, sgLkb1+sg1.Nr4a1, vehicle treated, Nur77, replicate 1 & 2
  "sgLkb1_sg2Nr4a1_vehicle.bw", # LKR10 cells, sgLkb1+sg2.Nr4a1, vehicle treated, Nur77, replicate 1 & 2
  "control_IFNg.bw", # LKR10 cells, non-targeting control, IFNg treated, Nur77, replicate 1 & 2
  "sgLkb1_IFNg.bw", # LKR10 cells, sgLkb1, IFNg treated, Nur77, replicate 1 & 2
  "sgLkb1_sg1Nr4a1_IFNg.bw", # LKR10 cells, sgLkb1+sg1.Nr4a1, IFNg treated, Nur77, replicate 1 & 2
  "sgLkb1_sg2Nr4a1_IFNg.bw" # LKR10 cells, sgLkb1+sg2.Nr4a1, IFNg treated, Nur77, replicate 1 & 2
)

plot_title <- c("sgNT Veh",
                "sgLKB1 Veh",
                "sgAB Veh",
                "sgAD Veh",
                "sgNT IFNg",
                "sgLKB1 IFNg",
                "sgAB IFNg",
                "sgAD IFNg")

plot_name <- paste0(plot_title, ".jpeg")

# overlap peaks from sgLKB1 Veh and sgLKB1 IFNg, generated from mergePeaks -d given from peaks_sgLKB1_Veh and peaks_sgLKB1_IFNg.
diffpeaks <- read.delim("peaks_sgLKB1_Veh_IFNg.txt", header=T)
head(diffpeaks)
colnames(diffpeaks)
table(diffpeaks$Source)
dim(diffpeaks) #43272

diffpeaks$simple <- sapply(as.vector(diffpeaks$Annotation),function(x) strsplit(x," ",fixed=T)[[1]][1])
table(diffpeaks$simple)
diffpeaks <- diffpeaks[diffpeaks$simple=="promoter-TSS",]
dim(diffpeaks) #7527

ifng_unique_peaks <- diffpeaks[diffpeaks$Source == "peaks_sgLKB1_IFNg.txt",c(1:4)]
shared_peaks <- diffpeaks[diffpeaks$Source == "peaks_sgLKB1_Veh.txt|peaks_sgLKB1_IFNg.txt",c(1:4)]

dim(ifng_unique_peaks) #4099 , remove 1 NA
dim(shared_peaks) #2939

#--------------
# ifng_unique_peaks: Nr4a1-binding
#-------------

# colors - Nr4a1-binding
colors <- colorRampPalette(brewer.pal(9, "OrRd"))(100)
breaks <- seq(0, 0.4, length.out = 100)
col_fun2 <- circlize::colorRamp2(breaks, colors)

# for sgLKB1 Nr4a1 only -unique peaks
targets <- makeGRangesFromDataFrame(
  df = ifng_unique_peaks,na.rm=TRUE,
  seqnames.field = "Chr", start.field = "Start", end.field = "End")
targets
ExtendSize <- 3000
targets.extended  <- resize(targets, fix = "center", width = ExtendSize*2)

BigWig1 <- rtracklayer::import(bigwigfnm[6],
                               format = "BigWig",
                               selection = BigWigSelection(targets.extended))
normMatrix1 <- normalizeToMatrix(signal = BigWig1,
                                 target = resize(targets, fix = "center", width = 1),
                                 background = 0,
                                 keep = c(0, 0.99),
                                 target_ratio = 0,
                                 mean_mode = "w0",
                                 value_column = "score",
                                 extend = ExtendSize)

#add another bw
BigWig <- rtracklayer::import(bigwigfnm[2],
                              format = "BigWig",
                              selection = BigWigSelection(targets.extended))
normMatrix <- normalizeToMatrix(signal = BigWig,
                                target = resize(targets, fix = "center", width = 1),
                                background = 0,
                                keep = c(0, 0.99),
                                target_ratio = 0,
                                mean_mode = "w0",
                                value_column = "score",
                                extend = ExtendSize)

EH1 <- EnrichedHeatmap( mat = normMatrix1,
                        pos_line = FALSE,
                        border = FALSE,
                        col = col_fun2,
                        column_title = plot_title[6],
                        column_title_gp = gpar(fontsize = 8, fontfamily = "sans"),
                        rect_gp = gpar(col = "transparent"),
                        heatmap_legend_param = list(
                          legend_direction = "horizontal",
                          title = NULL,
                          legend_height = unit(3, "in")),
                        top_annotation = HeatmapAnnotation(
                          enriched = anno_enriched(
                            gp = gpar(col = "black", lty = 1, lwd=2),
                            col="black",
                            axis_param = list(side = "left"),
                            ylim = c(0.02, 0.5),
                            legend_width = unit(3, "cm"),
                            title_position = "topcenter")))

EH2 <- EnrichedHeatmap( mat = normMatrix,
                        pos_line = FALSE,
                        border = FALSE,
                        col = col_fun2,
                        column_title = plot_title[2],
                        column_title_gp = gpar(fontsize = 8, fontfamily = "sans"),
                        rect_gp = gpar(col = "transparent"),
                        heatmap_legend_param = list(
                          legend_direction = "horizontal",
                          title = NULL,
                          legend_height = unit(3, "in")),
                        show_heatmap_legend = FALSE,
                        top_annotation = HeatmapAnnotation(
                          enriched = anno_enriched(
                            gp = gpar(col = "black", lty = 1, lwd=2),
                            col="black",
                            ylim = c(0.02, 0.5),
                            axis = FALSE,
                            legend_width = unit(3, "cm"),
                            title_position = "topcenter")

                        ))

#add other samples
sgNT_matrix <- countmatrix(targets.extended, 5, 1)
sgAB_matrix <- countmatrix(targets.extended, 7, 3)
sgAD_matrix <- countmatrix(targets.extended, 8, 4)
#
sgNT <- plotpairheatmap2(sgNT_matrix, 5, 1, col_fun2, col_fun2, ylim=0.5)
sgAB <- plotpairheatmap2(sgAB_matrix, 7, 3, col_fun2, col_fun2, ylim=0.5)
sgAD <- plotpairheatmap2(sgAD_matrix, 8, 4, col_fun2, col_fun2, ylim=0.5)

# combine all
EH <- EH1 + EH2 + sgNT[[1]] + sgNT[[2]] + sgAB[[1]] + sgAB[[2]] +sgAD[[1]] + sgAD[[2]]
jpeg("Fig4C_IFNg_unique_peaks.jpeg", width=15, height=8,  units="in",res=600, type="cairo")
ComplexHeatmap::draw(EH,
                     heatmap_legend_side = "bottom",
                     annotation_legend_side = "right",
                     padding = unit(c(4, 4, 4, 4), "mm")
)
dev.off()

#---shared peaks in nr4a1.---------------------------------------
targets <- makeGRangesFromDataFrame(
  df = shared_peaks, na.rm=T,
  seqnames.field = "Chr", start.field = "Start", end.field = "End")
targets

ExtendSize <- 3000
targets.extended  <- resize(targets, fix = "center", width = ExtendSize*2)

BigWig1 <- rtracklayer::import(bigwigfnm[6],
                               format = "BigWig",
                               selection = BigWigSelection(targets.extended))
normMatrix1 <- normalizeToMatrix(signal = BigWig1,
                                 target = resize(targets, fix = "center", width = 1),
                                 background = 0,
                                 keep = c(0, 0.99),
                                 target_ratio = 0,
                                 mean_mode = "w0",
                                 value_column = "score",
                                 extend = ExtendSize)

#add another bw
BigWig <- rtracklayer::import(bigwigfnm[2],
                              format = "BigWig",
                              selection = BigWigSelection(targets.extended))
normMatrix <- normalizeToMatrix(signal = BigWig,
                                target = resize(targets, fix = "center", width = 1),
                                background = 0,
                                keep = c(0, 0.99),
                                target_ratio = 0,
                                mean_mode = "w0",
                                value_column = "score",
                                extend = ExtendSize)

EH1 <- EnrichedHeatmap( mat = normMatrix1,
                        pos_line = FALSE,
                        border = FALSE,
                        col = col_fun2,
                        column_title = plot_title[6],
                        column_title_gp = gpar(fontsize = 8, fontfamily = "sans"),
                        rect_gp = gpar(col = "transparent"),
                        heatmap_legend_param = list(
                          legend_direction = "horizontal",
                          title = NULL,
                          legend_height = unit(3, "in")),
                        top_annotation = HeatmapAnnotation(
                          enriched = anno_enriched(
                            gp = gpar(col = "black", lty = 1, lwd=2),
                            col="black",
                            ylim = c(0.02, 0.65),
                            axis_param = list(side = "left"),
                            legend_width = unit(3, "cm"),
                            title_position = "topcenter")

                        ))

EH2 <- EnrichedHeatmap( mat = normMatrix,
                        pos_line = FALSE,
                        border = FALSE,
                        col = col_fun2,
                        column_title = plot_title[2],
                        column_title_gp = gpar(fontsize = 8, fontfamily = "sans"),
                        rect_gp = gpar(col = "transparent"),
                        heatmap_legend_param = list(
                          legend_direction = "horizontal",
                          title = NULL,
                          legend_height = unit(3, "in")),
                        show_heatmap_legend = FALSE,
                        top_annotation = HeatmapAnnotation(
                          enriched = anno_enriched(
                            gp = gpar(col = "black", lty = 1, lwd=2),
                            col="black",
                            ylim = c(0.02, 0.65),
                            axis = FALSE,
                            legend_width = unit(3, "cm"),
                            title_position = "topcenter")

                        ))

#add other samples
sgNT_matrix <- countmatrix(targets.extended, 5, 1)
sgAB_matrix <- countmatrix(targets.extended, 7, 3)
sgAD_matrix <- countmatrix(targets.extended, 8, 4)
#
sgNT <- plotpairheatmap2(sgNT_matrix, 5, 1, col_fun2, col_fun2, ylim=0.65)
sgAB <- plotpairheatmap2(sgAB_matrix, 7, 3, col_fun2, col_fun2, ylim=0.65)
sgAD <- plotpairheatmap2(sgAD_matrix, 8, 4, col_fun2, col_fun2, ylim=0.65)

# combine all
EH <- EH1 + EH2 + sgNT[[1]] + sgNT[[2]] + sgAB[[1]] + sgAB[[2]] +sgAD[[1]] + sgAD[[2]]
jpeg("Fig4C_IFNg_shared_peaks.jpeg", width=15, height=8,  units="in",res=600, type="cairo")
ComplexHeatmap::draw(EH,
                     heatmap_legend_side = "bottom",
                     annotation_legend_side = "right",
                     padding = unit(c(4, 4, 4, 4), "mm")
)
dev.off()

# ---- Figure 4D ----

# 254 DEGs that are promoter-bound by NR4A1; Table S3
deg <- read.delim("tableS3.txt", header=T)
deg <- deg$Gene.List

mouse2 <- readRDS("run2_QC.RDS")
mouse2
colData(mouse2[[2]])
table(mouse2[[2]]$guide_KO, mouse2[[2]]$treatment)

#heatmaps of 254 DEGs that are promoter-bound by NR4A1 on RNAseq LKB1 and NT
dat <- assay(mouse2[[2]])[rownames(assay(mouse2[[2]])) %in% deg, ]
dim(dat) #254
head(dat)

df <- as.data.frame(colData(mouse2[[2]])[, c("guide_KO", "treatment")])
names(df) <- c("KO","Treatment")
head(df)
str(df)
rownames(df) %in% colnames(dat)
df$Treatment <- factor(df$Treatment, levels=c("v", "F"))
df$KO <- factor(df$KO, levels=c("sg1.NT1", "sg1.LKB1", "sg1.Nr4a1", "sg2.Nr4a1"))
dim(df)

dat <- dat[,rownames(df[order(df$Treatment,df$KO),])]
unique(df$KO)

colnames(dat)
dim(dat)

ann_colors <- list(
  KO = c("sg1.NT1" = "grey", "sg1.LKB1" = "#FF004D", "sg1.Nr4a1"="#CC00FF", "sg2.Nr4a1"="#00FF66"),
  Treatment = c( "v" = "turquoise","F" = "gold2"),
  cluster = c("1" = "#E41A1C",
              "2" = "#377EB8",
              "3" = "#4DAF4A",
              "4" = "#984EA3",
              "5" = "#FF7F00")
)

p1 <- pheatmap::pheatmap(dat,
                         cluster_rows = TRUE,
                         show_rownames = F,
                         cluster_cols = F,
                         annotation_col = df,
                         annotation_colors = ann_colors,
                         main = "Nr4a1-bind-promotor genes from IFNg-sgLKB1; All Sample\n(n=6320)",
                         scale = "row",
                         clustering_distance_rows = "euclidean",
                         clustering_method = "ward.D2",
                         treeheight_row = 25,
                         cutree_rows = 5,
                         cutree_cols = 4,
                         fontsize = 10,
                         color = colorRampPalette(c("blue", "white", "red"))(50))

my_gene_col <- sort(cutree(p1$tree_row, k=5))
my_gene_df <- data.frame(row.names=names(my_gene_col), cluster = as.factor(my_gene_col))

ann_colors <- list(
  KO = c("sg1.NT1" = "grey", "sg1.LKB1" = "#B5403F", "sg1.Nr4a1" = "#8A6BB1", "sg2.Nr4a1" = "#3DB54A"),
  Treatment = c("v" = "turquoise", "F" = "gold2"),
  cluster = c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A",
              "4" = "#984EA3", "5" = "#FF7F00")
)

k <- 5
cl_raw <- cutree(p1$tree_row, k = k)

## renumber clusters 1..k in top-to-bottom display order
disp_order <- p1$tree_row$order
my_gene_col <- match(cl_raw, unique(cl_raw[disp_order]))
names(my_gene_col) <- names(cl_raw)

#show gene names
my_gene_df <- data.frame(row.names = names(my_gene_col),
                         cluster = factor(my_gene_col, levels = 1:k))

p2 <- pheatmap::pheatmap(dat,
                         cluster_rows = TRUE,
                         show_rownames = TRUE,
                         fontsize_row  = 4,
                         cluster_cols = FALSE,
                         annotation_col = df,
                         annotation_colors = ann_colors,
                         annotation_row = my_gene_df,
                         annotation_names_row = FALSE,
                         main = "DEGs with NR4A1 bound at promoter (n=254)",
                         scale = "row",
                         clustering_distance_rows = "euclidean",
                         clustering_method = "ward.D2",
                         treeheight_row = 25,
                         cutree_rows = k,
                         fontsize = 10,
                         color = colorRampPalette(c("blue", "white", "red"))(50),
                         silent = TRUE)

## blocks in display order (now simply 1, 2, 3, 4, 5 from top to bottom)
cl_disp <- my_gene_col[p2$tree_row$order]
n       <- length(cl_disp)
r       <- rle(as.integer(cl_disp))
ends    <- cumsum(r$lengths)
starts  <- ends - r$lengths + 1

coord <- pheatmap:::find_coordinates(n, gaps = ends[-length(ends)])
ypos  <- unit(1, "npc") - (coord$coord - 0.5 * coord$size)

gt <- p2$gtable
ai <- which(gt$layout$name == "row_annotation")
gt$widths[gt$layout$l[ai]] <- unit(16, "pt")

for (b in seq_along(r$values)) {
  mid <- 0.5 * (ypos[starts[b]] + ypos[ends[b]])
  gt <- gtable_add_grob(gt,
                        textGrob(paste("Cluster", r$values[b]), x = 0.5, y = mid, rot = 90,
                                 gp = gpar(col = "white", fontface = "bold", fontsize = 8)),
                        t = gt$layout$t[ai], l = gt$layout$l[ai],
                        clip = "off", name = paste0("cluster_label_", b))
}

jpeg("Fig4D.jpeg", width = 6, height = 15, units = "in", res = 600, type="cairo")
grid.newpage(); grid.draw(gt)
dev.off()

## CSV now uses the same numbering as the figure
my_gene_df$gene_name <- rownames(my_gene_df)
write.csv(my_gene_df, "Fig4D_gene_list.csv", quote = FALSE, row.names = FALSE)
table(my_gene_df$cluster)
