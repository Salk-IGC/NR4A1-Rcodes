# FigS8A

# Figure S8. IFN𝛾 reshapes NR4A1 chromatin occupancy in LKB1-deficient lung cancer cells
# (A) CUT&RUN signal coverage heatmaps and histograms of all NR4A1 IFNg-dependent and –independent peaks in sgNT, sgLkb1, and
# sgLkb1+sgNr4a1 cells treated with vehicle or IFN𝛾 that colocalize with NR4A1-binding sites from IFN𝛾–treated sgLkb1 cells.

# bigwigs from GSE302652

library(EnrichedHeatmap)
library(rtracklayer)
library(circlize)
library(RColorBrewer)
library(data.table)
library(here)

# Shared helper functions
source(here("R", "cutandrun_helpers.R"))

# histogram centered at peaks
peakfnm <- c( "control_vehicle.txt", # LKR10 cells, non-targeting control, vehicle treated, Nur77
              "sgLkb1_vehicle.txt", # LKR10 cells, sgLkb1, vehicle treated, Nur77
              "sgLkb1_sg1Nr4a1_vehicle.txt", # LKR10 cells, sgLkb1+sg1.Nr4a1, vehicle treated, Nur77, replicate 1 & 2
              "sgLkb1_sg2Nr4a1_vehicle.txt", # LKR10 cells, sgLkb1+sg2.Nr4a1, vehicle treated, Nur77, replicate 1 & 2
              "control_IFNg.txt", # LKR10 cells, non-targeting control, IFNg treated, Nur77, replicate 1 & 2
              "sgLkb1_IFNg.txt", # LKR10 cells, sgLkb1, IFNg treated, Nur77, replicate 1 & 2
              "sgLkb1_sg1Nr4a1_IFNg.txt", # LKR10 cells, sgLkb1+sg1.Nr4a1, IFNg treated, Nur77, replicate 1 & 2
              "sgLkb1_sg2Nr4a1_IFNg.txt" # LKR10 cells, sgLkb1+sg2.Nr4a1, IFNg treated, Nur77, replicate 1 & 2
)

bigwigfnm <- c(
  "control_vehicle.bw",
  "sgLkb1_vehicle.bw",
  "sgLkb1_sg1Nr4a1_vehicle.bw",
  "sgLkb1_sg2Nr4a1_vehicle.bw",
  "control_IFNg.bw",
  "sgLkb1_IFNg.bw",
  "sgLkb1_sg1Nr4a1_IFNg.bw",
  "sgLkb1_sg2Nr4a1_IFNg.bw"
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

# colors
colors <- colorRampPalette(brewer.pal(9, "OrRd"))(100)
breaks <- seq(0, 0.8, length.out = 100)
col_fun1 <- circlize::colorRamp2(breaks, colors)

colors <- colorRampPalette(brewer.pal(9, "Blues"))(100)
breaks <- seq(0, 0.8, length.out = 100)
col_fun2 <- circlize::colorRamp2(breaks, colors)

# ---- sgLKB1 IFNg vs sgLKB1 vehicle ----

tmp <- fread(peakfnm[6], header = T, data.table = FALSE)
targets <- makeGRangesFromDataFrame(
  df = tmp,
  seqnames.field = "Chr", start.field = "Start", end.field = "End")
targets
ExtendSize <- 3000
targets.extended  <- resize(targets, fix = "center", width = ExtendSize*2)

plotpairheatmap <- function(normMatrix, ifng_i, veh_j) {
  EH_plot <- list()
  #ifng
  EH_plot[[1]] <- EnrichedHeatmap( mat = normMatrix[[1]],
                                   pos_line = FALSE,
                                   border = FALSE,
                                   col = col_fun1,
                                   column_title = plot_title[ifng_i],
                                   column_title_gp = gpar(fontsize = 8, fontfamily = "sans"),
                                   rect_gp = gpar(col = "transparent"),
                                   show_heatmap_legend = FALSE,
                                   top_annotation = HeatmapAnnotation(
                                     enriched = anno_enriched(
                                       gp = gpar(col = "black", lty = 1, lwd=2),
                                       col="black",
                                       axis = FALSE,
                                       ylim = c(0.02, 0.6),
                                       legend_width = unit(3, "cm"),
                                       title_position = "topcenter"))
  )

  #veh
  EH_plot[[2]] <- EnrichedHeatmap( mat =  normMatrix[[2]],
                                   pos_line = FALSE,
                                   border = FALSE,
                                   col = col_fun2,
                                   column_title = plot_title[veh_j],
                                   column_title_gp = gpar(fontsize = 8, fontfamily = "sans"),
                                   rect_gp = gpar(col = "transparent"),
                                   heatmap_legend_param = list(
                                     legend_direction = "vertical",
                                     title = NULL,
                                     legend_height = unit(3, "in")),
                                   show_heatmap_legend = FALSE,
                                   top_annotation = HeatmapAnnotation(
                                     enriched = anno_enriched(
                                       gp = gpar(col = "black", lty = 1, lwd=2),
                                       col="black",
                                       ylim = c(0.02, 0.6),
                                       axis = FALSE,
                                       legend_width = unit(3, "cm"),
                                       title_position = "topcenter")

                                   )
  ) #

  names(EH_plot) <- c(plot_title[ifng_i], plot_title[veh_j])
  return(EH_plot)
}

sgNT_matrix <- countmatrix(targets.extended, 5, 1)
sgAB_matrix <- countmatrix(targets.extended, 7, 3)
sgAD_matrix <- countmatrix(targets.extended, 8, 4)

# colors
colors <- colorRampPalette(brewer.pal(9, "OrRd"))(100)
breaks <- seq(0, 0.3, length.out = 100)
col_fun1 <- circlize::colorRamp2(breaks, colors)

colors <- colorRampPalette(brewer.pal(9, "Blues"))(100)
breaks <- seq(0, 0.3, length.out = 100)
col_fun2 <- circlize::colorRamp2(breaks, colors)

#
sgNT <- plotpairheatmap(sgNT_matrix, 5, 1)
sgAB <- plotpairheatmap(sgAB_matrix, 7, 3)
sgAD <- plotpairheatmap(sgAD_matrix, 8, 4)

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
                        col = col_fun1,
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
                            ylim = c(0.02, 0.6),
                            legend_width = unit(3, "cm"),
                            title_position = "topcenter"))
)

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
                        top_annotation = HeatmapAnnotation(
                          enriched = anno_enriched(
                            gp = gpar(col = "black", lty = 1, lwd=2),
                            col="black",
                            ylim = c(0.02, 0.6),
                            axis = FALSE,
                            legend_width = unit(3, "cm"),
                            title_position = "topcenter")

                        )
) #

# combine all
EH <- EH1 + EH2 + sgNT[[1]] + sgNT[[2]] + sgAB[[1]] + sgAB[[2]] +sgAD[[1]] + sgAD[[2]]
jpeg("FigS8A_topPanel.jpeg", width=15, height=8,  units="in",res=600, type="cairo")
draw(EH,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "right",
     padding = unit(c(4, 4, 4, 4), "mm")
)
dev.off()

# ---- shared and unique peak panels ----
# overlap peaks from sgLKB1 Veh and sgLKB1 IFNg, generated from mergePeaks -d given from peaks_sgLKB1_Veh and peaks_sgLKB1_IFNg.
diffpeaks <- read.delim("peaks_sgLKB1_Veh_IFNg.txt", header=T)
head(diffpeaks)
colnames(diffpeaks)
table(diffpeaks$Source)

ifng_unique_peaks <- diffpeaks[diffpeaks$Source == "peaks_sgLKB1_IFNg.txt",c(1:4)]
shared_peaks <- diffpeaks[diffpeaks$Source == "peaks_sgLKB1_Veh.txt|peaks_sgLKB1_IFNg.txt",c(1:4)]

#--------------
# ifng_unique_peaks: Nr4a1-binding
#-------------

# modified plot function

colors <- colorRampPalette(brewer.pal(9, "OrRd"))(100)
breaks <- seq(0, 0.4, length.out = 100)
col_fun2 <- circlize::colorRamp2(breaks, colors)

# for sgLKB1 Nr4a1 only -unique peaks
targets <- makeGRangesFromDataFrame(
  df = ifng_unique_peaks,
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
sgNT <- plotpairheatmap2(sgNT_matrix, 5, 1, col_fun2, col_fun2, ylim=0.5, legend_direction = "horizontal")
sgAB <- plotpairheatmap2(sgAB_matrix, 7, 3, col_fun2, col_fun2, ylim=0.5, legend_direction = "horizontal")
sgAD <- plotpairheatmap2(sgAD_matrix, 8, 4, col_fun2, col_fun2, ylim=0.5, legend_direction = "horizontal")

# combine all
EH <- EH1 + EH2 + sgNT[[1]] + sgNT[[2]] + sgAB[[1]] + sgAB[[2]] +sgAD[[1]] + sgAD[[2]]
jpeg("FigS8A_unique_peaks.jpeg", width=15, height=8,  units="in",res=600, type="cairo")
draw(EH,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "right",
     padding = unit(c(4, 4, 4, 4), "mm")
)
dev.off()

#---shared peaks in nr4a1.---------------------------------------
targets <- makeGRangesFromDataFrame(
  df = shared_peaks,
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
sgNT <- plotpairheatmap2(sgNT_matrix, 5, 1, col_fun2, col_fun2, ylim=0.65, legend_direction = "horizontal")
sgAB <- plotpairheatmap2(sgAB_matrix, 7, 3, col_fun2, col_fun2, ylim=0.65, legend_direction = "horizontal")
sgAD <- plotpairheatmap2(sgAD_matrix, 8, 4, col_fun2, col_fun2, ylim=0.65, legend_direction = "horizontal")

# combine all
EH <- EH1 + EH2 + sgNT[[1]] + sgNT[[2]] + sgAB[[1]] + sgAB[[2]] +sgAD[[1]] + sgAD[[2]]
jpeg("FigS8A_shared_peaks.jpeg", width=15, height=8,  units="in",res=600, type="cairo")
draw(EH,
     heatmap_legend_side = "bottom",
     annotation_legend_side = "right",
     padding = unit(c(4, 4, 4, 4), "mm")
)
dev.off()
