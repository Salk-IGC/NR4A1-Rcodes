# CUT&RUN helper functions shared by Fig4 and FigS8
#
# Required packages: EnrichedHeatmap (loads ComplexHeatmap), rtracklayer, GenomicRanges
# The functions read these objects from the calling script:
#   countmatrix():      bigwigfnm, plot_title, targets, ExtendSize
#   plotpairheatmap2(): plot_title

# countmatrix(targets.extended, ifng_i, veh_j)
# Imports the IFNg bigwig (index ifng_i) and the vehicle bigwig (index veh_j) over the
# extended target peaks and returns the two normalizeToMatrix() coverage matrices,
# named by plot_title.
countmatrix <- function(targets.extended, ifng_i, veh_j) {
  normMatrix <- list()
  #ifng
  BigWig1 <- rtracklayer::import(bigwigfnm[ifng_i],
                                 format = "BigWig",
                                 selection = BigWigSelection(targets.extended))
  normMatrix[[1]] <- normalizeToMatrix(signal = BigWig1,
                                       target = resize(targets, fix = "center", width = 1),
                                       background = 0,
                                       keep = c(0, 0.99),
                                       target_ratio = 0,
                                       mean_mode = "w0",
                                       value_column = "score",
                                       extend = ExtendSize)

  #veh
  BigWig <- rtracklayer::import(bigwigfnm[veh_j],
                                format = "BigWig",
                                selection = BigWigSelection(targets.extended))
  normMatrix[[2]]  <- normalizeToMatrix(signal = BigWig,
                                        target = resize(targets, fix = "center", width = 1),
                                        background = 0,
                                        keep = c(0, 0.99),
                                        target_ratio = 0,
                                        mean_mode = "w0",
                                        value_column = "score",
                                        extend = ExtendSize)

  names(normMatrix) <- c(plot_title[ifng_i], plot_title[veh_j])
  return(normMatrix)
}

# plotpairheatmap2(normMatrix, ifng_i, veh_j, color1, color2, ylim, legend_direction)
# Builds the IFNg and vehicle EnrichedHeatmap objects (without heatmap legends) from
# the output of countmatrix(). Returns a list named by plot_title.
plotpairheatmap2 <- function(normMatrix, ifng_i, veh_j, color1, color2, ylim=0.5,
                             legend_direction = "vertical") {
  EH_plot <- list()
  #ifng
  EH_plot[[1]] <- EnrichedHeatmap( mat = normMatrix[[1]],
                                   pos_line = FALSE,
                                   border = FALSE,
                                   col = color1,
                                   column_title = plot_title[ifng_i],
                                   column_title_gp = gpar(fontsize = 8, fontfamily = "sans"),
                                   rect_gp = gpar(col = "transparent"),
                                   heatmap_legend_param = list(
                                     legend_direction = legend_direction,
                                     title = NULL,
                                     legend_height = unit(3, "in")),
                                   show_heatmap_legend = FALSE,
                                   top_annotation = HeatmapAnnotation(
                                     enriched = anno_enriched(
                                       gp = gpar(col = "black", lty = 1, lwd=2),
                                       col="black",
                                       ylim = c(0.02, ylim),
                                       axis = FALSE,
                                       legend_width = unit(3, "cm"),
                                       title_position = "topcenter")

                                   )
  )

  #veh
  EH_plot[[2]] <- EnrichedHeatmap( mat =  normMatrix[[2]],
                                   pos_line = FALSE,
                                   border = FALSE,
                                   col = color2,
                                   column_title = plot_title[veh_j],
                                   column_title_gp = gpar(fontsize = 8, fontfamily = "sans"),
                                   rect_gp = gpar(col = "transparent"),
                                   heatmap_legend_param = list(
                                     legend_direction = legend_direction,
                                     title = NULL,
                                     legend_height = unit(3, "in")),
                                   show_heatmap_legend = FALSE,
                                   top_annotation = HeatmapAnnotation(
                                     enriched = anno_enriched(
                                       gp = gpar(col = "black", lty = 1, lwd=2),
                                       col="black",
                                       ylim = c(0.02, ylim),
                                       axis = FALSE,
                                       legend_width = unit(3, "cm"),
                                       title_position = "topcenter")

                                   )

  ) #

  names(EH_plot) <- c(plot_title[ifng_i], plot_title[veh_j])
  return(EH_plot)
}
