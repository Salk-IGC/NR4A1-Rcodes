# Plot helper functions shared by the figure scripts
#
# Required packages: ggplot2, ggrepel (ggrepel is attached by EnhancedVolcano)

# save_pheatmap_jpeg(x, filename, width, height)
# Saves a pheatmap object as a 600 dpi JPEG.
save_pheatmap_jpeg <- function(x, filename, width = 7, height = 7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  jpeg(filename, width = width, height = height, units = "in",
       res = 600, type = "cairo")
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

# gVolcano4review(rl, geneList, name)
# Volcano plot of a DESeq2 results table: log2 fold change vs -log10 adjusted
# p-value, coloured by significance (|log2FC| > 1 and padj < 0.05). The genes in
# geneList (the top 30 by padj when geneList is NULL) are highlighted and labelled.
# Returns the ggplot object.
gVolcano4review <- function(rl = NULL, geneList = NULL, name){

  # handle zero padj
  if(min(rl$padj, na.rm = TRUE) == 0){
    xxx <- subset(rl, padj != 0)
    aaa <- min(xxx$padj)
    rl[!is.na(rl$padj) & rl$padj == 0, ]$padj <- (aaa / 10000000)
  }

  minx <- min(rl$log2FoldChange, na.rm = TRUE)
  maxx <- max(rl$log2FoldChange, na.rm = TRUE)
  maxy <- -log10(min(rl$padj, na.rm = TRUE))

  rl$sigDif <- "NS"
  rl$sigDif[rl$log2FoldChange > 1  & rl$padj < 0.05] <- "UP"
  rl$sigDif[rl$log2FoldChange < -1 & rl$padj < 0.05] <- "DOWN"
  rl$sigDif <- factor(rl$sigDif, levels = c("DOWN", "NS", "UP"))
  print(table(rl$sigDif))

  cp  <- c(`DOWN` = "#00AFBB", `NS` = "grey", `UP` = "#bb0c00")
  lbs <- c(`DOWN` = "Downregulated", `NS` = "Not significant", `UP` = "Upregulated")

  # highlighted genes (significant or not)
  if(is.null(geneList)){
    top <- head(rl[order(rl$padj), "gene_name"], 30)
    hl <- rl[rl$gene_name %in% top, ]
  } else {
    hl <- rl[rl$gene_name %in% geneList, ]
  }

  # report any requested genes missing from the results
  if(!is.null(geneList)){
    missing <- setdiff(geneList, rl$gene_name)
    if(length(missing) > 0) message("Not found in results: ", paste(missing, collapse = ", "))
  }

  myvolcanoplot <- ggplot(rl, aes(x = log2FoldChange, y = -log10(padj), col = sigDif)) +
    geom_vline(xintercept = c(-1, 1), col = "gray", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), col = "gray", linetype = "dashed") +
    # background points (made slightly transparent so highlights pop)
    geom_point(size = 2, alpha = 0.6) +
    # highlighted points: big, black outline, drawn on top
    geom_point(data = hl, aes(fill = sigDif),
               shape = 21, size = 4.5, stroke = 1.2, color = "black",
               show.legend = FALSE) +
    # highlighted labels: bold, pushed up and away from the crowd
    geom_label_repel(data = hl, aes(label = gene_name),
                     color = "black", fontface = "bold", size = 4,
                     fill = "white", label.size = 0.6,
                     box.padding = 1, point.padding = 0.8,
                     min.segment.length = 0, segment.color = "black",
                     segment.size = 0.6,
                     nudge_y = maxy * 0.15,
                     max.overlaps = Inf, show.legend = FALSE) +
    scale_color_manual(values = cp, labels = lbs) +
    scale_fill_manual(values = cp) +
    coord_cartesian(ylim = c(0, maxy + 5), xlim = c(minx, maxx)) +
    labs(color = "Significance",
         x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-adj")) +
    ggtitle(name) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      panel.grid.major = element_line(color = "grey80"),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white")
    )

  return(myvolcanoplot)
}
