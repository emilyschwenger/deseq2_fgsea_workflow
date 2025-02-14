## Custom plotEnrichment function for fGSEA output ##

plotEnrichment_Publication <-
  function (pathway, stats, gseaParam = 1, ticksSize = 0.2)
  {
    rnk <- rank(-stats)
    ord <- order(rnk)
    statsAdj <- stats[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
    statsAdj <- statsAdj/max(abs(statsAdj))
    pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
    pathway <- sort(pathway)
    gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,
                            returnAllExtremes = TRUE)
    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops
    n <- length(statsAdj)
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
    diff <- (max(tops) - min(bottoms))/8
    x = y = NULL
    g <- ggplot(toPlot, aes(x = x, y = y)) +
      geom_point(color = "green", size = 0.1) +
      geom_hline(yintercept = max(tops), size = 0.75, colour = "darkcyan",linetype = "dashed") +
      geom_hline(yintercept = min(bottoms), size = 0.75, colour = "darkcyan", linetype = "dashed") +
      geom_hline(yintercept = 0, colour = "black") +
      geom_line(size = 1, color = "#CCCC00") + # "#D55E00" "#7fc97f" "#CCCC00" "#E69F00
      geom_segment(data = data.frame(x = pathway),
                   mapping = aes(x = x, y = -diff/2, xend = x, yend = diff/2),
                   size = ticksSize) +
      xlab("Rank") +
      ylab("Enrichment Score") +
      theme_classic() +
      theme(plot.title = element_text(face="bold", hjust = 0.5),
            axis.title.x = element_text(face="bold"),
            axis.title.y = element_text(face="bold"),
            panel.border = element_rect(colour = "black", fill=NA, size=1))
    g
  }
