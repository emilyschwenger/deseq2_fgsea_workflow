## Custom barplots for fGSEA output ##

# Color by TRUE/FALSE pCutoff
plotBarEnrichment_pCutoff_Publication <-
  function(fgseaRes, pCutoff) {
    myColors <-  c("#E69F00", "grey33")
    names(myColors) <- rev(levels(factor(c("TRUE", "FALSE"))))
    colScale <- scale_fill_manual(name = paste0("padj < ", pCutoff) ,values = myColors)
    ggplot(fgseaRes, aes(reorder(pathway, NES), NES)) +
      geom_col(aes(fill=factor(padj<pCutoff))) +
      coord_flip() +
      colScale +
      #scale_y_discrete(labels = wrap_format(40)) +
      labs(x="Pathway", y="Normalized Enrichment Score",
           title=paste0(title, " pathways NES from GSEA")) +
      theme_classic() +
      theme(plot.title = element_text(face = "bold",  hjust = 0.5),
            axis.title = element_text(face = "bold"),
            panel.border = element_rect(colour = "black", fill=NA, size=1))
  }

# Color by numeric adjusted p-value
plotBarEnrichment_Publication <-
  function(fgseaRes) {
    ggplot(fgseaRes, aes(reorder(pathway, NES), NES)) +
      geom_col(aes(fill=padj)) +
      coord_flip() +
      scale_fill_viridis(option="viridis") + #inferno, magma
      labs(x="Pathway", y="Normalized Enrichment Score",
           title=paste0(title, " pathways NES from GSEA")) +
      theme_classic() +
      theme(plot.title = element_text(face = "bold", hjust = 0.5),
            axis.title = element_text(face = "bold"),
            panel.border = element_rect(colour = "black", fill=NA, size=1))
  }

