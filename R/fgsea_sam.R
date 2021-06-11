### FGSEA ###
## Source: https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html

## Load packages
library(fgsea)
library(here)
library(tidyverse)
library(data.table)
library(ggplot2)
library(viridis)
#library(grafify)
library("colorblindr")
library("scales")
#library("enrichplot")
#library(egg)

## Colors for plots
pub.colors = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff66","#6600cc",
               "#cccccc","#666666","#336633","#CCCC00","000000", "turquoise4")
show_col(palette_OkabeIto)
show_col(pub.colors)

## Function to convert first letter of string to caps
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

## Function to nicely wrap title for faceted plots to avoid cutting off title
wrapper <- function(x, ...)
{
  paste(strwrap(x, ...), collapse = "\n")
}
## Example usage
# my_title <- "This is a really long title of a plot that I want to nicely wrap and "
# p1 +
#   geom_smooth() +
#   ggtitle(wrapper(my_title, width = 60))

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

## Custom plotEnrichment function for fGSEA output
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
      geom_hline(yintercept = max(tops), size = 0.75, colour = "turquoise4",linetype = "dashed") +
      geom_hline(yintercept = min(bottoms), size = 0.75, colour = "turquoise4", linetype = "dashed") +
      geom_hline(yintercept = 0, colour = "black") +
      geom_line(size = 1.5, color = "#CCCC00") + #"#D55E00" "#7fc97f"
      geom_segment(data = data.frame(x = pathway),
                   mapping = aes(x = x, y = -diff/2, xend = x, yend = diff/2),
                   size = ticksSize) +
      xlab("Rank") +
      ylab("Enrichment Score") +
      theme_classic() +
      theme(plot.title = element_text(face="bold"),
            axis.title.x = element_text(face="bold"),
            axis.title.y = element_text(face="bold"),
            panel.border = element_rect(colour = "black", fill=NA, size=1))
    g
  }

## Function to round to the nearest 5
# mround <- function(x,base){
#   base*round(x/base)
# }

#####################################################################################################
#############################################################################
## Load pre-ranked data from DESeq2 output and format for fgsea
# Specify cell line
# cl <- "molm13"
cl <- "ure"

# Specify input and output directories
dir.in <- "~/Desktop/Sam/deseq2/"
dir.out <- "~/Desktop/Sam/gsea/"
dir.create(paste0("~/Desktop/Sam/gsea/", cl))
dir.create(paste0("~/Desktop/Sam/gsea/", cl, "/output_tables"))
dir.create(paste0("~/Desktop/Sam/gsea/", cl, "/fgsea_tables"))
dir.create(paste0("~/Desktop/Sam/gsea/", cl, "/barplots"))
dir.create(paste0("~/Desktop/Sam/gsea/", cl, "/enplots"))

# Read in Sam's data
ranks <- read.table(paste0(dir.in, "tables/DEseq2_", cl, "_res_db2115_vs_veh.rnk"), sep = "\t", header = F, colClasses = c("character", "numeric"))
head(ranks)
colnames(ranks) <- c("ID", "t")
ranks <- setNames(ranks$t, ranks$ID) # can also use deframe here to convert dataframe into named list
str(ranks)
all(!is.na(ranks))

# find the length of the longest label to determine plot width
gmt.all <- here("data/msigdb_genesets/msigdb.v7.4.symbols.gmt")
pathways <- gmtPathways(gmt.all)
max_length <- max(str_length(names(pathways))) # used for padding y-axis labels to ensure consistent plot size

## Load pathways
gmt.files <- list.files(here("data/msigdb_genesets"), full.names = T, recursive = F, include.dirs = F)
gmt.files <- gmt.files[!grepl(gmt.all, gmt.files)]
# gmt.files <- gmt.files[grep("c2", gmt.files)]
# gmt.files <- here("data/msigdb_genesets/h.all.v7.4.symbols.gmt")
gmt.files <- here("data/msigdb_genesets/custom/chip_pu1_promoter_ure.gmt")
gmt.files
pathways <- gmtPathways(gmt.files)

#for(i in 1:2){
for(i in 1:length(gmt.files)){
  tryCatch({
    message("Beginning analysis for gene set  ", gmt.files[i], "...")
    pathways <- gmtPathways(gmt.files[i])
    name <- gsub("\\.v7*.*", "", gmt.files[i])
    name <- gsub(".*/", "", name)

    message("Running fgseaMultilevel for gene set  ", gmt.files[i], "...")
    fgseaRes <- fgseaMultilevel(pathways, ranks, minSize=15, maxSize=500)
    message("fgseaMultilevel for gene set  ", gmt.files[i], " complete at: ", date())

    fgseaResOut <- fgseaRes %>%
      as_tibble() %>%
      arrange(padj) %>%
      slice_head(n = 200) %>%
      mutate(leadingEdge = sapply(leadingEdge, toString))
    fgseaResTidy <- fgseaRes %>%
      as_tibble() %>%
      arrange(padj) %>%
      slice_head(n = 40) %>%
      mutate(pathway = str_pad(pathway, max_length, side = "left"))   # pad the shorter labels to ensure consistent plot width
    print(head(fgseaResOut))
    write.table(fgseaResOut, file = paste0(dir.out, cl, "/output_tables/", cl, "_", name, "_fgsea_output.tsv"), sep = "\t", quote = F, row.names = F)

    ## Generate barplots ##
    message("Generating barplot for  ", gmt.files[i], "...")
    title <- toupper(gsub("\\.", " ", name))
    message("Title is: ", title)

    p <- plotBarEnrichment_pCutoff_Publication(fgseaResTidy, 0.05)
    q <- plotBarEnrichment_Publication(fgseaResTidy)

    ggsave(paste0(dir.out, cl, "/barplots/", cl, "_", name, "_gsea_barplot_pval0.05.png"), plot = p, dpi = 300, width = 15, height = 6)
    ggsave(paste0(dir.out, cl, "/barplots/", cl, "_", name, "_gsea_barplot.png"), plot = q, dpi = 300, width = 15, height = 6)

    # Plot awesome table w/ plots
    if(nrow(fgseaRes[ES > 0])>=20){num.pos = 20}
    else{num.pos = nrow(fgseaRes[ES > 0])}
    if(nrow(fgseaRes[ES < 0])>=20){num.neg = 20}
    else{num.neg = nrow(fgseaRes[ES < 0])}
    topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=num.pos), pathway]
    topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=num.neg), pathway]
    topPathways <- c(topPathwaysUp, rev(topPathwaysDown))

    pdf(paste0(dir.out, cl, "/fgsea_tables/", cl, "_", name, "_gsea_top20_table.pdf"), width = 20, height = 10)
    plotGseaTable(pathways[topPathways], ranks, fgseaRes,
                  gseaParam = 0.5)
    dev.off()

    # Generate enplots for top 20 gene lists
    for(j in 1:num.pos){
      assign(paste0("p",j), plotEnrichment_Publication(pathways[[topPathwaysUp[j]]], ranks) +
               ggtitle(wrapper(gsub("_", " ", topPathwaysUp[j]), width = 28)) )
    }

    for(j in 1:num.neg){
      assign(paste0("q",j), plotEnrichment_Publication(pathways[[topPathwaysDown[j]]], ranks) +
               ggtitle(wrapper(gsub("_", " ", topPathwaysDown[j]), width = 28)) )
    }

    p <- cowplot::plot_grid(plotlist = mget(paste0("p", 1:num.pos)), ncol=4) #, labels=LETTERS[1:20])
    q <- cowplot::plot_grid(plotlist = mget(paste0("q", 1:num.neg)), ncol=4) #, labels=LETTERS[1:20])

    ggsave(paste0(dir.out, cl, "/enplots/", cl, "_", name, "_top",num.pos, "Up_fgsea_enplots.png"), plot = p, dpi = 300, width = 15, height = 3.5*ceiling(num.pos/4) )
    ggsave(paste0(dir.out, cl, "/enplots/", cl, "_", name, "_top", num.neg, "Down_fgsea_enplots.png"), plot = q, dpi = 300, width = 15, height = 3.5*ceiling(num.pos/4) )
    rm(list = ls(pattern = "[p-q][1-9]"))

    message("Analysis of  ", name, " completed at: ", date())
    message("=====================================================================================================================================\n")
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")}
  )
}

###################################################################################################
########################################################################################
################################################################################
## Run one at a time ##

gmt.files <- here("data/msigdb_genesets/custom/chip_cutntag_all.gmt")
pathways <- gmtPathways(gmt.files[2])

# pathways <- gmtPathways(system.file(
#   "extdata", "mouse.reactome.gmt", package="fgsea"))
str(head(pathways))

# Show the first few pathways, and within those, show only the first few genes.
# pathways %>%
#   head() %>%
#   lapply(head)

## Run fgsea
fgseaRes <- fgsea(pathways, ranks, minSize=15, maxSize=500) #, nperm = 10000) only use nperm for fgseaSimple (NOT recommended, multilevel recommended)
fgseaRes <- fgseaMultilevel(pathways, ranks, minSize=15, maxSize=500) # default; default sample size = 101 (do not change)
# fgseaRes <- fgseaSimple(pathways, ranks, minSize=15, maxSize=500)
# fgseaRes[grep("stem", fgseaRes@pathway, ignore.case = T),]
# fgseaRes.all <- fgseaRes
# top 6 enriched pathways
head(fgseaRes[order(pval), ], 100)

# number of significant pathways at padj < 0.01
sum(fgseaRes[, padj < 0.01])
# [1] 1116

fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  #arrange(desc(NES)) %>%
  arrange(padj) %>%
  slice_head(n = 40)

fgseaResTidy %>%
  mutate(b = sapply(leadingEdge, toString))

fgseaRes <- fgseaRes %>%
  as_tibble() %>%
  arrange(padj) %>%
  mutate(leadingEdge = cat(unlist(leadingEdge), sep =','))
write.table(fgseaRes, file = paste0(dir.out, cl, "/output_tables/", cl, "_", name, "_fgsea_output.tsv"), sep = "\t", quote = F, row.names = F)


# Show in a nice table:
# fgseaResTidy %>%
#   dplyr::select(-leadingEdge, -ES) %>%
#   arrange(padj) %>%
#   DT::datatable()

# Barplot
#enrichDF2enrichResult(fgseaRes)
#emapplot(fgseaRes, showCategory=20)

p <- ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=factor(padj<0.05, levels = c("TRUE", "FALSE")))) +
  #geom_col(aes(fill=padj)) +
  coord_flip() +
  scale_fill_manual("padj < 0.05", values = c("#E69F00", "grey33")) + #palette_OkabeIto
  #scale_fill_OkabeIto() +
  #scale_fill_viridis(option="viridis") + #inferno, magma
  labs(x="Pathway", y="Normalized Enrichment Score",
       title=paste0(name, " pathways NES from GSEA")) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold",  hjust = 0.5),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
p

q <- ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  #geom_col(aes(fill=padj<0.05)) +
  geom_col(aes(fill=padj)) +
  coord_flip() +
  #scale_fill_manual(values = pub.colors) +
  #scale_fill_OkabeIto() +
  scale_fill_viridis(option="viridis") + #inferno, magma
  labs(x="Pathway", y="Normalized Enrichment Score",
       title=paste0(name, " pathways NES from GSEA")) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
q
ggsave(paste0("~/Desktop/Sam/gsea/barplots/molm13_", name, "_gsea_barplot_pval0.05.png"), width = 7.5, height = 6, plot = p, dpi = 300)
ggsave(paste0("~/Desktop/Sam/gsea/barplots/molm13_", name, "_gsea_barplot.png"), width = 7.5, height = 6, plot = q, dpi = 300)

# plot the most significantly enriched pathway
plotEnrichment(pathways[[head(fgseaRes[order(pval), ], 1)$pathway]],
               ranks) + labs(title=head(fgseaRes[order(pval), ], 1)$pathway)

pl <- "HALLMARK_COAGULATION"


pl1 <- "RGAGGAARY_PU1_Q6"
pl2 <- "ETS2_B"
pl3 <- "ETS_Q4"
pl4 <- "WGGAATGY_TEF1_Q6"

p1 <- plotEnrichment(pathways[[pl1]],
                     ranks) + labs(title=pl1)
p2 <- plotEnrichment(pathways[[pl2]],
                     ranks) + labs(title=pl2)
p3 <- plotEnrichment(pathways[[pl3]],
                     ranks) + labs(title=pl3)
p4 <- plotEnrichment(pathways[[pl4]],
                     ranks) + labs(title=pl4)

p <- cowplot::plot_grid(p1, p2, p3, p4, ncol=2, labels=LETTERS[1:4])

# dir.create("~/Desktop/Sam/gsea/enplots")
ggsave(paste0("~/Desktop/Sam/gsea/enplots/molm13_", name, "_gsea_enplots.png"), plot = p, dpi = 300) #width = 5, height = 4,


# Plot awesome table w/ plots
topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=20), pathway]
topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=20), pathway]
topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
dev.off()

dir.create("~/Desktop/Sam/gsea/tables")
pdf(paste0("~/Desktop/Sam/gsea/tables/molm13_", name, "_gsea_top20_table.pdf"))
plotGseaTable(pathways[topPathways], ranks, fgseaRes,
              gseaParam = 0.5)
dev.off()
