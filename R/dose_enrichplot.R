


library(DOSE)
data(geneList)
head(geneList)

de <- names(geneList)[abs(geneList) > 2] # given a vector of genes, this function will return the enrichment NCG categories with FDR control
de <- names(ranks) # given a vector of genes, this function will return the enrichment NCG categories with FDR control

edo <- enrichDGN(de, )
edo2@result

## Barplot ##
library(enrichplot)
barplot(edo, showCategory=20)


## Dotplot ##
edo2 <- gseNCG(geneList, nPerm=10000)
p1 <- dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA")
p2 <- dotplot(edo2, showCategory=30) + ggtitle("dotplot for GSEA")
plot_grid(p1, p2, ncol=2)


## Gene concept network ##
## convert gene ID to Symbol
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=geneList)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
p3 <- cnetplot(edox, foldChange=geneList, circular = TRUE, colorEdge = TRUE)
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))