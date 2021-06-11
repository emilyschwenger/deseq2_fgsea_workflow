### DESeq2 Analysis of db2115-Treated Cells vs. Control ###

## Renv commands to ensure reproducibility of code
# renv::init() # initialize r environment and automatically install packages declared in the previously stored lockfile
# renv::snapshot() # saves exact versions of libraries to ensure reproducibility
# renv::history() # view past versions of renv.lock
# renv::revert() # to pull out an old version of renv based on a previous commit then use...
# renv::restore() # to restore library from that state

## R usethis package allows for setup of Git user name and email from within R
library(usethis)
library(gitcreds)
library(credentials)

## Configure GitHub accounts
use_git_config(user.name = "emilyschwenger", user.email = "emilyschwenger@gmail.com")
usethis::create_github_token()
usethis::use_git()

## Alternatively can configure GitHub account in command line...
# git config --global user.name 'Jane Doe'
# git config --global user.email 'jane@example.com'
# git config --global --list

## You may need to specify your editor of choice...
# git config --global core.editor "emacs" # "vim"

## Set up personal access token for 2-factor authentification
gitcreds_set()
gitcreds_get()
set_github_pat()
# credentials::git_credential_forget() # clears a credential

## Load libraries
library(DESeq2)
library(here)
library(EnhancedVolcano)
library("readr")

## Designate input and output directories
dir.in <- "~/Desktop/Sam/rna-seq/data/raw"
dir.out <- "~/Desktop/Sam/deseq2/"

## Designate which cell lines to look at + file names
# cl <- "molm13"
cl <- "ure"
files <- sort(list.files(dir.in, pattern = cl, ignore.case = T))
names <- gsub("\\_ReadsPerGene.out.tab","",files)
names <- tolower(names)
names

## Read in .tab files directly (STAR output)
for(i in 1:length(files)){
 assign(names[i],read.table(paste0(dir,files[i]), sep = "\t", header = F))
 tab <- get(names[i])
 assign(names[i], tab[-c(1:4),])
}

## Use tximport if Salmon output
# txi <- tximport(paste0(dir, files.ure), type="salmon", tx2gene=tx2gene)

## Read in Boris's R object
# FileName <- "URE_counts.Rds"
# FileName2 <- "MOLM13_counts.Rds"
# dds.ure <- readRDS(paste0(FilePath, FileName))
# dds.molm <- readRDS(paste0(FilePath, FileName2))
# dds <- dds.ure

## Prepare count file from STAR tab file
# Sanity check: rownames the same across samples
all(ure2115_a$V1 %in% ure2115_b$V1)
all(ure2115_a$V1 %in% ure2115_v$V1)
all(mol2115_a$V1 %in% molveh_b$V1)

## Function to format STAR .tab output into one matrix for DESeq2 analysis
# Extract 4th column (reverse-stranded) based on strand-specific protocol
starFormat <- function(x) {
  x <- get(x)
  x.sub <- data.frame(x[,4])
  row.names(x.sub) <- x[,1]
  #print(dim(x.sub))
  #print(head(x.sub))
  return(x.sub)
}
# test <- starFormat("ure2115_a")
# dim(test)

countFile <- lapply(names, starFormat)
countFile <- as.data.frame(do.call(cbind, countFile))
dim(countFile)
colnames(countFile) <- names
tail(countFile)
countData <- data.matrix(countFile)

## Prepare coldata
colData = data.frame(colnames(countData),c(rep("DB2115",3),rep("Vehicle",3)), rep(c("B1","B1","B2"),2))
colnames(colData) = c("Sample","Treatment", "Batch")
rownames(colData) = colData$Sample
colData

all(rownames(colData) %in% colnames(countData))

### Design factors ###

# The simplest design formula for differential expression would be ∼ condition, where condition is a
# column in colData(dds) which specifies which of two (or more groups) the samples belong to. For
# the parathyroid experiment, we will specify ∼ patient + treatment, which means that we want to
# test for the effect of treatment (the last factor), controlling for the effect of patient (the first factor).

# ddsFull <- DESeqDataSet( se, design = ~ patient + treatment )
# dds2 <- DESeqDataSetFromMatrix(countData = countData11,
#                              colData = colData11, design = ~ Treatment + Stimulation)
# dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~ Genotype + Treatment +
# Genotype:Treatment)

dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData, design = ~ Batch + Treatment)
dds$Treatment <- relevel(dds$Treatment, ref = "Vehicle")
dds

#Minimal pre-filtering to remove rows that have only 0 or 1 read
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds
dds <- DESeq(dds)
res <- results(dds)
res

# Show size factors using the "median ratio method" described by Equation 5 in Anders and Huber (2010)
sizeFactors(dds)
# mol2115_a mol2115_b mol2115_v  molveh_a  molveh_b  molveh_v
# 0.8161427 0.8201502 1.6494480 0.7103518 0.8415737 1.5586588
## NOTE: batch2 ("v") has overall much more RNA expression

# ure2115_a ure2115_b ure2115_v  ureveh_a  ureveh_b  ureveh_v
# 1.0375677 1.0202100 1.1053600 0.9868409 0.8297385 1.0552639
# esf <- estimateSizeFactors(dds)
normCt <- counts(dds, normalized=TRUE)

#Extracting transformed values
rld <- rlog(dds, blind=FALSE)
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
vsd.fast <- vst(dds, blind=FALSE)
head(assay(rld), 3)
head(assay(vsd))

## Plot distribution of counts ##
hist(as.numeric(countData[,1]),breaks=10000,xlim=c(0,5000),ylim=c(0,1000))
hist(as.numeric(normCt[,4]),breaks=10000,xlim=c(0,5000),ylim=c(0,1000))

## Plot dispersion of data ##
plotDispEsts(dds)

## Plot filter stats ##
# We can visualize the optimization by plotting the filterNumRej attribute of the results object.
# The results function maximizes the number of rejections (adjusted p value less than a significance level),
# over the quantiles of a filter statistic (the mean of normalized counts).
# The threshold chosen (vertical line) is the lowest quantile of the filter for which the number of rejections is
# within 1 residual standard deviation to the peak of a curve fit to the number of rejections over the filter quantiles:

metadata(res)$alpha
metadata(res)$filterThreshold
# 25.20408%
# 3.301892
plot(metadata(res)$filterNumRej,
     type="b", ylab="number of rejections",
     xlab="quantiles of filter")
lines(metadata(res)$lo.fit, col="red")
abline(v=metadata(res)$filterTheta)
##########################################################
# Order our results table by the smallest adjusted p value
resOrdered <- res[order(res$padj),]
summary(res)

# How many adjusted p-values were less than 0.1?
sum(res$padj < 0.1, na.rm=TRUE)

res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)
# out of 24428 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2692, 11%
# LFC < 0 (down)     : 3875, 16%
# outliers [1]       : 0, 0%
# low counts [2]     : 6631, 27%
# (mean count < 4)

dir.create(paste0(dir.out,"plots"))
pdf(paste0(dir.out, "plots/", cl, "_db2115_vs_veh_plotMA.pdf"))
plotMA(res, main="DESeq2", ylim=c(-2,2))
dev.off()

#After calling plotMA, one can use the function identify to interactively detect the
#row number of individual genes by clicking on the plot. One can then recover the
#gene identifiers by saving the resulting indices:

# idx <- identify(res$baseMean, res$log2FoldChange)
# rownames(res)[idx]

#resMLE <- results(dds, addMLE=TRUE)
#head(resMLE, 4)
#One can make an MA-plot of the unshrunken estimates like so:
#plotMA(resMLE, MLE=TRUE, main="unshrunken LFC", ylim=c(-2,2))

#Plot Counts by Treatment
# plotCounts(dds, gene=which.min(res05$padj), intgroup="Treatment")

#For customized plotting, an argument returnData specifies that the function should
#only return a data.frame for plotting with ggplot.
# d <- plotCounts(dds, gene=which.min(res$padj), intgroup="MDMX.Tg",
#                 returnData=TRUE)
# library("ggplot2")
# ggplot(d, aes(x=MDMX.Tg, y=count)) +
#   geom_point(position=position_jitter(w=0.1,h=0)) +
#   scale_y_log10(breaks=c(25,100,400))

#More info on results column
# mcols(res)$description
# res[1:6]

## Specify species for ensembl id to gene conversion
# species = "human"
species = "mouse"

## Convert ensembl id to gene symbol
library(biomaRt)
if(species=="human"){dat = "hsapiens_gene_ensembl"; sym = "hgnc_symbol"}
if(species=="mouse"){dat = "mmusculus_gene_ensembl"; sym = "mgi_symbol"}
ensembl <- useMart("ensembl", dataset=dat)
annot <-getBM(c("ensembl_gene_id", sym, "chromosome_name", "strand", # use listAttributes(ensembl) for more gene identifiers
                "start_position", "end_position","gene_biotype","description"), mart=ensembl)
colnames(annot)[grep(sym, colnames(annot))] <- "gene_symbol"

# Save normalized counts based on geometric mean of row + division by median (of all rows) and eliminate number after "."
normCt <- as.data.frame(counts(dds, normalized=TRUE))
normCt$ensembl_gene_id <- gsub("*\\..*","", rownames(normCt))
res.write = data.frame(resOrdered)
res.write$ensembl_gene_id <- gsub("*\\..*","", rownames(res.write))

# Merge results file with gene annotation from BioMart and normalized counts
res.write <- res.write %>%
  as_tibble() %>%
  left_join(annot, by = 'ensembl_gene_id') %>%
  left_join(normCt, by = 'ensembl_gene_id')

# Exporting results
dir.create(paste0(dir.out, "tables"))
write.table(res.write, file=paste0(dir.out, "tables/DEseq2_", cl, "_res_db2115_vs_veh.tsv"), quote=F, sep="\t", row.names = F)

# Eliminate rows with no row name
res.gsea <- res.write[!(res.write$gene_symbol == ""),]

# Reformat results for GSEA preranked
res.gsea$rank <- (-log(x = res.gsea$padj,base = 10)*sign(res.gsea$log2FoldChange))

# Assign infinity values to the numeric maximum
rnk.max <- max(res.gsea$rank[!abs(res.gsea$rank)==Inf], na.rm = T)
pos.max <- grep(rnk.max, res.gsea$rank)
pos <- grep(Inf, res.gsea$rank)
pos.neg <- which(res.gsea$rank==-Inf)
res.gsea$rank[pos] <- rev(seq(from=rnk.max+0.1, by=0.1, length.out=length(pos)))
res.gsea$rank[pos.neg] <- -res.gsea$rank[pos.neg]

# Eliminate NAs and average duplicates
res.gsea <- res.gsea %>%
  dplyr::select(gene_symbol, rank) %>%
  mutate(gene_symbol = toupper(gene_symbol)) %>%
  na.omit() %>%
  distinct() %>%
  group_by(gene_symbol) %>%
  summarize(rank=mean(rank))
res.gsea <- res.gsea[order(res.gsea$rank, decreasing = T),]
res.gsea

write.table(res.gsea, file=paste0(dir.out, "tables/", "DEseq2_", cl, "_res_db2115_vs_veh.rnk"), quote = F, sep="\t", row.names = F, col.names = F)

## Generate publication-ready volcano plot
library(EnhancedVolcano)
title = paste0(toupper(cl), "-/-: DB2115 vs. Vehicle")
# myColors <- c("grey30", "#E69F00", "#56B4E9", "#009E73", "#F0E442") # c("grey30", "#FF0000", "royalblue", "#00A08A")
pdf(paste0(dir.out, "/plots/", "DEseq2_", cl, "_db2115_vs_veh_volcano_gene_symbols2.pdf"), width = 8, height = 6)
EnhancedVolcano(res.write,
                lab = res.write$gene_symbol,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                #drawConnectors = T,
                xlim = c(-max(abs(res.write$log2FoldChange)), max(abs(res.write$log2FoldChange))),
                title = title,
                subtitle = bquote(italic("Differential expression")),
                col = c("grey30", pub.colors[2:4])) +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(face = "bold",  hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
# c('grey30', 'forestgreen', 'royalblue', 'red2').
dev.off()

## Effects of transformations on variance
library("vsn")
notAllZero <- (rowSums(counts(dds))>0)
meanSdPlot(log2(counts(dds,normalized=TRUE)[notAllZero,] + 1))
meanSdPlot(assay(rld[notAllZero,]))
meanSdPlot(assay(vsd[notAllZero,]))

## Heatmap of the count matrix
# Function to save pheatmap output to pdf
save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
library("pheatmap")
for(i in c(50, 100, 200)) {
  num <- i
  res.plot <- res.write
  res.plot = res.plot[1:num, c('gene_symbol', colnames(res.plot)[grep(cl, colnames(res.plot))])]
  # is.na(res.plot$gene_symbol) %>% # sanity checks
  #   table()
  # is.infinite(res.plot$ureveh_v) %>%
  #   table()
  res.plot <- res.plot[!is.na(res.plot$gene_symbol),]
  res.plot <- data.frame(res.plot)
  rownames(res.plot) <- make.unique(res.plot$gene_symbol)
  res.plot <- res.plot[,-1]
  # all(rowSums(res.plot)>0)
  # all(colSums(res.plot)>0)
  rm(p)
  p <- pheatmap(res.plot,
                cluster_rows=TRUE,
                cluster_cols=FALSE,
                show_rownames=TRUE,
                scale = "row",
                border_color = NA,
                annotation_col=colData)
  save_pheatmap_pdf(p, paste0(dir.out, "plots/", "DEseq2_", cl, "_db2115_vs_veh_pheatmap_", num, ".pdf"), width = 10, height = num/7.5)
}


## Heatmap of the sample-to-sample distances
library("RColorBrewer")
sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rld$Sample
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
p <- pheatmap(sampleDistMatrix,
              clustering_distance_rows=sampleDists,
              clustering_distance_cols=sampleDists,
              col=colors)
save_pheatmap_pdf(p, paste0(dir.out, "plots/", "DEseq2_", cl, "_db2115_vs_veh_sampleDist.pdf"))

## Custom PCA plot with ggplot
data <- plotPCA(rld, intgroup=c('Treatment', 'Batch'), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
p <- ggplot(data, aes(PC1, PC2, color=Treatment, shape=Batch)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  coord_fixed() +
  scale_color_manual(values = pub.colors) +
  ggtitle(title) +
  theme_classic() +
  theme(plot.title = element_text(face = "bold",  hjust = 0.5),
        axis.title = element_text(face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1))
  #theme_Publication() +
  #scale_colour_Publication()
p
ggsave(paste0(dir.out, "plots/", "DEseq2_", cl, "_db2115_vs__veh_pca2.png"), p, dpi = 300, width = 5,height = 4)
# ggsave("pcaPlot_pub2.png", plot = g, device = NULL, path = "~/Desktop/Koki/for_paper/",
#        scale = 1, width = 10, height = 8, units = c("in", "cm", "mm"),
#        dpi = 2000, limitsize = TRUE)

## PCAtools to generate loadings plot ##
# library("PCAtools")
# #dat <-
# # res.write[,c(14:19)] %>%
# #   as_tibble() %>%
#
# dat <- res.write
# dat <- dat[!is.na(res.write$hgnc_symbol),]
# dat
# rownames(dat) <- make.unique(res.write$hgnc_symbol)
# p <- pca(countFile, metadata = colData, removeVar = 0.1)
# p <- pca(countFile, metadata = colData, removeVar = 0.1)
# p
#
# # Scree plot
# screeplot(p)
#
# # Biplot
# biplot(p)
#
# #colDataB$Disease <- toupper(colDataB$Disease)
# #colDataB$Disease <- gsub(" ","", colDataB$Disease)
# pdf("~/Desktop/Brad_Tricomi/pcaTools_biplot_updated_noLabels.pdf",width = 8,height = 6)
# biplot(p,colby = 'Diagnosis',colkey = c("MDS"="#386cb0","MF"="#fdb462",#"MDS/AML"="#7fc97f","MF/MPN"="#ef3b2c",
#                                         "Control"="grey69"), legendPosition = 'right', lab = NULL)
# dev.off()
# # shape = as.character('Anemic'),
# # values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
#
# # Pairs plot
# pairsplot(p)
#
# # Loadings plot
# plotloadings(p)
#
# # Eigencore plot
# eigencorplot(p1,
#              metavars = c('Age','Distant.RFS','ER','GGI','Grade','Size','Time.RFS'))
#
# eigencorplot(p,
#              metavars = c('HEP')) #,'Anemic')) #,'sample','background'))





#VARIATIONS TO THE STANDARD WORKFLOW

#Wald test individual steps
#The function DESeq runs the following functions in order:
# dds <- estimateSizeFactors(dds)
# dds <- estimateDispersions(dds)
#
#
#
# chip <- readxl::read_xlsx("~/Desktop/Sam/chip-seq_pos_rna-seq_pos.xlsx")
# #ovl <- chip$`All genes`[chip$`All genes` %in% chip$KRT80]
# chip$`All genes`[chip$`All genes` %in% res.write$hgnc_symbol[res.write$padj<0.05]]
# ovl <- chip$`All genes`[chip$`All genes` %in% res.write$hgnc_symbol[res.write$padj<0.05 & res.write$log2FoldChange>0]]
# #res.write[res.write$padj<0.05 & res.write$log2FoldChange>0,]
# write.table(ovl, "~/Desktop/Sam/molm13_chip-seq_pos_rna-seq_pos_overlap.txt", sep = "\t", quote = F, row.names = F)
#
# dds <- nbinomWaldTest(dds)
