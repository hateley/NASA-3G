# source("3G1G_DESeq2.R")
# This is DESeq2 pipeline for NASA 3G samples. This is the final analysis.

library("DESeq2")
library("vsn")
library("RColorBrewer")
library("gplots")

dataPrep <- function(counts) {
    # Read in source data.
    # I used a table created earlier with merged eXpress gene TPM for all 8 samples.
    print("Creating source data")
    countsTable <<- read.delim(counts,header=TRUE,row.names=1)

    # Create metadata
    # The two conditions are 1x gravity and 3x gravity. 4 samples each.
    hyperG_design <<- data.frame(
        row.names = colnames(countsTable),
        condition = c('1G','1G','1G','1G','3G','3G','3G','3G'))

    #create dataset, using our single condition of gravity variation.
    dds <<- DESeqDataSetFromMatrix(
        countData = countsTable,
        colData = hyperG_design,
        design = ~ condition)
    colData(dds)$condition <<- factor(colData(dds)$condition,
    	levels=c('1G','3G'))

    #prep blind data using regularized log transformation
    rld <<- rlog(dds)
    vsd <<- varianceStabilizingTransformation(dds)

}

runDESeq <- function() {
    # running differential expression steps
    # DESeq performs a bundled analysis, res generates table of the results.
    print("Performing DE-Seq")
    dds_de <<- DESeq(dds)
    res <<- results(dds_de, addMLE=TRUE)
    res <<- res[order(res$padj),]
#    print("Outputting DESeq results to file 'deseq2_degenes_3G1G.csv'")
#    write.csv(as.data.frame(res),file="deseq2_degenes_3G1G.csv")
    print("Here is a summary of the DE results:")
    summary(res, alpha = 0.05)

    #output significant genes with padj < 0.05
    resSig <<- subset(res, padj < 0.05)
#    print("Outputting significantly DE genes with padj < 0.05 to file 'deseq2_sig.05_3G1G.csv'")
#    write.csv(as.data.frame(resSig),file="deseq2_sig.05_3G1G.csv")

}

plots <- function() {
	print("Making MA plot")
	plotMA(res, main = "MA plot", ylim=c(-2,2))
        print("plotting dispersion estimates")
        plotDispEsts(dds_de)

}

QA <- function() {
	print("Performing quality assurance and outputting QA plots")
        #plot the standard deviation over the mean for shifted log and transformed data
        par(mfrow=c(1,3))
	notAllZero <- (rowSums(counts(dds_de))>0)
	meanSdPlot(log2(counts(dds_de,normalized=TRUE)[notAllZero,] + 1))
	meanSdPlot(assay(rld[notAllZero]), xlab = 'rld')
	meanSdPlot(assay(vsd[notAllZero]), xlab ='vsd')
  #plot heat maps of 30 highest expressed genes
  select <<- order(rowMeans(counts(dds_de,normalized=TRUE)),decreasing=TRUE)[1:30]
  hmcol <<- colorRampPalette(brewer.pal(9,"GnBu"))(100)
  heatmap.2(counts(dds_de,normalized=TRUE)[select,], col=hmcol,Colv = FALSE, scale = "none", dendrogram='none',trace="none", main = 'normalized counts')
  heatmap.2(assay(rld)[select,], col=hmcol, Colv = FALSE, dendrogram='none',trace="none", main = 'rld')
  heatmap.2(assay(vsd)[select,], col=hmcol, Colv = FALSE, dendrogram='none',trace="none", main = 'vsd')
  #plot sample-to-sample Euclidean distance matrix using rld
  distsRL <- dist(t(assay(rld)))
  mat <- as.matrix(distsRL)
  hc <- hclust(distsRL)
  heatmap.2(mat, Rowv=as.dendrogram(hc), symm=TRUE, trace="none", col = rev(hmcol), margin=c(13, 13))
  plotPCA(rld)
}

plotPCAWithSampleNames = function(x, intgroup="condition", ntop=500)

{
  library("genefilter")
  library("lattice")

  rv = rowVars(assay(x))
  select = order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  pca = prcomp(t(assay(x)[select,]))

  # extract sample names
  names = colnames(x)

  fac = factor(apply( as.data.frame(colData(x)[, intgroup, drop=FALSE]), 1, paste, collapse=" : "))
  if( nlevels(fac) >= 3 )
	  colours = brewer.pal(nlevels(fac), "Dark2")
	else
	  colours = c( "mediumpurple1", "darkgoldenrod1" )

  xyplot(
		PC2 ~ PC1, groups=fac, data=as.data.frame(pca$x), xlab=list("PC1: 56% variance", cex=1.4), ylab=list("PC2: 16% variance", cex=1.4), pch=16, cex=3,
		panel=function(x, y, ...) {
			panel.xyplot(x, y, ...);
			ltext(x=x, y=y, labels=names, pos=1, offset=0.1, cex=1.2)
		},
		aspect = "iso", col=colours,
		main = "PCA: 3 g vs 1 g")
}


dataPrep("merged_counts_eXpress_TPM100_genelevel.tsv")
runDESeq()
plots()
#QA()
plotPCAWithSampleNames(rld)
sessionInfo()
