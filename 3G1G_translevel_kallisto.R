# source("source_if_on_network")
# source("3G1G_translevel_kallisto.R")

library("DESeq")
library("gplots")
library("RColorBrewer")

dataPrep <- function(counts) {
    # Create source data
    print("Creating source data")
    hyperG_countsTable <<- read.delim(counts,header=TRUE,row.names=1)

    # Create metadata
    hyperG_design <<- data.frame(
        row.names = colnames(hyperG_countsTable),
        condition = c('1G','1G','1G','1G','3G','3G','3G','3G'))

    # Prepare data
    print("Preparing data")
    hyperG_cds <- newCountDataSet(hyperG_countsTable, hyperG_design$condition)
    hyperG_cds_size <<- estimateSizeFactors(hyperG_cds)
    hyperG_cds_disp <<- estimateDispersions(hyperG_cds_size)

    # Filter data
    print("Filtering data")
    rs = rowSums(counts(hyperG_cds_disp))
    data_filter <<- (rs > quantile(rs, probs = 0.3))
    cdsFilt <<- hyperG_cds_disp[data_filter, ]

    # Prepare blind data
    print("Blind data prep")
    hyperG_cds_blind <<- estimateDispersions(hyperG_cds_size, method="blind")
    hyperG_vsd <<- varianceStabilizingTransformation(hyperG_cds_blind)
}

binomDist <- function() {
    # Graph binomial distribution
    print("calculating binomial distributions")
    hyperG_binom <<- nbinomTest( hyperG_cds_disp, "1G","3G" )
    write.csv(hyperG_binom[ order(hyperG_binom$pval), ], file="hyperG_binom.csv" )

    print("Creating p-value histogram and MA plot")
    hist(hyperG_binom$pval, breaks=100, col="skyblue", border="slateblue", main="p-values, unfiltered")
    plotMA(hyperG_binom, xlab = "mean of unfiltered normalized counts")

    # List significant transcripts
    print("Outputting significantly different transcripts")
    hyperG_sig = hyperG_binom[hyperG_binom$padj < 0.1,]
    write.csv( hyperG_sig[ order(hyperG_sig$pval), ], file="hyperG_sig_transcripts_0.1.csv" )
}

deseqPCA <- function() {
    # Principle Component Analysis plot
    print("DESeq PCA unfiltered")
    print(plotPCA(hyperG_vsd, intgroup = 'condition', ntop = 500))
}

qualControl <- function(){
    print("Quality control: plotting dispersion")
    plotDispEsts(hyperG_cds_disp, xlab = 'dispersion of unfiltered normalized counts')
}

filterBinomDist <- function(){
    # Graph binomial distribution
    print("Binomial distributions of filtered data")
    filter_binom <<- nbinomTest( cdsFilt, "1G","3G" )
    write.csv(filter_binom[ order(filter_binom$pval), ], file="filter_binom.csv" )

    print("Creating p-value histogram and MA plot")
    hist(filter_binom$pval, breaks=100, col="skyblue", border="slateblue", main="p-value of filtered counts")
    plotMA(filter_binom, xlab = 'mean of filtered normalized counts')

    # List significant transcripts
    print("Significantly different transcripts in filtered normalized counts")
    filter_sig <<- filter_binom[filter_binom$padj < 0.1,]
    write.csv( filter_sig[ order(filter_sig$pval), ], file="filter_significant_transcripts_0.1.csv" )
}

filterQualControl <- function(){
	# Compare p-values of filtered data
	print("p-values of filtered data")
	h1 = hist(hyperG_binom$pval, breaks=50, plot=FALSE)
	h2 = hist(filter_binom$pval, breaks=50, plot=FALSE)
	colori = c('do not pass'="khaki",'pass'="powderblue")
	barplot(height = rbind(h1$counts, h2$counts), beside = FALSE, col = colori, space = 0, main = "p values before and after filtering", ylab="frequency")
	text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)), adj = c(0.5,1.7), xpd=NA)
	legend("topright", fill=rev(colori), legend=rev(names(colori)))
}

heatmap <- function(){
    # make heatmap of transformed counts to look for patterns
    print("building heatmap")
    select = order(rowMeans(counts(hyperG_cds_size)), decreasing=TRUE)[1:30]
    hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
    heatmap.2(exprs(hyperG_vsd)[select,], col = hmcol, trace="none",margin=c(10,6))
}

similarityPlot <- function() {
    # Display similarity heatmap
    print("Similarity matrix")
    dists = dist( t( exprs(hyperG_vsd) ) )
    mat = as.matrix( dists )
    hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
    rownames(mat) = colnames(mat) = with(pData(hyperG_cds_blind), paste(condition))
    heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
}

outputSession <- function(){
    print("outputting session info")
    writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
}


dataPrep("transcript_TPM100.txt")
binomDist()
deseqPCA()
qualControl()
filterBinomDist()
filterQualControl()
heatmap()
similarityPlot()
outputSession()
