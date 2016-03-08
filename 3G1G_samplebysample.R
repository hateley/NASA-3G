#Comparing the expression from different chromosomes to determine
#male to female sample ratio

# source("source_if_on_network")
# source("3G1G_samplebysample.R")

library("DESeq")
library("gplots")
#library("RColorBrewer")

dataPrep <- function(counts) {
    # Create source data
    print("Creating source data")
    countsTable <<- read.delim(counts,header=TRUE,row.names=2)

    # Create metadata
    counts_design <<- data.frame(
        row.names = colnames(countsTable),
        condition = c('chr', '1G','1G','1G','1G','3G','3G','3G','3G'))

    # Prepare data
    print("Preparing data")
    cds <- newCountDataSet(countsTable[2:9], counts_design$condition[2:9])
    cds_size <- estimateSizeFactors(cds)
    sizes <- sizeFactors(cds_size)
    counts_norm <- sweep(countsTable[2:9],MARGIN=2,sizes,'/')
    counts_norm$chrom <- countsTable[,1]
    sized_counts <<- data.frame(counts_norm[9],counts_norm[2:4], counts_norm[7],counts_norm[5:6],counts_norm[1],counts_norm[8])
}

subsetData <- function() {
    # Subset data into chromosome sets
    print("subsetting data")
    chro2L <<- sized_counts[sized_counts$chrom == '2L',]
    chro2R <<- sized_counts[sized_counts$chrom == '2R',]
    chro3R <<- sized_counts[sized_counts$chrom == '3R',]
    chro3L <<- sized_counts[sized_counts$chrom == '3L',]
    chroX <<- sized_counts[sized_counts$chrom == 'X',]
    chro4 <<- sized_counts[sized_counts$chrom == '4',]
    chroY <<- sized_counts[sized_counts$chrom  == 'Y',]
    chroother <<- sized_counts[sized_counts$chrom != '2L' & sized_counts$chrom != '2R' & sized_counts$chrom != '3R' & sized_counts$chrom != '3L' & sized_counts$chrom != 'X' & sized_counts$chrom != 'Y' & sized_counts$chrom != '4',]
}

sumyvalues <- function() {
    yvals <<- c(sum(chroY$X1GR2), sum(chroY$X1GR3),sum(chroY$X1GR4),sum(chroY$switch1G),sum(chroY$X3GR1),sum(chroY$X3GR2),sum(chroY$switch3G),sum(chroY$X3GR4))
}
sumxvalues <- function() {
    xvals <<- c(sum(chroX$X1GR2), sum(chroX$X1GR3),sum(chroX$X1GR4),sum(chroX$switch1G),sum(chroX$X3GR1),sum(chroX$X3GR2),sum(chroX$switch3G),sum(chroX$X3GR4))
}

matrixscatterplot <- function(){
    pairs(~X1GR2+X1GR3+X1GR4+switch1G+X3GR1+X3GR2+switch3G+X3GR4,main="whole sample scatterplot")
}

#pairs(~log2(X1GR2)+log2(X1GR3)+log2(X1GR4)+log2(switch1G)+log2(X3GR1)+log2(X3GR2)+log2(switch3G)+log2(X3GR4),main="Y chro scatterplot")

outputSession <- function(){
    print("outputting session info")
    writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
}


dataPrep("gene_counts_loc.tsv")
subsetData()
matrixscatterplot()
#outputSession()
