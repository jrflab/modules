#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("TitanCNA"))
suppressPackageStartupMessages(library("rtracklayer"))
suppressPackageStartupMessages(library("doMC"))
suppressPackageStartupMessages(library("hwriter"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("stringr"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList <- list(
        make_option("--outPrefix", default = NULL, type = "character", action = "store", help ="output prefix (required)"),
        make_option("--plotPrefix", default = NULL, type = "character", action = "store", help ="plot output prefix (required)"),
        make_option("--gcWig", default = NULL, type = "character", action = "store", help ="GC wig (required)"),
        make_option("--mapWig", default = NULL, type = "character", action = "store", help ="mappability wig (required)"),
        make_option("--numCores", default = 1, type = "integer", action = "store", help ="number of cores [default = %default]"),
        make_option("--numClusters", default = 5, type = "integer", action = "store", help ="number of clusters [default = %default]"),
        make_option("--ploidyPrior", default = 2, type = "double", action = "store", help ="ploidy prior [default = %default]"),
        make_option("--txnExpLen", default = 1e10, type = "double", action = "store", help ="self-transition probability [default = %default]"),
        make_option("--txnZstrength", default = 5e5, type = "double", action = "store", help ="clonal-cluster transition probability [default = %default]"),
        make_option("--tumorWig", default = NULL, type = "character", action = "store", help ="tumor wig (required)"),
        make_option("--normalWig", default = NULL, type = "character", action = "store", help ="normal wig (required)"),
        make_option("--genomeStyle", default = "NCBI", type = "character", action = "store", help ="genome style: NCBI (no chr) or UCSC (include chr) [default %default]"),
        make_option("--includeY", default = F, action = "store_true", help ="include Y chromosome"),
        make_option("--targetBed", default = NULL, type = "character", action = "store", help ="targeted interval bed"))

parser <- OptionParser(usage = "%prog [options] [tumour allele count file]", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need tumour allele count file\n\n")
    print_help(parser);
    stop();
} else if (is.null(opt$outPrefix)) {
    cat("Need output prefix\n\n")
    print_help(parser);
    stop();
} else if (is.null(opt$plotPrefix)) {
    cat("Need plot output prefix\n\n")
    print_help(parser);
    stop();
} else if (is.null(opt$gcWig)) {
    cat("Need gc wig file\n\n")
    print_help(parser);
    stop();
} else if (is.null(opt$mapWig)) {
    cat("Need mappability wig file\n\n")
    print_help(parser);
    stop();
} else if (is.null(opt$tumorWig)) {
    cat("Need tumour wig file\n\n")
    print_help(parser);
    stop();
} else if (is.null(opt$normalWig)) {
    cat("Need normal wig file\n\n")
    print_help(parser);
    stop();
}

# fix this function so that it works on data frames
filterData <- function (data, chrs = NULL, minDepth = 10, maxDepth = 200, positionList = NULL,
                        map = NULL, mapThres = 0.9, centromeres = NULL, centromere.flankLength = 0) {
    if (!is.null(map)) {
        keepMap <- map >= mapThres
    }
    else {
        keepMap <- !logical(length = length(data$refOriginal))
    }
    if (!is.null(positionList)) {
        chrPosnList <- paste(positionList[, 1], positionList[,
                             2], sep = ":")
        chrPosnData <- paste(data$chr, data$posn, sep = ":")
        keepPosn <- is.element(chrPosnData, chrPosnList)
    }
    else {
        keepPosn <- !logical(length = length(data$chr))
    }
    keepTumDepth <- data$tumDepth <= maxDepth & data$tumDepth >=
        minDepth
    if (is.null(chrs)) {
        keepChrs <- logical(length = length(data$chr))
    }
    else {
        keepChrs <- is.element(data$chr, chrs)
    }
    cI <- keepChrs & keepTumDepth & !is.na(data$logR) & keepMap & keepPosn
    data <- data[cI, ]
    if (!is.null(centromeres)) {
        colnames(centromeres)[1:3] <- c("space", "start", "ends")
        data <- removeCentromere(data, centromeres, flankLength = centromere.flankLength)
    }
    return(data)
}


#options(cores = opt$numCores)
registerDoMC(opt$numCores)

#pg <- openPage(paste(opt$outPrefix, '_titan_report_', opt$numClusters, '.html', sep = ''), title = 'TITAN Plots')

fn <- arguments$args[1]
Data <- data.frame(loadAlleleCounts(fn, header = F, genomeStyle = opt$genomeStyle), stringsAsFactors = F)
Data %<>% group_by(chr) %>% filter(n() > 1) %>% ungroup
chroms <- unique(Data$chr)
params <- loadDefaultParameters(copyNumber=5, numberClonalClusters=opt$numClusters, symmetric=TRUE, data = Data)
params$ploidyParams$phi_0 <- opt$ploidyPrior

if (!is.null(opt$targetBed)) {
    targetGr <- import(opt$targetBed)
    targets <- as.data.frame(targetGr)[,1:3]
    colnames(targets) <- c("chr", "start", "stop")
    cnData <- correctReadDepth(opt$tumorWig, opt$normalWig, opt$gcWig, opt$mapWig, targetedSequence = targets, genomeStyle = opt$genomeStyle)
} else {
    cnData <- correctReadDepth(opt$tumorWig, opt$normalWig, opt$gcWig, opt$mapWig, genomeStyle = opt$genomeStyle)
}

logR <- getPositionOverlap(Data$chr, Data$posn, cnData)
Data$logR <- log(2^logR)
rm(logR, cnData)

Data <- filterData(Data, chroms, minDepth = 10, maxDepth = 250)
mScore <- as.data.frame(wigToRangedData(opt$mapWig))
mScore <- getPositionOverlap(Data$chr, Data$posn, mScore[,-4])
Data <- filterData(Data, chroms, minDepth = 10, maxDepth = 250, map = mScore, mapThres = 0.8)
Data %<>% group_by(chr) %>% filter(n() > 1) %>% ungroup
chroms <- unique(Data$chr)

convergeParams <- runEMclonalCN(Data, gParams=params$genotypeParams, nParams=params$normalParams,
                                pParams=params$ploidyParams, sParams=params$cellPrevParams,
                                maxiter=20, maxiterUpdate=1500, txnExpLen=opt$txnExpLen,
                                txnZstrength=opt$txnZstrength, useOutlierState=FALSE,
                                normalEstimateMethod="map", estimateS=TRUE, estimatePloidy=TRUE)


optimalPath <- viterbiClonalCN(Data, convergeParams)

fn <- paste(opt$outPrefix, '.titan.txt', sep = "")
if (opt$numClusters <= 2) {
    results <- outputTitanResults(Data, convergeParams, optimalPath, filename = fn, posteriorProbs = F, subcloneProfiles = T)
} else {
    results <- outputTitanResults(Data, convergeParams, optimalPath, filename = fn, posteriorProbs = F)
}

fn <- paste(opt$outPrefix, '.params.txt', sep = "")
outputModelParameters(convergeParams, results, fn)

# plots
norm <- convergeParams$n[length(convergeParams$n)]
ploidy <- convergeParams$phi[length(convergeParams$phi)]

#library(SNPchip)  ## use this library to plot chromosome idiogram (optional)
outplot <- paste(opt$plotPrefix, '.titan.png', sep = '')
png(outplot,width=1200,height=1000,res=100, type = 'cairo-png')
if (opt$numClusters <= 2) { 
    par(mfrow=c(4,1))
} else {
    par(mfrow=c(3,1))
}
plotCNlogRByChr(results, chr = NULL, ploidy=ploidy, geneAnnot=NULL, spacing=4,ylim=c(-2,2),cex=0.5)
plotAllelicRatio(results, chr = NULL, geneAnnot=NULL, spacing=4, ylim=c(0,1),cex=0.5)
plotClonalFrequency(results, chr = NULL, normal=tail(convergeParams$n,1), geneAnnot=NULL, spacing=4,ylim=c(0,1),cex=0.5)
if (opt$numClusters <= 2){ 
    plotSubcloneProfiles(results, chr = NULL, cex = 2, spacing=6)
}
#pI <- plotIdiogram(chr,build="hg19",unit="bp",label.y=-4.25,new=FALSE,ylim=c(-2,-1))
null <- dev.off()

for (chr in intersect(results$Chr, chroms)) {
    outplot <- paste(opt$plotPrefix, '.titan.', chr, ".png", sep = '')
    #hwriteImage(basename(outplot), pg, br = T)
    png(outplot,width=1200,height=1000,res=100, type = 'cairo-png')
    if (opt$numClusters <= 2) { 
        par(mfrow=c(4,1))
    } else {
        par(mfrow=c(3,1))
    }
    plotCNlogRByChr(results, chr, ploidy=ploidy, geneAnnot=NULL, spacing=4,ylim=c(-2,2),cex=0.5,main=  chr)
    plotAllelicRatio(results, chr, geneAnnot=NULL, spacing=4, ylim=c(0,1),cex=0.5,main= chr)
    plotClonalFrequency(results, chr, normal=tail(convergeParams$n,1), geneAnnot=NULL, spacing=4,ylim=c(0,1),cex=0.5,main= chr)
    if (opt$numClusters <= 2){ 
        plotSubcloneProfiles(results, chr, cex = 2, spacing=6, main=chr)
    }
    #pI <- plotIdiogram(chr,build="hg19",unit="bp",label.y=-4.25,new=FALSE,ylim=c(-2,-1))
    null <- dev.off()
}

#closePage(pg)

fn <- paste(opt$outPrefix, '.titan.Rdata', sep = '')
save(Data, results, convergeParams, optimalPath, file = fn)

