function (sampleName = "", targetsName = "", referenceName = "", 
    destDir = "TEQCreport", reads = get.reads(), targets = get.targets(), 
    Offset = 0, pairedend = FALSE, genome = c(NA, "hg19", "hg18"), 
    genomesize, CovUniformityPlot = FALSE, CovTargetLengthPlot = FALSE, 
    CovGCPlot = FALSE, duplicatesPlot = FALSE, baits = get.baits(), 
    WigFiles = FALSE, saveWorkspace = FALSE) 
{
    if (!file.exists(destDir)) 
        dir.create(destDir, recursive = TRUE)
    wd <- getwd()
    setwd(destDir)
    Path <- getwd()
    setwd(wd)
    print(paste("results and report are saved in folder", Path))
    imgDir <- file.path(destDir, "image")
    if (!file.exists(imgDir)) 
        dir.create(imgDir)
    if (WigFiles) {
        wigDir <- file.path(destDir, "wiggle")
        if (!file.exists(wigDir)) 
            dir.create(wigDir)
    }
    genome <- match.arg(genome)
    if (missing(genomesize) & is.na(genome)) 
        stop("either 'genome' or 'genomesize' has to be specified")
    if (CovGCPlot) {
        if (missing(baits)) 
            stop("if 'CovGCPlot = TRUE', a 'baits' table has to be specified")
        else if (data.class(baits) != "RangedData") 
            stop("the 'baits' table has to be of class 'RangedData'")
    }
    print("reading data...")
    if (missing(reads)) 
        stop("'reads' have to be specified")
    if (data.class(reads) != "RangedData") 
        stop("the 'reads' table has to be of class 'RangedData'")
    if (missing(targets)) 
        stop("'targets' have to be specified")
    if (data.class(targets) != "RangedData") 
        stop("the 'targets' table has to be of class 'RangedData'")
    n.reads <- nrow(reads)
    n.targets <- nrow(targets)
    if (pairedend) {
        print("collapsing reads to pairs...")
        readpairs <- reads2pairs(reads)
        if (is.list(readpairs)) {
            n.pairs <- nrow(readpairs$readpairs)
            n.singles <- nrow(readpairs$singleReads)
        }
        else {
            n.pairs <- nrow(readpairs)
            n.singles <- 0
        }
    }
    ft <- fraction.target(targets, Offset = Offset, genome = genome, 
        genomesize = genomesize)
    if (pairedend) {
        print("calculating fraction of on-target read pairs")
        fr <- fraction.reads.target(readpairs, targets, Offset = Offset)
    }
    else {
        print("calculating fraction of on-target reads")
        fr <- fraction.reads.target(reads, targets, Offset = Offset)
    }
    enr <- as.character(round(fr/ft))
    print("calculating coverage...")
    Coverage <- coverage.target(reads, targets, Offset = Offset)
    avgcov <- data.frame(round(Coverage$avgTargetCoverage, 2), 
        round(Coverage$targetCoverageSD, 2), matrix(Coverage$targetCoverageQuantiles, 
            ncol = 5))
    names(avgcov) <- c("avgTargetCoverage", "targetCoverageSD", 
        paste(names(Coverage$targetCoverageQuantiles), "quantile"))
    print("counting reads per target...")
    targetcov0 <- Coverage$targetCoverages
    targetcov0 <- readsPerTarget(reads, targetcov0, Offset = Offset)
    targetcov <- as.data.frame(targetcov0)
    write.table(targetcov, file = file.path(destDir, "target_coverage.txt"), 
        sep = "\t", row.names = F, quote = F)
    if (nrow(targetcov) > 20) 
        targetcov <- rbind(apply(targetcov[1:20, ], 2, as.character), 
            "...")
    sensi <- round(covered.k(Coverage$coverageTarget) * 100, 
        2)
    N <- paste(">=", names(sensi), "X", sep = "")
    sensi <- paste(sensi, "%", sep = "")
    names(sensi) <- N
    print("generating figures...")
    values <- list(SAMPLE = sampleName, NREADS = as.character(nrow(reads)), 
        TARGETS = targetsName, NTARGETS = as.character(nrow(targets)), 
        REFERENCE = referenceName, OFFSET = as.character(Offset), 
        SPECIFICITY = hwrite(paste(round(fr * 100, 2), "%", sep = "")), 
        ENRICHMENT = hwrite(enr), CHROM_BARPLOT = htmlChromBarplot(destDir, 
            reads, targets), AVGCOV = hwrite(avgcov), COVTARG = hwrite(targetcov), 
        SENSITIVITY = hwrite(sensi), COV_HIST = htmlCoverageHist(destDir, 
            Coverage$coverageTarget, covthreshold = 8))
    if (pairedend) 
        values <- c(values, list(NPAIRS = as.character(n.pairs), 
            NSINGLES = as.character(n.singles), ISIZEHIST = htmlInsertSizeHist(destDir, 
                readpairs)))
    if (CovUniformityPlot) 
        values <- c(values, list(COV_UNIFORM = htmlCovUniformity(destDir, 
            Coverage)))
    if (CovTargetLengthPlot) 
        values <- c(values, list(COV_TARGLEN = htmlCovTargetLength(destDir, 
            targetcov0)))
    if (CovGCPlot) 
        values <- c(values, list(COV_GC = htmlCovGC(destDir, 
            Coverage$coverageAll, baits)))
    if (duplicatesPlot) {
        print("duplicates analysis...")
        if (pairedend) 
            values <- c(values, list(DUPLICATES = htmlDuplicatesBarplot(destDir, 
                readpairs, targets, ylab = "Fraction of read pairs")))
        else values <- c(values, list(DUPLICATES = htmlDuplicatesBarplot(destDir, 
            reads, targets)))
    }
    if (WigFiles) {
        make.wigfiles(Coverage$coverageAll, filename = file.path(destDir, 
            "wiggle", "Coverage"))
        chroms <- names(Coverage$coverageAll)
        wignames <- paste("Coverage_", chroms, ".wig", sep = "")
        values <- c(values, list(WIG = hwrite(matrix(wignames), 
            link = file.path(".", "wiggle", wignames))))
    }
    print("generating html report...")
    make.report(destDir = destDir, values = values, pairedend = pairedend, 
        CovUniformityPlot = CovUniformityPlot, CovTargetLengthPlot = CovTargetLengthPlot, 
        CovGCPlot = CovGCPlot, duplicatesPlot = duplicatesPlot, 
        WigFiles = WigFiles)
    if (saveWorkspace) {
        print("saving workspace...")
        if (pairedend) 
            save(reads, targets, Coverage, readpairs, file = file.path(destDir, 
                "results.RData"))
        else save(reads, targets, Coverage, file = file.path(destDir, 
            "results.RData"))
    }
}
