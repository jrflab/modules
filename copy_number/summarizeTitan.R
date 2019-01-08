#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("TitanCNA"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("inflection"))


options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
        make_option("--outDir", default = NULL, type = "character", action = "store", help ="copy optimal files to this directory (required)"))

parser <- OptionParser(usage = "%prog [options] [sample params files]", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need sample params files\n\n")
    print_help(parser);
    stop();
} else if (is.null(opt$outDir)) {
    cat("Need output prefix dir\n\n")
    print_help(parser);
    stop();
}

fns <- arguments$args
dir.create(opt$outDir, showWarnings = F)


Data <- list()
filenames <- list()
for (fn in fns) {
    s <- sub(".*/", "", sub('\\..*', '', fn))
    i <- as.integer(str_match(fn, '\\.z([0-9])\\.params\\.txt')[,2])
    inlist <- strsplit(readLines(fn), ":\t")
    params <- lapply(inlist, tail, n = -1)
    names(params) <- lapply(inlist, head, n = 1)
    params <- lapply(params, function(x) as.numeric(unlist(strsplit(x, "[[:space:]]+"))))
    paramList <- list()
    paramList[["densBw"]] <- params[["S_Dbw dens.bw (Both)"]]
    paramList[["scat"]] <- params[["S_Dbw scat (Both)"]]
    paramList[["SDbwIndex"]] <- params[["S_Dbw validity index (Both)"]]
    paramList[["normalContaminationEstimate"]] <- params[["Normal contamination estimate"]]
    paramList[["avgTumorPloidyEstimate"]] <- params[["Average tumour ploidy estimate"]]

    Data[[s]][[i]] <- paramList
    filenames[[s]][[as.character(i)]] <- fn
}

sdbwM <- matrix(NA, nrow = length(Data), ncol = length(Data[[1]]), dimnames = list(names(Data)))
scatM <- matrix(NA, nrow = length(Data), ncol = length(Data[[1]]), dimnames = list(names(Data)))
densBwM <- matrix(NA, nrow = length(Data), ncol = length(Data[[1]]), dimnames = list(names(Data)))
avgTumorPloidyEstM <- matrix(NA, nrow = length(Data), ncol = length(Data[[1]]), dimnames = list(names(Data)))
normalContamEstM <- matrix(NA, nrow = length(Data), ncol = length(Data[[1]]), dimnames = list(names(Data)))
for (s in names(Data)) {
    for (i in 1:length(Data[[1]])) {
        sdbwM[s, i] <- Data[[s]][[i]]$SDbwIndex
        densBwM[s, i] <- Data[[s]][[i]]$densBw
        scatM[s, i] <- Data[[s]][[i]]$scat
        avgTumorPloidyEstM[s, i] <- Data[[s]][[i]]$avgTumorPloidyEstimate
        normalContamEstM[s, i] <- Data[[s]][[i]]$normalContaminationEstimate
    }
}
sdbwMin <- apply(sdbwM, 1, which.min)

maxScale <- 100
sc <- 1:maxScale
ratio <- rep(NA, length(sc))
for (i in 1:length(sc)){
    SDbw <- sc[i] * densBwM + scatM
    counts <- table(apply(SDbw, 1, which.min))
    ratio[i] <- counts[1] / sum(counts[-1])  # clust1 / sum <- C (clustC)
}
fn <- paste(opt$outDir, '/sdbw_scale.png', sep = '')
png(fn, type = 'cairo-png')
plot(sc, ratio, type="o", pch=19, xaxt="n", las=2, ylab="Clonal:Subclonal Ratio",
            xlab="Scale", main="S_Dbw = scale * dens + scat")
axis(1)

if (sum(ratio) == Inf) {
    cat("Error: Can't find inflection point\n")
    q(save = 'no')
}

## find and plot inflection plot using EDE (R package inflection)
inflectPt <- findiplist(x=as.matrix(sc), y=ratio, index=1)
abline(v=inflectPt[2,1], col="red")
mtext(text=sc[inflectPt[2,1]],side=3,at=inflectPt[2,1])
null <- dev.off()

SDbw <- sc[inflectPt[2,1]] * densBwM + scatM
optClust <- apply(SDbw, 1, which.min)

avgTumorPloidyEstOpt <- avgTumorPloidyEstM[cbind(1:nrow(sdbwM), optClust)]
normalContamEstOpt <- normalContamEstM[cbind(1:nrow(sdbwM), sdbwMin)]


results <- data.frame(row.names = names(optClust), optClust = optClust, avgTumorPloidyEst = avgTumorPloidyEstOpt, normalContamEst = normalContamEstOpt)


fn <- paste(opt$outDir, '/titan_summary.txt', sep = '')
write.table(results, file = fn, sep = '\t', quote = F)

exts <- c(".titan.png", ".titan.Rdata", ".titan.seg", ".titan.txt", ".titan_seg.txt", paste('.titan.chr', c(1:22, "X"), '.png', sep = ''))

optFns <- do.call('rbind', filenames)[cbind(names(optClust), optClust)]
optFns <- sub('\\.params\\.txt', '', optFns)
optFns <- as.vector(sapply(optFns, paste, exts, sep = ''))
null <- file.copy(optFns, opt$outDir, overwrite = T)

optFns2 <- sub('\\.z[0-9]', '', optFns)
optFns2 <- sub('.*/', paste(opt$outDir, '/', sep = ''), optFns2)
null <- file.copy(optFns, optFns2, overwrite = T)
