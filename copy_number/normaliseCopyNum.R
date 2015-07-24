#!/usr/bin/env Rscript
# MSK-IMPACT copy number procedure

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("RColorBrewer"));
suppressPackageStartupMessages(library("plyr"));
suppressPackageStartupMessages(library("dplyr"));
suppressPackageStartupMessages(library("tidyr"));
suppressPackageStartupMessages(library("stringr"));
suppressPackageStartupMessages(library("magrittr"));
suppressPackageStartupMessages(library("limma"));
suppressPackageStartupMessages(library("DNAcopy"));
suppressPackageStartupMessages(library("foreach"));

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList <- list(
                make_option("--sampleSetsFile", default = NULL, help = "sample sets file"),
                make_option("--outDir", default = NULL, help = "output file directory for plots"),
                make_option("--nucFile", default = NULL, help = "bedtools nuc output for 100bp window modified target bed file"),
                make_option("--minCov", default = 5, type = 'integer', help = "minimum coverage required in a window"),
                make_option("--alpha", default = 0.05, type = "double", action = "store", help ="sginficance levels for the test to accept change-points"),
                make_option("--outlierSDscale", default = 2.5, type = "double", action = "store", help = "the number of SDs from the median in the smoothing region where a smoothed point is positioned"),
                make_option("--trim", default = 0.05, type = "double", action = "store", help = "proportion of data to be trimmed for variance calculation for smoothing outliers and undoing splits based on SD"),
                make_option("--undoSD", default = 2, type = "double", action = "store", help = "the number of SDs between means to keep a split"),
                make_option("--centromereFile", help = "Centromere position table"));

parser <- OptionParser(usage = "%prog [options] [gatk DoC interval summary files]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$nucFile)) {
    cat("Need nuc file\n")
    print_help(parser);
    stop();
} else if (is.null(opt$outDir)) {
    cat("Need output dir\n")
    print_help(parser);
    stop();
} else if (is.null(opt$sampleSetsFile)) {
    cat("Need sample sets file\n")
    print_help(parser);
    stop();
} else if (is.null(opt$centromereFile)) {
    cat("Need centromere file\n")
    print_help(parser);
    stop();
} else {
    files <- arguments$args;
}

dir.create(opt$outDir, showWarnings = F, recursive = T)

sampleSets <- scan(opt$sampleSetsFile, what = 'character', sep = '\n') %>% str_split(' ')
samplePairs <- sampleSets %>%
    ldply(function(x) data.frame(x[1:(length(x)-1)], last(x))) %>%
    setNames(c("Tumor", "Normal"))

nuc <- read.table(opt$nucFile, sep = '\t', header = T, comment.char = '', as.is = T)
colnames(nuc)[1:3] <- c("chr", "start", "end")
colnames(nuc) <- str_replace(colnames(nuc), "X\\d+_", "")
nuc %<>% mutate(start = start + 1)

doc <- files %>%
    llply(read.delim, as.is = T) %>%
    llply(select, -one_of('total_coverage', 'average_coverage')) %>%
    join_all(type = 'inner')

pos <- str_match(doc[["Target"]], "(.*):(.*)-(.*)") %>%
    data.frame(stringsAsFactors = F) %>% 
    setNames(c("Target", "chr", "start", "end")) %>%
    mutate_each(funs(as.integer), start:end)

doc <- inner_join(pos, doc)
doc <- inner_join(nuc, doc, by = c("chr", "start", "end"))

doc %<>% filter(rowSums(select(., one_of(str_c(samplePairs[['Normal']], '_mean_cvg'))) < opt$minCov) == 0)

# step 1: square-root transformed
doc %<>% mutate_each(funs(sqrt), ends_with("total_cvg"))

# step 2: Loess normalization of depth based on GC content
loessNorm <- function(x, y) x + loessFit(x, y)$residual
doc_norm <- doc %>% select(ends_with('total_cvg')) %>%
    apply(2, loessNorm, doc$pct_gc)
#bugged: doc_norm <- doc %>% mutate_each(funs(loessNorm(., pct_gc)), ends_with('total_cvg'))
colnames(doc_norm) <- gsub("_total_cvg", "_norm_cvg", colnames(doc_norm))
doc <- bind_cols(doc, as.data.frame(doc_norm))

# step 3: filter out regions with normalised depth in the top or bottom 5% in >=20% of normal samples
low <- select(doc, one_of(str_c(samplePairs[["Normal"]], '_norm_cvg'))) %>%
    (colwise(function(x) { x < quantile(x, 0.05) })) %>% 
    apply(1, function(x) sum(x) >= length(x) * 0.2)

high <- select(doc, one_of(str_c(samplePairs[["Normal"]], '_norm_cvg'))) %>%
    (colwise(function(x) { x > quantile(x, 0.95) })) %>% 
    apply(1, function(x) sum(x) >= length(x) * 0.2)

doc_ft <- doc %>% filter(!(low | high))
doc_ft %<>% mutate(chr = str_replace(chr, 'X', 23))
doc_ft %<>% mutate(chr = str_replace(chr, 'Y', 24))
doc_ft %<>% mutate(chr = as.integer(chr))

# step 4: compute log ratios
genomedat <- foreach (i = 1:nrow(samplePairs), .combine = cbind) %do% {
    tumor <- str_c(samplePairs[i, "Tumor"], "_norm_cvg")
    normal <- str_c(samplePairs[i, "Normal"], "_norm_cvg")
    log(doc_ft[, tumor] / doc_ft[, normal], base = 2) %>% scale(scale = F)
}
colnames(genomedat) <- with(samplePairs, paste(Tumor, Normal, sep = '_'))

######################### THE REST SHOULD BE THE SAME (ALMOST THE SAME) AS VARSCAN COPYNUMBER #############
# step 5: segmentation
cna <- doc_ft %$% CNA(genomedat, chr, round(start + seq_len / 2),
           sampleid = colnames(genomedat))
smoothed.cna <- smooth.CNA(cna, outlier.SD.scale = opt$outlierSDscale, trim = opt$trim)
seg <- segment(smoothed.cna, undo.SD = opt$undoSD, alpha = opt$alpha, undo.splits = "sdundo", trim = opt$trim)
#smoothed.cna <- smooth.CNA(cna, outlier.SD.scale = 1, trim = 0.01)
#seg <- segment(smoothed.cna, undo.SD = 2, alpha = 0.05, undo.splits = "sdundo")


# step 6: plot (copied from existing script)
cen <- read.table(opt$centromereFile, sep = '\t')
for (i in colnames(genomedat)) {
    fn <- str_c(opt$outDir, '/', i, '.seg_plot.png')
    png(fn, type = 'cairo-png', height=400, width=2000)
    obj <- subset(seg, sample=i)

    objdat <- obj$data[which(!is.na(obj$data[,3])),]

    plot(objdat[,3], pch=20, xlab='Position', ylab="Copy number", xaxt='n', ylim=c(min(objdat[,3]), max(objdat[,3])+0.5), main = i)
    points(unlist(apply(obj$output, 1, function(x) {rep(x[6], x[5])})), pch = 20, col = 'blue')
    abline(v=cumsum(rle(as.vector(objdat$chrom))$lengths), col="red", lty=3)

        for (j in unique(cen[,1])) {
            pos <- cen[which(cen[,1]==j)[1],3]
            index <- which(objdat$chrom==j & objdat$maploc > pos)[1]
            if (!is.na(index)) {
                abline(v=index, col="darkgrey", lty=3)
            }
            text(cumsum(rle(as.vector(objdat$chrom))$lengths)-((rle(as.vector(objdat$chrom))$lengths)/2), max(objdat[,3])+0.5-0.25)
        }

    dev.off()
}

# step 7: write data
write.table(seg$output, file = paste(opt$outDir, "/seg.txt", sep=""), row.names = F, quote = F, sep = "\t")

segSplit <- split(seg$output, seg$output$ID)
for (i in colnames(genomedat)) {
    fn <- str_c(opt$outDir, '/', i, '.seg.txt')
    write.table(segSplit[[i]][, -1], file = fn, row.names = F, quote = F, sep = "\t")
}

