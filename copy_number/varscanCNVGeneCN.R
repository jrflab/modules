#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("facets"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("Cairo"))
suppressPackageStartupMessages(library("RMySQL"))
suppressPackageStartupMessages(library("rtracklayer"))
suppressPackageStartupMessages(library("CGHbase"))
suppressPackageStartupMessages(library("Biobase"))

optList <- list(
                make_option("--outFile", default = NULL, help = "output file"),
				make_option("--mysqlHost", default = '10.0.200.48', help = "MySQL server hostname"),
				make_option("--mysqlPort", default = 38493, help = "MySQL server port"),
				make_option("--mysqlUser", default = 'embl', help = "MySQL server username"),
				make_option("--mysqlPassword", default = NULL, help = "MySQL server password"),
				make_option("--mysqlDb", default = 'homo_sapiens_core_75_37', help = "MySQL server database"),
                make_option("--genesFile", default = NULL, help = "list of genes to include (hgnc symbols)"))
parser <- OptionParser(usage = "%prog [options] [facets files]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need varscanCNV output files\n")
    print_help(parser);
    stop();
} else if (is.null(opt$outFile)) {
    cat("Need output prefix\n")
    print_help(parser);
    stop();
} else {
    varscanCNVFiles <- arguments$args
}

connect <- function() dbConnect(MySQL(), host = opt$mysqlHost, port = opt$mysqlPort, user = opt$mysqlUser, password = opt$mysqlPassword, dbname = opt$mysqlDb)
cat('Connecting to ensembl ... ')
mydb <- connect()
on.exit(dbDisconnect(mydb))

query <- "select r.name as chrom,
g.seq_region_start as start,
g.seq_region_end as end,
x.display_label as hgnc,
k.band as band
from gene as g
join seq_region as r on g.seq_region_id = r.seq_region_id
join xref as x on g.display_xref_id = x.xref_id
left join karyotype k on g.seq_region_id = k.seq_region_id
and ((g.seq_region_start >= k.seq_region_start and g.seq_region_start <= k.seq_region_end)
or (g.seq_region_end >= k.seq_region_start and g.seq_region_end <= k.seq_region_end))
where x.external_db_id = 1100;"
repeat {
    rs <- try(dbSendQuery(mydb, query), silent = T)
    if (is(rs, "try-error")) {
        cat("Lost connection to mysql db ... ")
        mydb <- connect()
        cat("reconnected\n")
    } else {
        break
    }
}
genes <- dbFetch(rs, -1)
cat(paste("Found", nrow(genes), "records\n"))

genes %<>% filter(chrom %in% as.character(c(1:22, "X", "Y"))) %>%
    distinct(hgnc) %>%
    arrange(as.integer(chrom), start, end)

if (!is.null(opt$genesFile)) {
    g <- scan(opt$genesFile, what = 'character')
    genes %<>% filter(hgnc %in% g)
    absentGenes <- g[!g %in% genes$hgnc]
    if (length(absentGenes) > 0) {
        print("Unable to find", length(absentGenes), "in database\n");
        cat(absentGenes, sep = '\n');
    }
}

cat(paste("Filtering to", nrow(genes), "records\n"))
genesGR <- genes %$% GRanges(seqnames = chrom, ranges = IRanges(start, end), band = band, hgnc = hgnc)
save(genesGR, file="temp.Rdata")            
mm <- lapply(varscanCNVFiles, function(f) {
    tab <- read.delim(f, as.is=T)
    tab$Chromosome[which(tab$Chromosome==23)] <- "X"

    tabGR <- tab %$% GRanges(seqnames = Chromosome, ranges = IRanges(Start, End))
    mcols(tabGR) <- tab %>% select(nBins, log2Ratio)

    fo <- findOverlaps(tabGR, genesGR)

    df <- as.data.frame(cbind(mcols(genesGR)[subjectHits(fo),], mcols(tabGR)[queryHits(fo),]))
	print(head(df))
    df %<>% group_by(hgnc) %>% top_n(1, abs(log2Ratio))

    load(gsub("collapsed_seg.txt", "segment.Rdata", f, fixed=T))
    noise <- median(abs(assayDataElement(segmented,"copynumber")-assayDataElement(segmented,"segmented")))

    lrr <- sort(assayDataElement(segmented,"copynumber"))
    if (noise <= 0.2) { lrr <- lrr[round(0.25*length(lrr)):round(0.75*length(lrr))]
    } else if ( noise <= 0.3 ) { lrr <- lrr[round(0.275*length(lrr)):round(0.725*length(lrr))]
    } else { lrr <- lrr[round(0.3*length(lrr)):round(0.7*length(lrr))]}

    df$GL2 <- 0
    df$GL2[df$log2Ratio < median(lrr)-(2.5*sd(lrr))] <- -1
    df$GL2[df$log2Ratio < median(lrr)-(7*sd(lrr))] <- -2
    df$GL2[df$log2Ratio > median(lrr)+(2*sd(lrr))] <- 1
    df$GL2[df$log2Ratio > median(lrr)+(6*sd(lrr))] <- 2

    df %>% select(hgnc, GL2)
})
names(mm) <- varscanCNVFiles
for (f in varscanCNVFiles) {
    n <- sub('\\..*', '', sub('.*/', '', f))
    colnames(mm[[f]])[2] <- paste(n, c("LRR_threshold"), sep="_")
}

mm <- left_join(genes, join_all(mm, type = 'full', by="hgnc")) %>% arrange(as.integer(chrom), start, end)
write.table(mm, file=opt$outFile, sep="\t", row.names=F, na="", quote=F)


