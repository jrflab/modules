#!/usr/bin/env Rscript
#### turn segmented copy number data to gene-based copy number with findOverlaps
## define HomDel as TCN=0, loss as TCN<ploidy, gain as TCN>ploidy, amp as TCN>=ploidy+4
## where ploidy= mode of TCN
### some variant of the below, also need one for the breast panel, IMPACT310 and exome

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("RColorBrewer"));
suppressPackageStartupMessages(library("GenomicRanges"));
suppressPackageStartupMessages(library("plyr"));
suppressPackageStartupMessages(library("dplyr"));
suppressPackageStartupMessages(library("tidyr"));
suppressPackageStartupMessages(library("stringr"));
suppressPackageStartupMessages(library("magrittr"));
suppressPackageStartupMessages(library("facets"));
suppressPackageStartupMessages(library("foreach"));
suppressPackageStartupMessages(library("Cairo"));
suppressPackageStartupMessages(library("RMySQL"))
suppressPackageStartupMessages(library("rtracklayer"))

optList <- list(
                make_option("--outFile", default = NULL, help = "output file"),
                make_option("--genesFile", default = NULL, help = "list of genes to include (hgnc symbols)"))
parser <- OptionParser(usage = "%prog [options] [facets files]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
    cat("Need facets output files\n")
    print_help(parser);
    stop();
} else if (is.null(opt$outFile)) {
    cat("Need output prefix\n")
    print_help(parser);
    stop();
} else {
    facetsFiles <- arguments$args
}

connect <- function() dbConnect(MySQL(), host = "10.0.200.48", port = 38493, user = "embl", password = "embl", dbname = 'homo_sapiens_core_78_38')
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
join karyotype k on g.seq_region_id = k.seq_region_id and g.seq_region_start >= k.seq_region_start and g.seq_region_end <= k.seq_region_end
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
        X <- data.frame(hgnc = absentGenes, stringsAsFactors = F)
        genes %<>% full_join(X)
    }
}

cat(paste("Filtering to", nrow(genes), "records\n"))

genesGR <- genes %$% GRanges(seqnames = chrom, ranges = IRanges(start, end), band = band, hgnc = hgnc)
            
mm <- lapply(facetsFiles, function(f) {
    tab <- read.delim(f, as.is=T)
    tab$chrom[which(tab$chrom==23)] <- "X"

    tabGR <- tab %$% GRanges(seqnames = chrom, ranges = IRanges(loc.start, loc.end))
    mcols(tabGR) <- tab %>% select(cnlr.median:lcn.em)

    fo <- findOverlaps(tabGR, genesGR)

    df <- as.data.frame(cbind(mcols(genesGR)[subjectHits(fo),], mcols(tabGR)[queryHits(fo),]))
    df %<>% group_by(hgnc) %>% top_n(1, abs(cnlr.median))

    ploidy <- table(df$tcn)
    ploidy <- as.numeric(names(ploidy)[which.max(ploidy)])

    df$GL <- 0
    df$GL[df$tcn < ploidy] <- -1
    df$GL[df$tcn == 0] <- -2
    df$GL[df$tcn > ploidy] <- 1
    df$GL[df$tcn >= ploidy + 4] <- 2
    df %>% select(hgnc, GL)
})
names(mm) <- facetsFiles
for (f in facetsFiles) {
    n <- sub('\\..*', '', sub('.*/', '', f))
    colnames(mm[[f]])[2] <- n
}

mm <- inner_join(genes, join_all(mm)) %>% arrange(as.integer(chrom), start, end)
write.table(mm, file=opt$outFile, sep="\t", row.names=F, na="", quote=F)


