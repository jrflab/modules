#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("RMySQL"))
suppressPackageStartupMessages(library("biomaRt"));
suppressPackageStartupMessages(library("plyr"));
suppressPackageStartupMessages(library("dplyr"));
suppressPackageStartupMessages(library("tidyr"));
suppressPackageStartupMessages(library("stringr"));
suppressPackageStartupMessages(library("magrittr"));


if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }));
}
options(useFancyQuotes = F)

optList <- list(
                make_option("--ref", default = 'b37', help = "reference [default %default]"),
                make_option("--breakpointsFile", default = NULL, help = "integrate breakpoints file"),
                make_option("--sumFile", default = NULL, help = "integrate sum file"),
                make_option("--exonsFile", default = NULL, help = "integrate exons file"),
                make_option("--java", default = 'java', help = "java binary"),
                make_option("--mysqlHost", default = '10.0.200.48', help = "MySQL server hostname"),
                make_option("--mysqlPort", default = 38493, help = "MySQL server port"),
                make_option("--mysqlUser", default = 'embl', help = "MySQL server username"),
                make_option("--mysqlPassword", default = 'embl', help = "MySQL server password"),
                make_option("--mysqlDb", default = 'homo_sapiens_core_75_37', help = "MySQL server database"),
                make_option("--oncofuseJar", default = '~/share/usr/oncofuse-v1.0.6/Oncofuse.jar', help = "oncofuse jar"),
                make_option("--oncofuseTissueType", default = 'EPI', help = "oncofuse tissue type"),
                make_option("--outPrefix", default = NULL, help = "Output prefix"));

parser <- OptionParser(usage = "%prog [options]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;


if (is.null(opt$breakpointsFile)) {
    cat("Need breakpoints file\n");
    print_help(parser);
    stop();
} else if (is.null(opt$sumFile)) {
    cat("Need summary file\n");
    print_help(parser);
    stop();
} else if (is.null(opt$exonsFile)) {
    cat("Need exons file\n");
    print_help(parser);
    stop();
} else if (is.null(opt$outPrefix)) {
    cat("Need output prefix\n");
    print_help(parser);
    stop();
} 

breakpoints <- read.delim(opt$breakpointsFile, as.is = T, check.names = F)
summ <- read.delim(opt$sumFile, as.is = T, check.names = F)
exons <- read.delim(opt$exonsFile, as.is = T, check.names = F)
if (nrow(breakpoints) == 0 || nrow(summ) == 0 || nrow(exons) == 0) {
    fn <- str_c(opt$outPrefix, ".oncofuse.txt")
    x <- unlist(strsplit('x3p x5p chr1 rna_bk1 exon_bk1 chr2 rna_bk2 exon_bk2 wgs_bk1 wgs_bk2 fusion_candidate reciprocal tier type en_rna sp_rna splicings id x5p_transcript x5p_exon x5p_exon_strand x5p_exon_chr x5p_exon_start x5p_exon_end x5p_exon_seq x5p_exon_150 x3p_transcript x3p_exon x3p_exon_strand x3p_exon_chr x3p_exon_start x3p_exon_end x3p_exon_seq x3p_exon_150 x5p_strand x5p_start x5p_end x3p_strand x3p_start x3p_end FUSION_ID TISSUE SPANNING_READS ENCOMPASSING_READS X5_FPG_GENE_NAME X5_IN_CDS. X5_SEGMENT_TYPE X5_SEGMENT_ID X5_COORD_IN_SEGMENT X5_FULL_AA X5_FRAME X3_FPG_GENE_NAME X3_IN_CDS. X3_SEGMENT_TYPE X3_SEGMENT_ID X3_COORD_IN_SEGMENT X3_FULL_AA X3_FRAME FPG_FRAME_DIFFERENCE P_VAL_CORR DRIVER_PROB EXPRESSION_GAIN X5_DOMAINS_RETAINED X3_DOMAINS_RETAINED X5_DOMAINS_BROKEN X3_DOMAINS_BROKEN X5_PII_RETAINED X3_PII_RETAINED CTF G H K P TF', ' '))
    fc <- file(fn)
    writeLines(paste(x, collapse = '\t'), fc)
    close(fc)
    quit('no', 0)
}

colnames(breakpoints)[1:2] <- c('x5p','x3p')
colnames(breakpoints) <- tolower(colnames(breakpoints))
colnames(summ)[2:3] <- c('x5p','x3p')
colnames(summ) <- tolower(colnames(summ))
exons %<>%
    setNames(c("id", "x5p", "x3p",
               "x5p_transcript", "x5p_exon", "x5p_exon_strand", "x5p_exon_chr", "x5p_exon_start", "x5p_exon_end",
               "x5p_exon_seq", "x5p_exon_150",
               "x3p_transcript", "x3p_exon", "x3p_exon_strand", "x3p_exon_chr", "x3p_exon_start", "x3p_exon_end",
               "x3p_exon_seq", "x3p_exon_150"))

# split genes separated by '/' to create separate rows (with same info)
X <- cbind(breakpoints, summ[, -(2:3)])
X <- left_join(X, exons)
y <- strsplit(X[['x5p']], '/', fixed = T)
X <- data.frame(x5p = unlist(y), X[rep(1:nrow(X), sapply(y, length)), -1], stringsAsFactors = F)
y <- strsplit(X[['x3p']], '/', fixed = T)
X <- data.frame(x3p = unlist(y), X[rep(1:nrow(X), sapply(y, length)), -2], stringsAsFactors = F)
results <- X
results$chr1 <- as.character(results$chr1)
results$chr2 <- as.character(results$chr2)

if (opt$ref == "b37" || opt$ref == "hg19") {
    dbName <- 'homo_sapiens_core_75_37'
} else if (opt$ref == "mm10" || opt$ref == "GRCm38") {
    dbName <- 'mus_musculus_core_82_38'
} else {
    cat("reference unsupported:", opt$ref, '\n')
    quit('no', status = 1)
}
cat(opt$ref, ': Connecting to ensembl', dbName, '... ')
connect <- function() dbConnect(MySQL(), host = opt$mysqlHost, port = opt$mysqlPort, user = opt$mysqlUser, password = opt$mysqlPassword, dbname = opt$mysqlDb)
mydb <- connect()
on.exit(dbDisconnect(mydb))
cat('done\n')

syms <- paste(sQuote(unique(c(results$x5p, results$x3p))), collapse = ',')
query1 <- paste("SELECT DISTINCT G.seq_region_strand AS strand, ES.synonym AS gene_symbol, SR.name AS chrom, G.seq_region_start AS start, G.seq_region_end AS end
               FROM gene G
               JOIN seq_region SR ON (SR.seq_region_id = G.seq_region_id)
               JOIN xref X ON (G.display_xref_id = X.xref_id)
               LEFT JOIN external_synonym ES using (xref_id)
               WHERE ES.synonym in (", syms, ");")
query2 <- paste("SELECT DISTINCT G.seq_region_strand AS strand, X.display_label AS gene_symbol, SR.name AS chrom, G.seq_region_start AS start, G.seq_region_end AS end
               FROM gene G
               JOIN seq_region SR ON (SR.seq_region_id = G.seq_region_id)
               JOIN xref X ON (G.display_xref_id = X.xref_id)
               LEFT JOIN external_synonym ES using (xref_id)
               WHERE X.display_label in (", syms, ");")
cat("Looking up strand for genes ...\n")
sendQuery <- function(query) {
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
    fetch(rs, -1)
}
queryResults <- unique(rbind(sendQuery(query1), sendQuery(query2)))
queryResults %<>% group_by(chrom, strand, gene_symbol) %>% summarize(start = min(start), end = max(end)) %>% ungroup
cat(paste("Found", nrow(queryResults), "records\n"))

results %<>% left_join(queryResults, by = c("x5p" = "gene_symbol", "chr1" = "chrom")) %>% rename(x5p_strand = strand, x5p_start = start, x5p_end = end)
results %<>% mutate(x5p_strand = ifelse(pmax(x5p_start, x5p_end) - rna_bk1 > 0 & rna_bk1 - pmin(x5p_start, x5p_end) > 0, x5p_strand, NA))
results %<>% left_join(queryResults, by = c("x3p" = "gene_symbol", "chr2" = "chrom")) %>% rename(x3p_strand = strand, x3p_start = start, x3p_end = end)
results %<>% mutate(x3p_strand = ifelse(pmax(x3p_start, x3p_end) - rna_bk2 > 0 & rna_bk2 - pmin(x3p_start, x3p_end) > 0, x3p_strand, NA))
results %<>% mutate(x3p_strand = ifelse(x3p_strand == 1, '+', '-'))
results %<>% mutate(x5p_strand = ifelse(x5p_strand == 1, '+', '-'))
results %<>% mutate(x5p_exon_strand = ifelse(is.na(x5p_exon_strand),
                                             x5p_strand,
                                             x5p_exon_strand))
results %<>% mutate(x3p_exon_strand = ifelse(is.na(x3p_exon_strand),
                                             x3p_strand,
                                             x3p_exon_strand))
oncofuse <- results %$% data.frame(CHROM5p=str_c("chr", chr1), RNA_BK1 = rna_bk1, CHROM3p = str_c("chr", chr2),
                                   RNA_BK2 = rna_bk2, TISSUE_TYPE = opt$oncofuseTissueType, STRAND5p = x5p_exon_strand, STRAND3p = x3p_exon_strand)
x <- rowSums(is.na(oncofuse)) == 0
oncofuse <- oncofuse[x, ]
rmResults <- results[!x, ]
results <- results[x, ]
oncofuse %<>% mutate(RNA_BK1 = ifelse(STRAND5p == "+", RNA_BK1 + 1, RNA_BK1 - 1))
oncofuse %<>% mutate(RNA_BK2 = ifelse(STRAND3p == "+", RNA_BK2 - 1, RNA_BK2 + 1))

results$id <- paste("chr", results$chr1, ":", oncofuse$RNA_BK1, ">",
                    "chr", results$chr2, ":", oncofuse$RNA_BK2, sep = '')
rmResults$id <- NA

ifn <- str_c(opt$outPrefix, ".oncofuse.input.txt")
write.table(oncofuse[which(!is.na(oncofuse[,"STRAND5p"])),1:5], file=ifn, sep="\t", row.names=F, col.names=F, quote=F, na="")

ofn <- str_c(opt$outPrefix, ".oncofuse.output.txt")
cmd <- paste(opt$java, '-Xmx1G -jar', opt$oncofuseJar, ifn, "coord - ", ofn)
system(cmd, wait = T)

oncofuse_output <- read.delim(ofn, as.is=T)
results <- left_join(results, oncofuse_output, by = c("id" = "GENOMIC"))
results <- select(results, -SAMPLE_ID)

# add back in removed result rows
rmResults[, colnames(results)[!colnames(results) %in% colnames(rmResults)]] <- NA
results <- rbind(results, rmResults)

fn <- str_c(opt$outPrefix, ".oncofuse.txt")
write.table(results, file=fn, sep="\t", row.names=F, quote=F, na="")

