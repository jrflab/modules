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
                make_option("--breakpointsFile", default = NULL, help = "integrate breakpoints file"),
                make_option("--sumFile", default = NULL, help = "integrate sum file"),
                make_option("--exonsFile", default = NULL, help = "integrate exons file"),
                make_option("--java", default = 'java', help = "java binary"),
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

breakpoints <- read.delim(opt$breakpointsFile, as.is = T, check.names = F) %>%
    setNames(c("x5p", "x3p", "chr1", "rna_bk1", "exon_bk1", "chr2", "rna_bk2", "exon_bk2", "wgs_bk1", "wgs_bk2"))
summ <- read.delim(opt$sumFile, as.is = T, check.names = F) %>%
    setNames(c("id", "x5p", "x3p", "reciprocal", "tier", "type", "en_rna", "sp_rna", "splicings"))
exons <- read.delim(opt$exonsFile, as.is = T, check.names = F) %>%
    setNames(c("id", "x5p", "x3p",
               "x5p_transcript", "x5p_exon", "x5p_exon_strand", "x5p_exon_chr", "x5p_exon_start", "x5p_exon_end",
               "x5p_exon_seq", "x5p_exon_150",
               "x3p_transcript", "x3p_exon", "x3p_exon_strand", "x3p_exon_chr", "x3p_exon_start", "x3p_exon_end",
               "x3p_exon_seq", "x3p_exon_150"))

# split genes separated by '/' to create separate rows (with same info)
X <- inner_join(breakpoints, summ)
y <- strsplit(X[['x5p']], '/', fixed = T)
X <- data.frame(x5p = unlist(y), X[rep(1:nrow(X), sapply(y, length)), -1], stringsAsFactors = F)
y <- strsplit(X[['x3p']], '/', fixed = T)
X <- data.frame(x3p = unlist(y), X[rep(1:nrow(X), sapply(y, length)), -2], stringsAsFactors = F)
results <- X %>% full_join(exons)

cat('Connecting to ensembl ... ')
connect <- function() dbConnect(MySQL(), host = "10.0.200.48", port = 38493, user = "embl", password = "embl", dbname = 'homo_sapiens_core_75_37')
mydb <- connect()
on.exit(dbDisconnect(mydb))
cat('done\n')

syms <- paste(sQuote(unique(c(results$x5p, results$x3p))), collapse = ',')
query1 <- paste("SELECT DISTINCT G.seq_region_strand AS strand, ES.synonym AS gene_symbol
               FROM gene G
               JOIN xref X ON (G.display_xref_id = X.xref_id)
               LEFT JOIN external_synonym ES using (xref_id)
               WHERE ES.synonym in (", syms, ");")
query2 <- paste("SELECT DISTINCT G.seq_region_strand AS strand, X.display_label as gene_symbol
               FROM gene G
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
cat(paste("Found", nrow(queryResults), "records\n"))

results %<>% left_join(queryResults, by = c("x5p" = "gene_symbol")) %>% rename(x5p_strand = strand)
results %<>% left_join(queryResults, by = c("x3p" = "gene_symbol")) %>% rename(x3p_strand = strand)
results %<>% mutate(x5p_exon_strand = ifelse(is.na(x5p_exon_strand),
                                             ifelse(x5p_strand == 1, '+', '-'),
                                             x5p_exon_strand))
results %<>% mutate(x3p_exon_strand = ifelse(is.na(x3p_exon_strand),
                                             ifelse(x3p_strand == 1, '+', '-'),
                                             x3p_exon_strand))
oncofuse <- results %$% data.frame(CHROM5p=str_c("chr", chr1), RNA_BK1 = rna_bk1, CHROM3p = str_c("chr", chr2),
                                   RNA_BK2 = rna_bk2, TISSUE_TYPE = opt$oncofuseTissueType, STRAND5p = x5p_exon_strand, STRAND3p = x3p_exon_strand)
x <- rowSums(is.na(oncofuse)) == 0
oncofuse <- oncofuse[x, ]
results <- results[x, ]
oncofuse %<>% mutate(RNA_BK1 = ifelse(STRAND5p == "+", RNA_BK1 + 1, RNA_BK1 - 1))
oncofuse %<>% mutate(RNA_BK2 = ifelse(STRAND3p == "+", RNA_BK2 - 1, RNA_BK2 + 1))

results$id <- paste("chr", results$chr1, ":", oncofuse$RNA_BK1, ">",
                    "chr", results$chr2, ":", oncofuse$RNA_BK2, sep = '')

ifn <- str_c(opt$outPrefix, ".oncofuse.input.txt")
write.table(oncofuse[which(!is.na(oncofuse[,"STRAND5p"])),1:5], file=ifn, sep="\t", row.names=F, col.names=F, quote=F, na="")

ofn <- str_c(opt$outPrefix, ".oncofuse.output.txt")
cmd <- paste(opt$java, '-Xmx1G -jar', opt$oncofuseJar, ifn, "coord - ", ofn)
system(cmd, wait = T)

oncofuse_output <- read.delim(ofn, as.is=T)
results <- left_join(results, oncofuse_output, by = c("id" = "GENOMIC"))

fn <- str_c(opt$outPrefix, ".oncofuse.txt")
write.table(results, file=fn, sep="\t", row.names=F, quote=F, na="")

