#!/usr/bin/env Rscript
# run oncofuse on defuse output

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
                make_option("--java", default = 'java', help = "java binary"),
                make_option("--oncofuseJar", default = '~/share/usr/oncofuse-v1.0.6/Oncofuse.jar', help = "oncofuse jar"),
                make_option("--oncofuseTissueType", default = 'EPI', help = "oncofuse tissue type"),
                make_option("--outPrefix", default = NULL, help = "Output file"));

parser <- OptionParser(usage = "%prog [options] [defuse file]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;


if (is.null(opt$outPrefix)) {
    cat("Need output prefix\n");
    print_help(parser);
    stop();
} else if (length(arguments$args) < 1) {
    cat("Need input file\n");
    print_help(parser);
    stop();
}


fn <- arguments$args[1]
defuse <- read.table(fn, sep = '\t', header = T, stringsAsFactors = F)
upstreamChr <- defuse %$% ifelse(upstream_gene == gene_name1, gene_chromosome1, gene_chromosome2)
downstreamChr <- defuse %$% ifelse(downstream_gene == gene_name1, gene_chromosome1, gene_chromosome2)
upstreamStrand <- defuse %$% ifelse(upstream_gene == gene_name1, gene_strand1, gene_strand2)
downstreamStrand <- defuse %$% ifelse(downstream_gene == gene_name1, gene_strand1, gene_strand2)
upstreamGenomicBreakPos <- defuse %$% ifelse(upstream_gene == gene_name1, genomic_break_pos1, genomic_break_pos2)
downstreamGenomicBreakPos <- defuse %$% ifelse(downstream_gene == gene_name1, genomic_break_pos1, genomic_break_pos2)
upstreamPos <- defuse %$% ifelse(upstreamStrand == "+", upstreamGenomicBreakPos + 1, upstreamGenomicBreakPos - 1)
downstreamPos <- defuse %$% ifelse(downstreamStrand == "+", downstreamGenomicBreakPos - 1, downstreamGenomicBreakPos + 1)

defuse$id <- paste('chr', upstreamChr, ':', upstreamPos, '>', 'chr', downstreamChr, ':', downstreamPos, sep = '')
tissueType <- opt$oncofuseTissueType
oncofuseInput <- data.frame(paste('chr', upstreamChr, sep = ''), upstreamPos,
                            paste('chr', downstreamChr, sep = ''), downstreamPos,
                            tissueType)
ifn <- str_c(opt$outPrefix, ".oncofuse.input.txt")
write.table(oncofuseInput, file=ifn, sep="\t", row.names=F, col.names=F, quote=F, na="")

ofn <- str_c(opt$outPrefix, ".oncofuse.output.txt")
cmd <- paste(opt$java, '-Xmx1G -jar', opt$oncofuseJar, ifn, "coord - ", ofn)
system(cmd, wait = T)

oncofuseOutput <- read.delim(ofn, as.is=T)
results <- left_join(defuse, oncofuseOutput, by = c("id" = "GENOMIC"))
results <- select(results, -SAMPLE_ID)

fn <- str_c(opt$outPrefix, ".oncofuse.txt")
write.table(results, file=fn, sep="\t", row.names=F, quote=F, na="")


