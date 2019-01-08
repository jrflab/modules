#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("VariantAnnotation"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene"))
suppressPackageStartupMessages(library("BSgenome.Hsapiens.UCSC.hg19"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))
suppressPackageStartupMessages(library("RMySQL"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("rlist"))
suppressPackageStartupMessages(library("pipeR"))
suppressPackageStartupMessages(library("foreach"))
suppressPackageStartupMessages(library("doMC"))


if (!interactive()) {
    options(warn=-1, error=quote({ traceback(2); q('no', status=1) }))
}
options(useFancyQuotes=F)

optList <- list(
        make_option("--genome", default='b37', help="genome build [default %default]"),
        make_option("--aaTable", default='~/share/reference/aa_table.tsv', help="amino acid table (3-letter to 1 letterIUPAC) [default %default]"),
        make_option("--ensemblTxdb", default='~/share/reference/hsapiens_ensembl_biomart.sqlite', help="Ensembl TxDb SQLite [default %default]"),
        make_option("--mysqlHost", default='10.0.200.71', help="MySQL server hostname [default %default]"),
        make_option("--mysqlPort", default=38493, help="MySQL server port [default %default]"),
        make_option("--mysqlUser", default='embl', help="MySQL server username [default %default]"),
        make_option("--mysqlPassword", default=NULL, help="MySQL server password [default %default]"),
        make_option("--mysqlDb", default='homo_sapiens_core_75_37', help="MySQL server database [default %default]"),
        make_option("--provean", default='provean.sh', help="provean script [default %default]"),
        make_option("--numThreads", default=4, help="Number of provean threads [default %default]"),
        make_option("--memPerThread", default='1G', help="Amount of memory per thread [default %default]"),
        make_option("--qsub", default='modules/scripts/qsub.pl', help="qsub perl script [default %default]"),
        make_option("--qsubPriority", default=-800, help="qsub priority [default %default]"),
        make_option("--queue", default='jrf.q', help="qsub queue [default %default]"),
        make_option("--outFile", default=NULL, help="vcf output file [default %default]")
        )
parser <- OptionParser(usage="%prog vcf.file", option_list=optList);
arguments <- parse_args(parser, positional_arguments=T);
opt <- arguments$options;

registerDoMC(30)

if (is.null(opt$outFile)) {
    cat("Need output file\n");
    print_help(parser);
    stop();
} else if (length(arguments$args) != 1) {
    cat("Need vcf file\n");
    print_help(parser);
    stop();
}

proveanOpts <- paste('--num_threads', opt$numThreads)

aa <- read.table(opt$aaTable, sep='\t', stringsAsFactors=F)

fn <- arguments$args[1];
outfn <- opt$outFile
null <- suppressWarnings(file.remove(outfn))
out <- file(outfn, open='a')

gen <- BSgenome.Hsapiens.UCSC.hg19
seqlevelsStyle(gen) = "NCBI"

cat('Loading transcriptdb ... ')
if (is.null(opt$ensemblTxdb)) {
    txdb <- makeTranscriptDbFromBiomart(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl')
} else {
    txdb <- loadDb(opt$ensemblTxdb)
}
cat('done\n')
cdsByTx <- cdsBy(txdb, 'tx', use.names=T)

cat('Connecting to ensembl ... ')
connect <- function() dbConnect(MySQL(), host = opt$mysqlHost, port = opt$mysqlPort, user = opt$mysqlUser, password = opt$mysqlPassword, dbname = opt$mysqlDb)
mydb <- connect()
on.exit(dbDisconnect(mydb))
cat('done\n')

cat('Reading vcf header ... ')
# create new header
vcfHeader <- scanVcfHeader(fn)
hinfo <- apply(as.data.frame(info(vcfHeader)), 2, as.character)
rownames(hinfo) <- rownames(info(vcfHeader))
hinfo <- rbind(hinfo, provean_embl_id=c("A", "String", "provean query embl ID"))
hinfo <- rbind(hinfo, provean_hgvsp=c("A", "String", "provean query HGVSp"))
hinfo <- rbind(hinfo, provean_score=c("A", "Float", "provean score"))
hinfo <- DataFrame(hinfo, row.names=rownames(hinfo))
hlist <- header(vcfHeader)
hlist$INFO <- hinfo
newVcfHeader <- new("VCFHeader", samples=vcfHeader@samples, header=hlist)
cat('done\n')

cat('Indexing vcf ... ')
temp <- tempfile()
zipped <- bgzip(fn, temp)
idx <- indexTabix(temp, "vcf")
cat('done\n')

tab <- TabixFile(zipped, idx, yieldSize=2000)
open(tab)

cat('Processing vcf by chunk\n')
i <- 1
while(nrow(vcf <- readVcf(tab, genome=opt$genome))) {
    oldwd <- getwd()
    # replace header
    metadata(vcf)$header <- newVcfHeader

    cat(paste('Chunk', i, "\n"))
    i <- i + 1
    passIds <- which(rowRanges(vcf)$FILTER == "PASS" & seqnames(rowRanges(vcf)) %in% c(1:22, "X", "Y"))
    if (length(passIds) == 0) {
        cat("No unfiltered variants\n")
    } else {
        vcf <- vcf[passIds, ]
        cat("Predicting coding from reference...\n")
        if (any(!seqlevels(vcf) %in% seqlevels(txdb))) {
            seqlevels(vcf, force=T) <- seqlevels(vcf)[-which(!seqlevels(vcf) %in% seqlevels(txdb))]
        }
        ann <- info(vcf)$ANN %>>% list.map(f(x, i) ~ paste(i, x, sep = "|"))
        ann <- ann %>>% list.map(f(x) ~ strsplit(x, '\\|')) %>>% list.flatten()
        rowIds <- ann %>>% list.map(.[1]) %>>% unlist %>>% as.integer
        ids <- ann %>>% list.map(.[8]) %>>% unlist
        hgvsp <- ann %>>% list.map(.[12]) %>>% unlist
        hgvsp <- sub('^p.', '', hgvsp)
        for (i in 1:nrow(aa)) {
            hgvsp <- gsub(aa[i, 1], aa[i, 2], hgvsp)
        }
        Df <- data.frame(rowId=rowIds, refseqId=ids, HGVSp=hgvsp, stringsAsFactors=F)
        Df <- Df %>% filter(HGVSp != "" & !grepl('fs', HGVSp))

        if (sum(Df$HGVSp == "") == 0) {
            cat("No variants to query\n")
        } else {
            query <- paste("SELECT X.display_label AS refseqId, T.stable_id as emblId FROM object_xref as O JOIN xref as X ON O.xref_id=X.xref_id JOIN transcript as T ON O.ensembl_id=T.transcript_id WHERE X.display_label IN (", paste(sQuote(Df$refseqId), collapse=','), ");")
            cat(paste(query, "\n", sep=""));
            cat("Looking up ensembl transcript IDs ...\n")
            repeat {
                rs <- try(dbSendQuery(mydb, query), silent=F)
                if (is(rs, "try-error")) {
                    dbDisconnect(mydb)
                    cat("Lost connection to mysql db ... ")
                    Sys.sleep(10)
                    mydb <- connect()
                    cat("reconnected\n")
                } else {
                    break
                }
            }
            results <- fetch(rs, -1)
            dbDisconnect(mydb)
            cat(paste("Found", nrow(results), "records\n"))
            if (nrow(results) > 0 && ncol(results) > 0) {
                Df <- right_join(Df, results)
                cds <- cdsByTx[Df$emblId]
                cdsSeqs <- extractTranscriptSeqs(gen, cds)
                protSeqs <- translate(cdsSeqs, if.fuzzy.codon='solve')
                names(protSeqs) <- Df$refseqId
                sDf <- split(Df, Df$refseqId)
                proveanOutputs <- foreach(i=1:length(sDf)) %dopar% {
                    tmpFasta <- tempfile()
                    tmpVar <- tempfile()
                    pout <- tempfile()
                    writeXStringSet(protSeqs[names(sDf)[i]], tmpFasta)
                    write(sDf[[i]][, "HGVSp"], file=tmpVar)
                    cmd <- paste('echo "', opt$provean, ' ', proveanOpts, ' -q ', tmpFasta, ' -v ', tmpVar,
                                 ' | sed \\"1,/PROVEAN scores/d\\" > ', pout, '" | ', opt$qsub,
                                 ' -- -V -b n -o /dev/null -p ', opt$qsubPriority, ' -j y -q ', opt$queue,
                                 ' -pe smp ', opt$numThreads, ' -l h_vmem=', opt$memPerThread, sep='')
                    cat("Running:", cmd, "\n")
                    system(cmd)
                    if (file.exists(pout) && length(readLines(pout)) > 0) {
                        proveanOutput <- read.table(pout, sep = '\t', stringsAsFactors=F)
                        colnames(proveanOutput) <- c('HGVSp', 'proveanScore')
                        proveanOutput$refseqId <- names(protSeqs)[i]
                        proveanOutput
                    } else {
                        NA
                    }
                }
                proveanOutputs <- proveanOutputs[!is.na(proveanOutputs)]
                if (length(proveanOutputs) > 0) {
                    X <- do.call('rbind', proveanOutputs)
                    Df <- left_join(Df, X)
                } else {
                    Df$proveanScore <- NA
                }
            } else {
                Df$proveanScore <- NA
            }
            if (!is.null(results) && nrow(results) > 0) {
                cat("Merging provean results ... ")
                infodf <- info(vcf)
                infodf[Df$rowId,"provean_score"] <- Df$proveanScore
                infodf[Df$rowId,"provean_hgvsp"] <- Df$HGVSp
                infodf[Df$rowId,"provean_embl_id"] <- Df$emblId
                info(vcf) <- infodf
                cat("done\n")
            } else {
                cat("No results from provean\n")
            }
        }
    }

    #fix sample genotype order
    if ("GT" %in% names(geno(vcf))) {
        x <- which(names(geno(vcf)) == "GT")
        ord <- c(x, (1:length(geno(vcf)))[-x])
        geno(vcf) <- geno(vcf)[ord]
    }

    cat("Appending vcf chunk to", opt$outFile, "... ")
    setwd(oldwd)
    writeVcf(vcf, out)
    cat("done\n")
}

if (i == 1) {
    cat("No entries, creating empty vcf file\n")
    writeVcf(vcf, out)
}

close(tab)
close(out)


