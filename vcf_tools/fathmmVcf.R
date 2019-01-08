#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("VariantAnnotation"))
suppressPackageStartupMessages(library("data.table"))
suppressPackageStartupMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene"))
suppressPackageStartupMessages(library("BSgenome.Hsapiens.UCSC.hg19"))
suppressPackageStartupMessages(library("org.Hs.eg.db"))
suppressPackageStartupMessages(library("RMySQL"))
suppressPackageStartupMessages(library("foreach"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(2); q('no', status = 1) }))
}
options(useFancyQuotes = F)

optList <- list(
        make_option("--genome", default = 'b37', help = "genome build [default %default]"),
        make_option("--fathmmDir", default = '~/share/usr/fathmm', help = "fathmm dir"),
        make_option("--fathmmConfig", help = "fathmm config file", default = NULL),
        make_option("--fathmmAlg", default = 'Cancer', help = "fathmm algorithm [default %default]"),
        make_option("--fathmmOnt", default = 'DO', help = "fathmm ontology [default %default]"),
        make_option("--ensemblTxdb", default = '~/share/reference/hsapiens_ensembl_biomart.sqlite', help = "Ensembl TxDb SQLite"),
        make_option("--mysqlHost", default = '10.0.200.48', help = "MySQL server hostname"),
        make_option("--mysqlPort", default = 38493, help = "MySQL server port"),
        make_option("--mysqlUser", default = 'embl', help = "MySQL server username"),
        make_option("--mysqlPassword", default = NULL, help = "MySQL server password"),
        make_option("--mysqlDb", default = 'homo_sapiens_core_75_37', help = "MySQL server database"),
        make_option("--ref", default = '~/share/reference/GATK_bundle/2.3/human_g1k_v37.fasta', help = "Reference fasta file"),
        make_option("--python", default = 'python', help = "python executable [default %default]"),
        make_option("--outFile", default = NULL, help = "vcf output file [default %default]")
        )
parser <- OptionParser(usage = "%prog vcf.file", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$fathmmDir)) {
    cat("Need fathmm dir\n");
    print_help(parser);
    stop();
} else if (is.null(opt$outFile)) {
    cat("Need output file\n");
    print_help(parser);
    stop();
} else if (is.null(opt$ref)) {
    cat("Need reference fasta file\n");
    print_help(parser);
    stop();
} else if (length(arguments$args) != 1) {
    cat("Need vcf file\n");
    print_help(parser);
    stop();
}

fn <- arguments$args[1];
outfn <- opt$outFile
null <- suppressWarnings(file.remove(outfn))
out <- file(outfn, open = 'a')

cat('Loading transcriptdb ... ')
if (is.null(opt$ensemblTxdb)) {
    txdb <- makeTranscriptDbFromBiomart(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl')
} else {
    txdb <- loadDb(opt$ensemblTxdb)
}
if (opt$genome == 'hg19') {
    seqlevels(txdb) <- sub('^', 'chr', seqlevels(txdb))
}
cat('done\n')

ref <- FaFile(opt$ref)

cat('Reading vcf header ... ')
# create new header
vcfHeader <- scanVcfHeader(fn)
hinfoprime <- as.data.frame(info(vcfHeader))
hinfoprime <- rbind(hinfoprime, fathmm_query = c("Number" = "A", "Type" = "String", "Description" = "fathmm query"))
hinfoprime <- rbind(hinfoprime, fathmm_pred = c("Number" = "A", "Type" = "String", "Description" = "fathmm prediction"))
hinfoprime <- rbind(hinfoprime, fathmm_score = c("Number" = "A", "Type" = "Float", "Description" = "fathmm score"))
hinfoprime <- DataFrame(hinfoprime, row.names = rownames(hinfoprime))
hlist <- header(vcfHeader)
hlist$INFO <- hinfoprime
colnames(hlist$INFO) <- c('Number', 'Type', 'Description')
newVcfHeader <- new("VCFHeader", samples = vcfHeader@samples, header = hlist)
cat('done\n')


cat('Indexing vcf ... ')
temp <- tempfile()
zipped <- bgzip(fn, temp)
idx <- indexTabix(temp, "vcf")
cat('done\n')

tab <- TabixFile(zipped, idx, yieldSize = 2000)
open(tab)

cat('Connecting to ensembl ... ')
connect <- function() dbConnect(MySQL(), host = opt$mysqlHost, port = opt$mysqlPort, user = opt$mysqlUser, password = opt$mysqlPassword, dbname = opt$mysqlDb)
mydb <- connect()
on.exit(dbDisconnect(mydb))
cat('done\n')

cat('Processing vcf by chunk\n')
i <- 1
while(nrow(vcf <- readVcf(tab, genome = opt$genome))) {
    oldwd <- getwd()
    # replace header
    metadata(vcf)$header <- newVcfHeader

    cat(paste('Chunk', i, "\n"))
    i <- i + 1

    passIds <- which(rowRanges(vcf)$FILTER == "PASS" & (seqnames(rowRanges(vcf)) %in% c(1:22, "X", "Y") | seqnames(rowRanges(vcf)) %in% paste('chr', c(1:22, "X", "Y"), sep = '')))
    if (length(passIds) == 0) {
        cat("No unfiltered variants\n")
    } else {
        cat(length(passIds), "variants pass\n")

        cat("Predicting coding from reference...\n")
        vcfPass <- vcf[passIds, ]
        if (any(!seqlevels(vcfPass) %in% seqlevels(txdb))) {
            seqlevels(vcfPass, force = T) <- seqlevels(vcfPass)[-which(!seqlevels(vcfPass) %in% seqlevels(txdb))]
        }
        predCod <- predictCoding(vcfPass, txdb, ref)
        cat(" done\n")

        if (sum(predCod$CONSEQUENCE == "nonsynonymous") == 0) {
            cat("No non-syn variants\n")
        } else {
            cat(sum(predCod$CONSEQUENCE == "nonsynonymous"), "non-syn variants\n")

            predCod <- subset(predCod, CONSEQUENCE == "nonsynonymous")

            # retrieve transcript ids
            x <- transcripts(txdb, filter = list(tx_id = predCod$TXID), columns = c('tx_id', 'tx_name'))
            enstIds <- x$tx_name
            names(enstIds) <- x$tx_id
            aa = cbind(queryId = passIds[predCod$QUERYID], aa = paste(as.character(predCod$REFAA), lapply(predCod$PROTEINLOC, function(x) x[1]), as.character(predCod$VARAA), sep = ''))
            rownames(aa) <- predCod$TXID

            splitEnstIds <- split(enstIds, ceiling(seq_along(enstIds) / 50))
            idList <- list()
            for (ss in names(splitEnstIds)) {
                #cat(ss, "\n")
                query <- paste("SELECT P.stable_id AS peptide_id, T.stable_id AS transcript_id
                            from transcript as T JOIN translation as P ON T.transcript_id = P.transcript_id
                            where T.stable_id in (", paste(sQuote(unlist(splitEnstIds[[ss]])), collapse = ','), ");")
                cat(paste(query, "\n", sep = ""));
                cat("Looking up ensembl peptide IDs ...\n")
                repeat {
                    rs <- try(dbSendQuery(mydb, query), silent = T)
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
                idList[[ss]] <- fetch(rs, -1)
                dbClearResult(rs)
            }
            ids <- do.call('rbind', idList)
            cat(paste("Found", nrow(ids), "records\n"))
            #ids <- getBM(filters = 'ensembl_transcript_id', attributes = c('ensembl_transcript_id', 'ensembl_peptide_id'), values = enstIds, mart = ensembl)
            if (nrow(ids) > 0 && ncol(ids) > 0) {
                rownames(ids) <- names(enstIds)[match(ids$transcript_id, enstIds)]
                m <- match(rownames(aa), rownames(ids))
                ids <- cbind(aa[!is.na(m), , drop = F], ids[m[!is.na(m)], , drop = F])
                cat("done\n")

                fathmmInput <- subset(ids, peptide_id != "", select = c('peptide_id', 'aa'))

                cat("Calling fathmm: ")
                tmp1 <- tempfile()
                tmp2 <- tempfile()
                setwd(paste(opt$fathmmDir, '/cgi-bin', sep = ''))
                cmd <- paste(opt$python, 'fathmm.py')
                if (!is.null(opt$fathmmConfig)) {
                    cmd <- paste(cmd, '-c', opt$fathmmConfig)
                }
                cmd <- paste(cmd, '-w', opt$fathmmAlg, '-p', opt$fathmmOnt, tmp1, tmp2)
                write.table(subset(ids, peptide_id != "", select = c('peptide_id', 'aa')), file = tmp1, quote = F, sep = ' ', row.names = F, col.names = F)
                #cmd <- paste('python fathmm.py -w Cancer', tmp1, tmp2)
                system(cmd)
                cat("\ndone\n")

                cat("Reading results ... ")
                results <- read.table(tmp2, sep = '\t', header = T, comment.char = '', row.names = 1, quote = '')
                cat("done\n")
                results <- merge(ids, results, by.x = c('aa', 'peptide_id'), by.y = c('Substitution', 'Protein.ID'))

                split.results <- split(results, factor(results$queryId))
                cat("Selecting minimum scores ... ")
                results <- rbindlist(lapply(split.results, function(x) x[which.min(x$Score), ]))
                cat("done\n")
            } else {
                results <- NULL
            }
            if (!is.null(results) && nrow(results) > 0) {
                cat("Merging fathmm results ... ")
                infodprime <- info(vcf)
                infodprime[as.integer(as.character(results$queryId)),"fathmm_query"] <- with(results, paste(peptide_id, aa, sep = "_"))
                infodprime[as.integer(as.character(results$queryId)),"fathmm_pred"] <- as.character(results$Prediction)
                infodprime[as.integer(as.character(results$queryId)),"fathmm_score"] <- results$Score
                info(vcf) <- infodprime
                cat("done\n")
            } else {
                cat("No results from fathmm\n")
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


