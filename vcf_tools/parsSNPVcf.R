#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("VariantAnnotation"))
suppressPackageStartupMessages(library("nnet"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(2); q('no', status = 1) }))
}
options(useFancyQuotes = F)

optList <- list(
        make_option("--genome", default = 'b37', help = "genome build [default %default]"),
        make_option("--parsnpRdata", default = NULL, help = "genome build [default %default]"),
        make_option("--outFile", default = NULL, help = "vcf output file [default %default]")
        )
parser <- OptionParser(usage = "%prog vcf.file", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$parsnpRdata)) {
    cat("Need parsnp Rdata\n");
    print_help(parser);
    stop();
} else if (length(arguments$args) != 1) {
    cat("Need vcf file\n");
    print_help(parser);
    stop();
}

fn <- arguments$args[1];

cat('Reading vcf header ... ')
# create new header
vcfHeader <- scanVcfHeader(fn)
hinfo <- apply(as.data.frame(info(vcfHeader)), 2, as.character)
rownames(hinfo) <- rownames(info(vcfHeader))
hinfo <- rbind(hinfo, parssnp_score = c("A", "Float", "parsSNP score"))
hinfo <- rbind(hinfo, parssnp_pred = c("A", "String", "parsSNP prediction"))
hinfo <- DataFrame(hinfo, row.names = rownames(hinfo))
hlist <- header(vcfHeader)
hlist$INFO <- hinfo
newVcfHeader <- new("VCFHeader", samples = vcfHeader@samples, header = hlist)
cat('done\n')


vcf <- readVcf(fn, genome = opt$genome)
metadata(vcf)$header <- newVcfHeader

load(opt$parsnpRdata)


print("Input and supporting files read.")
flush.console()

#First, we replace any missing values from the exonic function column with values from the more general function column.
func <- sapply(info(vcf)$ExonicFunc_refGene, function(x) x[1])
if (any(func == '.')) {
    func[func == '.'] <- sapply(info(vcf)[func == '.', 'Func_refGene'], function(x) x[1])
    func[func == '.'] <- NA
}

#In some cases, the gene and exonic function will contain two or more terms separated by a semicolon. These are removed.

#Now we edit the amino acid changes, which currently have many alternative transcripts annotated.
#Only the first transcript is retained; then the actual protein change is parsed out. 
aachange <- sapply(info(vcf)$AAChange_refGene, function(x) x[1])
aachange <- gsub(".+p\\.(.+)", "\\1", aachange)

#Now we add some new descriptors. First, any genes that are in the dataset, but which do not have gene-level 
#descriptors, are assigned mean values in the gene.anno table. 
gene.anno[setdiff(unlist(info(vcf)$Gene_refGene), rownames(gene.anno)),]<-c(NA, NA, NA, NA)
gene.anno<-apply(gene.anno,2,function(i) {i[is.na(i)]<-mean(i,na.rm=T)
	return(i)})
gene.anno<-as.data.frame(gene.anno)

#Gene level annotations are appended to the dataset.
g <- sapply(info(vcf)$Gene_refGene, function(x) x[1])
d <- data.frame(func, AAChange = aachange, gene.anno[g, ], stringsAsFactors = F)

#ParsSNP has access to indicator variables that show if the mutations does NOT affect protein (e.g. silent), or is a truncation event.
d$VarClassS<-c(0, 1)[1 + (d$func %in% c("ncRNA_exonic", "ncRNA_splicing", "synonymous SNV", "UTR3", "UTR5"))]
d$VarClassT<-c(0, 1)[1 + (d$func %in% c("frameshift deletion", "frameshift insertion", "splicing", "stopgain", "stoploss"))]

#Now the B62 score is added. 
#First, the amino acid positions are further parsed out. 
aa<-gsub("(\\D+)(\\d+)(\\D+)", "\\1 \\3 \\2", d$AAChange, perl=T)
aa[is.na(aa) | aa == '.']<-"UK UK UK"
aa<-as.data.frame(do.call("rbind", strsplit(aa, " ")), stringsAsFactors=F)
colnames(aa)<-c("Old_AA", "New_AA", "Peptide_Position")

#UK stands for "unknown", as in the mutation alteration is not simply one codon into another.
aa[["Old_AA"]][!(aa[["Old_AA"]] %in% c(toupper(letters), "UK"))]<-"UK"
aa[["New_AA"]][!(aa[["New_AA"]] %in% c(toupper(letters), "UK"))]<-"UK"
aa[["Peptide_Position"]]<-as.numeric(aa[["Peptide_Position"]])
aa[["Old_AA"]]<-as.character(aa[["Old_AA"]])
aa[["New_AA"]]<-as.character(aa[["New_AA"]])

#Then the B62 matrix is appended to ensure that stop codons "X" and mutations without defineable codons "UK" can be assigned a score.
b62<-b62[-(21:23), -(21:23)]
colnames(b62)[21]<-"X"
rownames(b62)[21]<-"X"
b62["UK",]<-0
b62[,"UK"]<-0


d$b62<-apply(aa[,c("Old_AA", "New_AA")], 1, function(x) b62[x[1], x[2]])
d$Norm_Position<-as.numeric(as.character(aa[["Peptide_Position"]]))/d$Prot_Length
d$Norm_Position[d$Norm_Position>1]<-1
d<-cbind(d, aa)

scoreCols <- c("SIFT_score",
"Polyphen2_HDIV_score",
"Polyphen2_HVAR_score",
"LRT_score",
"MutationTaster_score",
"MutationAssessor_score",
"FATHMM_score",
"RadialSVM_score",
"LR_score",
"VEST3_score",
"CADD_raw",
"CADD_phred",
"phyloP46way_placental",
"phyloP100way_vertebrate",
"GERP++_RS",
"SiPhy_29way_logOdds")

X <- as.data.frame(info(vcf)[, scoreCols])
scores <- structure(apply(X, 2, function(x) as.numeric(sapply(x, function(y) y[1]))), dim = dim(X))
colnames(scores) <- make.names(scoreCols)

d <- cbind(d, scores)


#Next, we set up the x matrix, which has columns corresponding to scaling, FIS.imp, missing, etc. 
x<-d[rownames(scaling)]


#We ensure that X is comletely numeric.
x<-as.data.frame(structure(apply(x, 2, function(i) as.numeric(as.character(i))), dim = dim(x)))
colnames(x) <- rownames(scaling)


#First, mutations that truncate the protrein (presumed to be loss-of-function changes, "lof"), or do not affect the protein (silent) are marked.
lof<-d$func %in% c("frameshift deletion", "frameshift insertion", "splicing", "stopgain", "stoploss", "nonframeshift deletion", "nonframeshift insertion")
sil<-d$func %in% c("ncRNA_exonic", "ncRNA_splicing", "synonymous SNV", "UTR3", "UTR5") 


#For each of the 16 impact scores, we replace the necessary values.
x$SIFT_score[lof & is.na(x$SIFT_score)]<-FIS.imp$sift[1]
x$SIFT_score[sil & is.na(x$SIFT_score)]<-FIS.imp$sift[2]

x$Polyphen2_HDIV_score[lof & is.na(x$Polyphen2_HDIV_score)]<-FIS.imp$hdiv[1]
x$Polyphen2_HDIV_score[sil & is.na(x$Polyphen2_HDIV_score)]<-FIS.imp$hdiv[2]

x$Polyphen2_HVAR_score[lof & is.na(x$Polyphen2_HVAR_score)]<-FIS.imp$hvar[1]
x$Polyphen2_HVAR_score[sil & is.na(x$Polyphen2_HVAR_score)]<-FIS.imp$hvar[2]

x$LRT_score[lof & is.na(x$LRT_score)]<-FIS.imp$lrt[1]
x$LRT_score[sil & is.na(x$LRT_score)]<-FIS.imp$lrt[2]

x$MutationTaster_score[lof & is.na(x$MutationTaster_score)]<-FIS.imp$mt[1]
x$MutationTaster_score[sil & is.na(x$MutationTaster_score)]<-FIS.imp$mt[2]

x$MutationAssessor_score[lof & is.na(x$MutationAssessor_score)]<-FIS.imp$ma[1]
x$MutationAssessor_score[sil & is.na(x$MutationAssessor_score)]<-FIS.imp$ma[2]

x$FATHMM_score[lof & is.na(x$FATHMM_score)]<-FIS.imp$fathmm[1]
x$FATHMM_score[sil & is.na(x$FATHMM_score)]<-FIS.imp$fathmm[2]

x$RadialSVM_score[lof & is.na(x$RadialSVM_score)]<-FIS.imp$svm[1]
x$RadialSVM_score[sil & is.na(x$RadialSVM_score)]<-FIS.imp$svm[2]

x$LR_score[lof & is.na(x$LR_score)]<-FIS.imp$lr[1]
x$LR_score[sil & is.na(x$LR_score)]<-FIS.imp$lr[2]

x$VEST3_score[lof & is.na(x$VEST3_score)]<-FIS.imp$vest[1]
x$VEST3_score[sil & is.na(x$VEST3_score)]<-FIS.imp$vest[2]

x$CADD_raw[lof & is.na(x$CADD_raw)]<-FIS.imp$caddraw[1]
x$CADD_raw[sil & is.na(x$CADD_raw)]<-FIS.imp$caddraw[2]

x$CADD_phred[lof & is.na(x$CADD_phred)]<-FIS.imp$caddphred[1]
x$CADD_phred[sil & is.na(x$CADD_phred)]<-FIS.imp$caddphred[2]

x$GERP.._RS[lof & is.na(x$GERP.._RS)]<-FIS.imp$gerp[1]
x$GERP.._RS[sil & is.na(x$GERP.._RS)]<-FIS.imp$gerp[2]

x$phyloP46way_placental[lof & is.na(x$phyloP46way_placental)]<-FIS.imp$pp46[1]
x$phyloP46way_placental[sil & is.na(x$phyloP46way_placental)]<-FIS.imp$pp46[2]

x$phyloP100way_vertebrate[lof & is.na(x$phyloP100way_vertebrate)]<-FIS.imp$pp100[1]
x$phyloP100way_vertebrate[sil & is.na(x$phyloP100way_vertebrate)]<-FIS.imp$pp100[2]

x$SiPhy_29way_logOdds[lof & is.na(x$SiPhy_29way_logOdds)]<-FIS.imp$sp29[1]
x$SiPhy_29way_logOdds[sil & is.na(x$SiPhy_29way_logOdds)]<-FIS.imp$sp29[2]

#Then we scale X into [0,1] for each descriptor.
x<-sweep(x, 2, scaling[,2], "-")
x<-sweep(x, 2, scaling[,1]-scaling[,2], "/")

#Finally, missing values are replaced.
x[is.na(x)]<-do.call("rbind", rep(missing, nrow(x)))[is.na(x)]

score <- as.numeric(predict(ParsSNP, x))
info(vcf)$parssnp_score <- score
info(vcf)$parssnp_pred <- ifelse(score <= 0.08, 'Passenger', ifelse(score < 0.16, 'Indeterminate', 'Driver'))

outfn <- opt$outFile
null <- suppressWarnings(file.remove(outfn))
out <- file(outfn, open = 'a')
writeVcf(vcf, out)
close(out)


