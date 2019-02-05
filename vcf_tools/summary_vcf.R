#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list = list(make_option("--sample_name", default = NA, type = 'character', help = "sample name"))
				  
parser = OptionParser(usage = "%prog", option_list = args_list)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

vcf = read.csv(file=paste0("cravat/", opt$sample_name, ".vcf"), header=FALSE, sep="\t", comment.char="#", stringsAsFactors=FALSE)
colnames(vcf) = c("Chromosome", "Position", "ID", "Reference", "Alternate", "Base_Quality", "Filter", "Info")
maf = read.csv(file=paste0("cravat/", opt$sample_name, ".maf"), header=TRUE, sep="\t", comment.char="#", stringsAsFactors=FALSE)
index = maf[,"Variant_Classification"] %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site")
vcf = vcf[index,,drop=FALSE]
maf = maf[index,,drop=FALSE]
clinvar = read.csv(file=paste0("cravat/", opt$sample_name, ".clinvar.var"), header=FALSE, sep="\t", comment.char="#", stringsAsFactors=FALSE)
colnames(clinvar) = c("UID", "Clinical_Significance", "Disease_Ref_Nums", "Disease_Names")
cosmic = read.csv(file=paste0("cravat/", opt$sample_name, ".cosmic.var"), header=FALSE, sep="\t", comment.char="#", stringsAsFactors=FALSE)
colnames(cosmic) = c("UID", "ID", "Variant_Count_Tissue", "Variant_Count", "Transcript", "Protein_Change")
dbsnp = read.csv(file=paste0("cravat/", opt$sample_name, ".dbsnp.var"), header=FALSE, sep="\t", comment.char="#", stringsAsFactors=FALSE)
colnames(dbsnp) = c("UID", "dbSNP_rsID")
gnomad = read.csv(file=paste0("cravat/", opt$sample_name, ".gnomad.var"), header=FALSE, sep="\t", comment.char="#", stringsAsFactors=FALSE)
colnames(gnomad) = c("UID", "AF_Total", "AF_African", "AF_American", "AF_Ashkenazi_Jewish", "AF_East_Asian", "AF_Finnish", "AF_Non_Fin_Eur", "AF_Other", "AF_South_Asian")
hgvs = read.csv(file=paste0("cravat/", opt$sample_name, ".hgvs.var"), header=FALSE, sep="\t", comment.char="#", stringsAsFactors=FALSE)
colnames(hgvs) = c("UID", "Genomic_Location", "HGVSc_Primary", "HGVSp_Primary", "HGVSc_All", "HGVSp_All")
res = data.frame(vcf[,c("Chromosome", "Position", "Reference", "Alternate"), drop=FALSE],
				 maf[,c("Hugo_Symbol", "Variant_Classification", "HGVSc", "HGVSp", "HGVSp_Short"), drop=FALSE],
				 stringsAsFactors=FALSE)

tmp = matrix("", nrow=nrow(res), ncol=ncol(clinvar))
for (i in 1:ncol(clinvar)) {
	tmp[clinvar$UID,i] = clinvar[,i]
}
colnames(tmp) = colnames(clinvar)
res = cbind(res, tmp[,-1,drop=FALSE], stringsAsFactors=FALSE)

tmp = matrix("", nrow=nrow(res), ncol=ncol(cosmic))
for (i in 1:ncol(cosmic)) {
	tmp[cosmic$UID,i] = cosmic[,i]
}
colnames(tmp) = colnames(cosmic)
res = cbind(res, tmp[,-1,drop=FALSE], stringsAsFactors=FALSE)

tmp = matrix("", nrow=nrow(res), ncol=ncol(dbsnp))
for (i in 1:ncol(dbsnp)) {
	tmp[dbsnp$UID,i] = dbsnp[,i]
}
colnames(tmp) = colnames(dbsnp)
res = cbind(res, tmp[,-1,drop=FALSE], stringsAsFactors=FALSE)

tmp = matrix("", nrow=nrow(res), ncol=ncol(gnomad))
for (i in 1:ncol(gnomad)) {
	tmp[gnomad$UID,i] = gnomad[,i]
}
colnames(tmp) = colnames(gnomad)
res = cbind(res, tmp[,-1,drop=FALSE], stringsAsFactors=FALSE)

tmp = matrix("", nrow=nrow(res), ncol=ncol(hgvs))
for (i in 1:ncol(hgvs)) {
	tmp[hgvs$UID,i] = hgvs[,i]
}
colnames(tmp) = colnames(hgvs)
res = cbind(res, tmp[,-1,drop=FALSE], stringsAsFactors=FALSE)
res[is.na(res)] = NA
res[res==""] = NA
sample_name = rep(opt$sample_name, nrow(res))
res = cbind("Sample_ID"=sample_name, res, stringsAsFactors=FALSE)
write.table(res, file=paste0("cravat/", opt$sample_name, ".txt"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, append=FALSE)
