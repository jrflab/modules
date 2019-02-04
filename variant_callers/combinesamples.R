#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--sample_set", default = NA, type = 'character', help = "sample names set"))
				  
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

sample_names = unlist(strsplit(opt$sample_set, split="_", fixed=TRUE))

all_vars = read.csv(file="summary/tsv/mutation_summary.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
tmp_vars = all_vars[all_vars$TUMOR_SAMPLE %in% sample_names,,drop=FALSE]
keys = paste0(tmp_vars$CHROM, ":", tmp_vars$POS, ":", tmp_vars$REF, ":", tmp_vars$ALT)
ukeys = unique(keys)
vars = NULL
for (i in 1:length(ukeys)) {
	index = which(keys==ukeys[i])
 	Chromosome = tmp_vars[index[1],"CHROM"]
	Position = tmp_vars[index[1],"POS"]
	Ref = tmp_vars[index[1],"REF"]
	Alt = tmp_vars[index[1],"ALT"]
	Variant_Caller = tmp_vars[index[1],"variantCaller"]
	Gene_Symbol = tmp_vars[index[1],"SYMBOL"]
	Variant_Classification = tmp_vars[index[1],"Variant_Classification"]
	HGVSp_Short = tmp_vars[index[1],"HGVSp_Short"]
	Fuentes = tmp_vars[index[1],"fuentes"]
	dgd = tmp_vars[index[1],"dgd"]
	OncoKB_Level = tmp_vars[index[1],"oncoKB_level"]
	OncoKB_Cancer_Type = tmp_vars[index[1],"oncoKB_cancer_type"]
	Cancer_Gene_Census = tmp_vars[index[1],"cancer_gene_census"]
	Kandoth = tmp_vars[index[1],"kandoth"]
	Lawrence = tmp_vars[index[1],"lawrence"]
	Hap_Insuf = tmp_vars[index[1],"hap_insuf"]
	ExAC_AF = tmp_vars[index[1],"ExAC_AF"]
	MutationTaster = tmp_vars[index[1],"MutationTaster_pred"]
	PROVEAN = tmp_vars[index[1],"PROVEAN_pred"]
	FATHMM = tmp_vars[index[1],"FATHMM_pred"]
	BRCA_Chasm = tmp_vars[index[1],"BRCA_chasm_pred"]
	Parssnp = tmp_vars[index[1],"parssnp_pred"]
	Pathogenicity = tmp_vars[index[1],"pathogenicity"]
	HOTSPOT = tmp_vars[index[1],"HOTSPOT"]
	HOTSPOT_INTERNAL = tmp_vars[index[1],"HOTSPOT_INTERNAL"]
	CMO_HOTSPOT = tmp_vars[index[1],"cmo_hotspot"]
	vars = rbind(vars, c("Chromosome"=Chromosome,
					   	 "Position"=Position,
					     "Ref"=Ref,
					     "Alt"=Alt,
					     "Variant_Caller"=Variant_Caller,
					     "Gene_Symbol"=Gene_Symbol,
					     "Variant_Classification"=Variant_Classification,
					     "HGVSp"=HGVSp_Short,
					     "Fuentes"=Fuentes,
					     "dgd"=dgd,
					     "OncoKB_Level"=OncoKB_Level,
					     "OncoKB_Cancer_Type"=OncoKB_Cancer_Type,
					     "Cancer_Gene_Census"=Cancer_Gene_Census,
					     "Kandoth"=Kandoth,
					     "Lawrence"=Lawrence,
					     "Hap_Insuf"=Hap_Insuf,
					     "ExAC"=ExAC_AF,
					     "MutationTaster"=MutationTaster,
					     "PROVEAN"=PROVEAN,
					     "FATHMM"=FATHMM,
					     "BRCA_Chasm"=BRCA_Chasm,
					     "Parssnp"=Parssnp,
					     "Pathogenicity"=Pathogenicity,
					     "HOTSPOT"=HOTSPOT,
					     "HOTSPOT_INTERNAL"=HOTSPOT_INTERNAL,
					     "HOTSPOT_CMO"=CMO_HOTSPOT))
}

normal_name = tmp_vars[1,"NORMAL_SAMPLE"]

VAF = DEPTH = LOH = CALLS = matrix(NA, nrow=length(ukeys), ncol=length(sample_names), dimnames=list(ukeys, sample_names))
for (j in 1:nrow(tmp_vars)) {
	sample_name = tmp_vars[j,"TUMOR_SAMPLE"]
	ukey = paste0(tmp_vars$CHROM[j], ":", tmp_vars$POS[j], ":", tmp_vars$REF[j], ":", tmp_vars$ALT[j])
	VAF[ukey,sample_name] = tmp_vars[j,"TUMOR_MAF"]
	VAF[ukey,normal_name] = tmp_vars[j,"NORMAL_MAF"]
	DEPTH[ukey,sample_name] = tmp_vars[j,"TUMOR_DP"]
	DEPTH[ukey,normal_name] = tmp_vars[j,"NORMAL_DP"]
	LOH[ukey,sample_name] = tmp_vars[j,"facetsLOHCall"]
	CALLS[ukey,sample_name] = 1
}
colnames(VAF) = paste0("MAF_", colnames(VAF))
colnames(DEPTH) = paste0("DP_", colnames(DEPTH))
colnames(LOH) = paste0("LOH_", colnames(LOH))
colnames(CALLS) = paste0("CALL_", colnames(CALLS))
CALLS[is.na(CALLS)] = 0
vars = cbind(vars, VAF, DEPTH, LOH, CALLS)
mutect = grepl("mutect", vars[,"Variant_Caller"])
main_indels = grepl("varscan", vars[,"Variant_Caller"]) & grepl("strelka", vars[,"Variant_Caller"])
other_indels = ((grepl("platypus", vars[,"Variant_Caller"]) & grepl("scalpel", vars[,"Variant_Caller"])) |
			   (grepl("platypus", vars[,"Variant_Caller"]) & grepl("lancet", vars[,"Variant_Caller"]))) &
			   (nchar(vars[,"Ref"])>3 | nchar(vars[,"Alt"])>3) &
			   !grepl("In_Frame", vars[,"Variant_Classification"])
index = mutect | main_indels | other_indels
vars = vars[index,,drop=FALSE]
index = vars[,"Variant_Classification"] %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Site")
vars = vars[index,,drop=FALSE]

write.table(vars, file=paste0("sufam/", opt$sample_set, ".txt"), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
