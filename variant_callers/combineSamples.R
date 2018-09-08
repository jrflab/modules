#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"));

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--patient", default = NA, type = 'character', help = "type of analysis"))
				  
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

all_vars = read.csv(file="summary/tsv/mutation_summary.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
tmp_vars = all_vars[all_vars$NORMAL_SAMPLE==opt$patient,,drop=FALSE]
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
VAF = DEPTH = LOH = CCF = Clonal_Status = matrix(NA, nrow=length(ukeys), ncol=length(unique(tmp_vars$TUMOR_SAMPLE))+1, dimnames=list(ukeys, c("N", unique(tmp_vars$TUMOR_SAMPLE))))
for (j in 1:nrow(tmp_vars)) {
	sample_name = tmp_vars[j,"TUMOR_SAMPLE"]
	ukey = paste0(tmp_vars$CHROM[j], ":", tmp_vars$POS[j], ":", tmp_vars$REF[j], ":", tmp_vars$ALT[j])
	VAF[ukey,sample_name] = tmp_vars[j,"TUMOR_MAF"]; VAF[ukey,"N"] = tmp_vars[j,"NORMAL_MAF"]
	DEPTH[ukey,sample_name] = tmp_vars[j,"TUMOR_DP"]; DEPTH[ukey,"N"] = tmp_vars[j,"NORMAL_DP"]
	LOH[ukey,sample_name] = tmp_vars[j,"facetsLOHCall"]
	CCF[ukey,sample_name] = tmp_vars[j,"ccf"]
	Clonal_Status[ukey,sample_name] = tmp_vars[j,"clonalStatus"]
}
vars = cbind(vars, VAF, DEPTH, LOH, CCF, Clonal_Status)
write.table(vars, file=paste0("sufam_multisample/", opt$patient, ".tsv"), col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

