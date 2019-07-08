include modules/Makefile.inc
include modules/vcf_tools/vcftools.mk

LOGDIR = log/absoluteSeq.$(NOW)
MEM := 4G
PE := 1
SHELL = modules/scripts/Rshell
.SHELLFLAGS = $(if $(findstring true,$(USE_CLUSTER)),-s) -m $(MEM) -p $(PE) -n $(@F) -l $(LOGDIR) -e 

.ONESHELL:
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: absolute absolute_rdata absolute_reviewed absolute_tables

PRIMARY_DISEASE ?= breast
PLATFORM ?= Illumina_WES
PYTHON_ENV ?= /ifs/e63data/reis-filho/usr/anaconda-envs/pyenv27-chasm

absolute : absolute/review/all.PP-calls_tab.txt absolute_rdata
absolute_rdata : $(foreach pair,$(SAMPLE_PAIRS),absolute/results/$(pair).ABSOLUTE.RData)
absolute_reviewed : absolute/reviewed/all.seq.ABSOLUTE.table.txt
absolute_tables : $(foreach pair,$(SAMPLE_PAIRS),absolute/tables/$(pair).absolute.txt)

USE_TITAN_COPYNUM ?= false
USE_TITAN_ESTIMATES ?= false
TITAN_RESULTS_DIR ?= titan/optclust_results_w10000_p2
TITAN_ESTIMATE_FILE ?= $(TITAN_RESULTS_DIR)/titan_summary.txt

USE_ONCOSCAN_COPYNUM ?= false

USE_FACETS_COPYNUM ?= true

define LIB_INIT
for (lib in c("dplyr","stringr","readr","ABSOLUTE")) {
    suppressPackageStartupMessages(library(lib, character.only=TRUE))
}
endef

absolute/tables/%.somatic.txt : tables/%.mutect.tab.txt tables/%.strelka_varscan_indels.tab.txt

	$(R_INIT)
	$(LIB_INIT)

	tn <- unlist(strsplit("$*", '_'))
	snvs <- read.delim("$<", header = T, sep = '\t', comment.char = '', as.is = T, check.names = F)
	snvs.tref <- c()
	snvs.talt <- c()

	if (nrow(snvs) > 0) {
		snvs.tref <- sapply(strsplit(snvs[[paste(tn[1], ".AD", sep = '')]], ','), function (x) x[1])
		snvs.talt <- sapply(strsplit(snvs[[paste(tn[1], ".AD", sep = '')]], ','), function (x) x[2])
	}

	indels <- read.delim("$(<<)", header = T, sep = '\t', comment.char = '', as.is = T, check.names = F)
	indels.tref <- c()
	indels.talt <- c()
	if (nrow(indels) > 0) {
		indels.tref <- sapply(strsplit(indels[[paste(tn[1], ".AD", sep = '')]], ','), function (x) x[1])
		indels.talt <- sapply(strsplit(indels[[paste(tn[1], ".AD", sep = '')]], ','), function (x) x[2])
	}

	chr       <- c(snvs[["CHROM"]], indels[["CHROM"]])
	chr       <- as.integer(sub('X', '23', chr))
	ref       <- c(snvs[["REF"]], indels[["REF"]])
	alt       <- c(snvs[["ALT"]], indels[["ALT"]])
	genes     <- c(snvs[['ANN[*].GENE']], indels[['ANN[*].GENE']])
	pos       <- c(snvs[["POS"]], indels[["POS"]])
	tRefCount <- c(snvs.tref, indels.tref)
	tAltCount <- c(snvs.talt, indels.talt)
	type      <- c(rep('snv', nrow(snvs)), rep('indel', nrow(indels)))

	Data <-
		data_frame( Sample           = "$*",
					Gene             = sub('\\|.*', '', genes),
					Chromosome       = chr,
					Position         = pos,
					Ref              = ref,
					Alt              = alt,
					tRefCount        = tRefCount,
					tAltCount        = tAltCount,
					Type             = type ) %>%
		filter(Gene!=".")

	x <- with(Data, paste(Chromosome, Position, sep=":"))
	Data <- subset(Data, !duplicated(x))
	write_tsv(Data, path="$@")


absolute/maf/%.maf.txt : absolute/tables/%.somatic.txt
	$(R_INIT)
	$(LIB_INIT)
	X <- read.delim("$(<)", header = T, sep = '\t', comment.char = '', as.is = T, check.names = F)
	Data <- with(X, data.frame(Tumor_Sample_Barcode = Sample, Hugo_Symbol = Gene, t_ref_count = tRefCount, t_alt_count = tAltCount, dbSNP_Val_Status = "validated", Chromosome = Chromosome, Start_position = Position, stringsAsFactors = F))
	write.table(Data, file = "$@", sep = '\t', quote = F, row.names = F)

ifeq ($(USE_TITAN_COPYNUM),true)
absolute/segment/%.seg.txt : $(TITAN_RESULTS_DIR)/%.z*.titan.seg
	$(R_INIT)
	$(LIB_INIT)
	X <- read.delim("$<", header = T, sep = '\t', as.is = T)[,-1]
	colnames(X) <- c('Chromosome', 'Start', 'End', 'Num_Probes', 'Segment_Mean')
	X[,1] <- sub('X', '23', X[,1])
	write.table(X, file = "$@", row.names = F, quote = F, sep = '\t')
else ifeq ($(USE_ONCOSCAN_COPYNUM),true)
absolute/segment/%.seg.txt : oncoscan/%.probes.txt oncoscan/%.segments.txt
	$(R_INIT)
	r <- system('unset PYTHONPATH && source $(PYTHON_ENV)/bin/activate $(PYTHON_ENV) && python <<EOF
	import pandas as pd
	import math
	p = pd.read_csv("$<", sep="\t")
	s = pd.read_csv("$(<<)", sep="\t")
	s["Num_Probes"] = s.apply(lambda x: len(p[(p.Start >= x.Start) & (p.End <= x.End)]), axis=1)
	s["Segment_Mean"] = s.Value
	del s["Value"]
	# Center segment mean around zero
	s.Segment_Mean = (2 - (s.Segment_Mean.apply(lambda x: 2**x * 2).mean() - s.Segment_Mean.apply(lambda x: 2**x * 2))).apply(lambda x: max(0.1, x)).apply(lambda x: math.log(x/2, 2))
	# Rename chromosomes
	s.Chromosome = s.Chromosome.str.replace("chrX","23").str.replace("chrY", "24").str.replace("chr","")
	s.to_csv("$@", index=False, sep="\t")
	EOF')
	q(status=r, save="no")
else ifeq ($(USE_FACETS_COPYNUM),true)
absolute/segment/%.seg.txt : facets/cncf/%.cncf.txt
	$(R_INIT)
	$(LIB_INIT)
	cncf <-
		read.delim("$<", stringsAsFactors=F, header=T, sep="\t") %>%
		select(ID, Chromosome=chrom, Start=loc.start, End=loc.end, seg, Num_Probes=num.mark, nhet, Segment_Mean=cnlr.median, mafR, segclust, mafR.clust, cnlr.median.clust, cf.em, tcn.em, lcn.em)

	cncf %>% write_tsv(path="$@")
else
absolute/segment/%.seg.txt : varscan/segment/%.collapsed_seg.txt
	$(R_INIT)
	$(LIB_INIT)
	X <- read.delim("$<", header = T, sep = '\t', as.is = T, check.names = F)
	colnames(X) <- c('Chromosome', 'Start', 'End', 'Num_Probes', 'Segment_Mean')
	write.table(X, file = "$@", row.names = F, quote = F, sep = '\t')
endif
	
ifeq ($(USE_TITAN_ESTIMATES),true)
$(info absolute/results/%.ABSOLUTE.RData : absolute/segment/%.seg.txt absolute/maf/%.maf.txt $(TITAN_ESTIMATE_FILE))
absolute/results/%.ABSOLUTE.RData : absolute/segment/%.seg.txt absolute/maf/%.maf.txt $(TITAN_ESTIMATE_FILE)
	$(R_INIT)
	$(LIB_INIT)
	titanResults <- read.delim("$(<<<)", row.names = 1, header = T, check.names = F)
	sigma.p <- 0
	max.sigma.h <- 0.07
	avgTumorPloidyEst <- titanResults["$*", "avgTumorPloidyEst"]
	min.ploidy <- max(0.95, avgTumorPloidyEst - 1)
	max.ploidy <- min(7, avgTumorPloidyEst + 1)
	primary.disease <- "$(PRIMARY_DISEASE)"
	sample.name <- "$*"
	platform <- "$(PLATFORM)"
	max.as.seg.count <- 3500
	copynum.type <- "total"
	max.neg.genome <- 0
	max.non.clonal <- 0
	min.mut.af <- 0
	seg.dat.fn <- "$<"
	results.dir <- "$(@D)"
	output.fn.base = "$*"
	maf.fn = "$(<<)"
	RunAbsolute(seg.dat.fn, sigma.p, max.sigma.h, min.ploidy, max.ploidy, primary.disease, platform,sample.name, results.dir, max.as.seg.count, max.non.clonal, max.neg.genome, copynum.type, maf.fn = maf.fn, min.mut.af = min.mut.af, output.fn.base = output.fn.base, verbose = T)
else
$(info absolute/results/%.ABSOLUTE.RData : absolute/segment/%.seg.txt absolute/maf/%.maf.txt)
absolute/results/%.ABSOLUTE.RData : absolute/segment/%.seg.txt absolute/maf/%.maf.txt
	$(R_INIT)
	$(LIB_INIT)
	sigma.p <- 0
	max.sigma.h <- 0.07
	min.ploidy <- 0.95 
	max.ploidy <- 7
	primary.disease <- "$(PRIMARY_DISEASE)"
	sample.name <- "$*"
	platform <- "$(PLATFORM)"
	max.as.seg.count <- 3500
	copynum.type <- "total"
	max.neg.genome <- 0
	max.non.clonal <- 0
	min.mut.af <- 0
	seg.dat.fn <- "$<"
	results.dir <- "$(@D)"
	output.fn.base = "$*"
	maf.fn = "$(<<)"
	RunAbsolute(seg.dat.fn, sigma.p, max.sigma.h, min.ploidy, max.ploidy, primary.disease, platform,sample.name, results.dir, max.as.seg.count, max.non.clonal, max.neg.genome, copynum.type, maf.fn = maf.fn, min.mut.af = min.mut.af, output.fn.base = output.fn.base, verbose = T)
endif

absolute/review/%.PP-calls_tab.txt absolute/review/%.PP-modes.data.RData : $(foreach pair,$(SAMPLE_PAIRS),absolute/results/$(pair).ABSOLUTE.RData)
	$(R_INIT)
	$(LIB_INIT)
	absolute.files <- qw("$^")
	indv.results.dir <- "$(@D)"
	copynum.type <- "total"
	CreateReviewObject(obj.name = "$*", absolute.files, indv.results.dir, copynum.type, plot.modes = T, verbose = T)

absolute/reviewed/all.seq.ABSOLUTE.table.txt : absolute/review/all.PP-calls_tab.reviewed.txt absolute/review/all.PP-modes.data.RData
	$(R_INIT)
	$(LIB_INIT)
	ExtractReviewedResults("$<", 'seq', "$(<<)", "absolute", "all", verbose = T, copy_num_type = "total")

absolute/tables/%.absolute.txt : absolute/reviewed/all.seq.ABSOLUTE.table.txt absolute/tables/%.somatic.txt
	$(R_INIT)
	$(LIB_INIT)
	fn <- 'absolute/reviewed/SEG_MAF/$*_ABS_MAF.txt'
	absData <- read.delim(fn, sep = '\t', header = T, as.is = T, check.names = F)
	tn <- absData[['sample']][1]
	somatic <- read.delim("$(<<)", header = T, sep = '\t', comment.char = '', as.is = T)
	absData[, "refSeq"] <- somatic[, "Ref"]
	absData[, "altSeq"] <- somatic[ ,"Alt"]
	absData[, "type"] = somatic[, "Type"]
	absData[, "mut"] <- paste(absData[, "refSeq"], absData[, "altSeq"], sep = ">")
	absData <- transform(absData, somaticSubclonal = Pr_subclonal > 0.5 & Pr_somatic > 0.95, somaticClonal = Pr_somatic_clonal > 0.5)
	absData <- subset(absData, somaticSubclonal | somaticClonal)
	write.table(absData, file = "$@", quote = F, row.names = F, sep = '\t')
