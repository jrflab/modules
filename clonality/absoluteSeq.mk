include modules/Makefile.inc
include modules/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR = log/absoluteSeq.$(NOW)
MEM := 2G
PE := 1
SHELL = scripts/Rshell
.SHELLFLAGS = -s -m $(MEM) -p $(PE) -n $(@F) -l $(LOGDIR) -e 

.ONESHELL:
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: absolute absolute_rdata absolute_reviewed absolute_tables

PRIMARY_DISEASE ?= breast
PLATFORM ?= Illumina_WES

absolute : absolute/review/all.PP-calls_tab.txt absolute_rdata
absolute_rdata : $(foreach pair,$(SAMPLE_PAIRS),absolute/results/$(pair).ABSOLUTE.RData)
absolute_reviewed : absolute/reviewed/all.seq.ABSOLUTE.table.txt
absolute_tables : $(foreach pair,$(SAMPLE_PAIRS),absolute/tables/$(pair).absolute.txt)

USE_TITAN_COPYNUM ?= true
USE_TITAN_ESTIMATES ?= false
TITAN_RESULTS_DIR ?= titan/optclust_results_w1000_p2
TITAN_ESTIMATE_FILE ?= $(TITAN_RESULTS_DIR)/titan_summary.txt

define LIB_INIT
library(ABSOLUTE)
endef


absolute/tables/%.somatic.txt : tables/%.$(call FILTER_SUFFIX,mutect).tab.txt tables/%.$(call FILTER_SUFFIX,strelka_indels).tab.txt tables/%.$(call FILTER_SUFFIX,scalpel_indels).tab.pass.txt
	$(R_INIT)
	$(LIB_INIT)
	tn <- unlist(strsplit("$*", '_'))
	snvs <- read.table("$<", header = T, sep = '\t', comment.char = '', as.is = T)
	snvs.tref <- c()
	snvs.talt <- c()
	if (nrow(snvs) > 0) {
		snvs.tref <- sapply(strsplit(snvs[[paste(tn[1], ".AD", sep = '')]], ','), function (x) x[1])
		snvs.talt <- sapply(strsplit(snvs[[paste(tn[1], ".AD", sep = '')]], ','), function (x) x[2])
	}
	indels <- read.table("$(<<)", header = T, sep = '\t', comment.char = '', as.is = T)
	indels.tref <- c()
	indels.talt <- c()
	if (nrow(indels) > 0) {
		indels.talt <- as.integer(sapply(strsplit(indels[["TUMOR.TIR"]], ','), function (x) x[1]))
		indels.tref <- indels[["TUMOR.DP"]] - indels.talt
	}
	indels2 <- read.table("$(<<<)", header = T, sep = '\t', comment.char = '', as.is = T)
	indels2.tref <- c()
	indels2.talt <- c()
	if (nrow(indels2) > 0) {
		indels2.tref <- sapply(strsplit(indels2[[paste(tn[1], ".AD", sep = '')]], ','), function (x) x[1])
		indels2.talt <- sapply(strsplit(indels2[[paste(tn[1], ".AD", sep = '')]], ','), function (x) x[2])
	}
	chr <- c(snvs[["X.CHROM"]], indels[["X.CHROM"]], indels2[["X.CHROM"]])
	chr <- as.integer(sub('X', '23', chr))
	ref <- c(snvs[["REF"]], indels[["REF"]], indels2[["REF"]])
	alt <- c(snvs[["ALT"]], indels[["ALT"]], indels2[["ALT"]])
	genes <- c(snvs[['EFF....GENE']], indels[['EFF....GENE']], indels2[['EFF....GENE']])
	pos <- c(snvs[["POS"]], indels[["POS"]], indels2[["POS"]])
	tRefCount <- c(snvs.tref, indels.tref, indels2.tref)
	tAltCount = c(snvs.talt, indels.talt, indels2.talt)
	type = c(rep('snv', nrow(snvs)), rep('indel', nrow(indels)), rep('indel', nrow(indels2)))
	Data <- data.frame(Sample = "$*", Gene = genes, Chromosome = chr, Position = pos, Ref = ref, Alt = alt, tRefCount = tRefCount, tAltCount = tAltCount, Type = type, stringsAsFactors = F)
	Data[["Gene"]] <- sub('\\|.*', '', Data[["Gene"]])
	Data <- subset(Data, Gene != ".")
	x <- with(Data, paste(Chromosome, Position, sep = ":"))
	Data <- subset(Data, !duplicated(x))
	write.table(Data, file = "$@", sep = '\t', quote = F, row.names = F)


absolute/maf/%.maf.txt : absolute/tables/%.somatic.txt
	$(R_INIT)
	$(LIB_INIT)
	X <- read.table("$(<)", header = T, sep = '\t', comment.char = '', as.is = T)
	Data <- with(X, data.frame(Tumor_Sample_Barcode = Sample, Hugo_Symbol = Gene, t_ref_count = tRefCount, t_alt_count = tAltCount, dbSNP_Val_Status = "validated", Chromosome = Chromosome, Start_position = Position, stringsAsFactors = F))
	write.table(Data, file = "$@", sep = '\t', quote = F, row.names = F)

ifeq ($(USE_TITAN_COPYNUM),true)
absolute/segment/%.seg.txt : $(TITAN_RESULTS_DIR)/%.z*.titan.seg
	$(R_INIT)
	$(LIB_INIT)
	X <- read.table("$<", header = T, sep = '\t', as.is = T)[,-1]
	colnames(X) <- c('Chromosome', 'Start', 'End', 'Num_Probes', 'Segment_Mean')
	X[,1] <- sub('X', '23', X[,1])
	write.table(X, file = "$@", row.names = F, quote = F, sep = '\t')
else
absolute/segment/%.seg.txt : varscan/segment/%.collapsed_seg.txt
	$(R_INIT)
	$(LIB_INIT)
	X <- read.table("$<", header = T, sep = '\t', as.is = T)
	colnames(X) <- c('Chromosome', 'Start', 'End', 'Num_Probes', 'Segment_Mean')
	write.table(X, file = "$@", row.names = F, quote = F, sep = '\t')
endif
	
ifeq ($(USE_TITAN_ESTIMATES),true)
absolute/results/%.ABSOLUTE.RData : absolute/segment/%.seg.txt absolute/maf/%.maf.txt $(TITAN_ESTIMATE_FILE)
	$(R_INIT)
	$(LIB_INIT)
	titanResults <- read.table("$(<<<)", row.names = 1, header = T)
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
	absData <- read.table(fn, sep = '\t', header = T, as.is = T)
	tn <- absData[['sample']][1]
	somatic <- read.table("$(<<)", header = T, sep = '\t', comment.char = '', as.is = T)
	absData[, "refSeq"] <- somatic[, "Ref"]
	absData[, "altSeq"] <- somatic[ ,"Alt"]
	absData[, "type"] = somatic[, "Type"]
	absData[, "mut"] <- paste(absData[, "refSeq"], absData[, "altSeq"], sep = ">")
	absData <- transform(absData, somaticSubclonal = Pr_subclonal > 0.5 & Pr_somatic > 0.95, somaticClonal = Pr_somatic_clonal > 0.5)
	absData <- subset(absData, somaticSubclonal | somaticClonal)
	write.table(absData, file = "$@", quote = F, row.names = F, sep = '\t')
