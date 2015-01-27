include ~/share/modules/Makefile.inc
include ~/share/modules/variant_callers/somatic/mutect.inc
include ~/share/modules/variant_callers/somatic/strelka.inc

LOGDIR = log/absoluteSeq.$(NOW)
MEM := 2G
PE := 1
SHELL = $(HOME)/share/scripts/Rshell
.SHELLFLAGS = -s -m $(MEM) -p $(PE) -n $(@F) -l $(LOGDIR) -e 

.ONESHELL:
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: all results

PRIMARY_DISEASE ?= breast
PLATFORM ?= Illumina_WES

all : absolute/review/all.PP-calls_tab.txt results
results : $(foreach pair,$(SAMPLE_PAIRS),absolute/results/$(pair).ABSOLUTE.RData)

USE_TITAN ?= false
TITAN_RESULTS_DIR ?= titan/optclust_results_w1000_p2
TITAN_ESTIMATE_FILE ?= $(TITAN_RESULTS_DIR)/titan_summary.txt

define LIB_INIT
library(ABSOLUTE)
endef

absolute/maf/%.maf.txt : tables/%.mutect.$(MUTECT_FILTER_SUFFIX).tab.txt tables/%.strelka_indels.$(STRELKA_FILTER_SUFFIX.strelka_indels).tab.txt
	$(R_INIT)
	$(LIB_INIT)
	tn <- unlist(strsplit("$*", '_'))
	snvs <- read.table("$<", header = T, sep = '\t', comment.char = '', as.is = T)
	snvs.tref <- sapply(strsplit(snvs[[paste(tn[1], ".AD", sep = '')]], ','), function (x) x[1])
	snvs.talt <- sapply(strsplit(snvs[[paste(tn[1], ".AD", sep = '')]], ','), function (x) x[2])
	indels <- read.table("$(<<)", header = T, sep = '\t', comment.char = '', as.is = T)
	indels.talt <- as.integer(sapply(strsplit(indels[["TUMOR.TIR"]], ','), function (x) x[1]))
	indels.tref <- indels[["TUMOR.DP"]] - indels.talt
	chr <- c(snvs[["X.CHROM"]], indels[["X.CHROM"]])
	chr <- as.integer(sub('X', '23', chr))
	Data <- data.frame(Tumor_Sample_Barcode = tn[1], Hugo_Symbol = c(snvs[['EFF....GENE']], indels[['EFF....GENE']]), t_ref_count = c(snvs.tref, indels.tref), t_alt_count = c(snvs.talt, indels.talt), dbSNP_Val_Status = "validated", Chromosome = chr, Start_position = c(snvs[["POS"]], indels[["POS"]]), stringsAsFactors = F)
	Data <- subset(Data, Hugo_Symbol != ".")
	write.table(Data, file = "$@", sep = '\t', quote = F, row.names = F)

ifeq ($(USE_TITAN),true)
absolute/segment/%.seg.txt : $(TITAN_RESULTS_DIR)/%.z*.titan.seg
	$(R_INIT)
	$(LIB_INIT)
	X <- read.table("$<", header = T, sep = '\t')[,-1]
	colnames(X) <- c('Chromosome', 'Start', 'End', 'Num_Probes', 'Segment_Mean')
	write.table(X, file = "$@", row.names = F, quote = F, sep = '\t')
	
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
absolute/segment/%.seg.txt : varscan/segment/%.collapsed_seg.txt
	$(R_INIT)
	$(LIB_INIT)
	X <- read.table("$<", header = T, sep = '\t')
	colnames(X) <- c('Chromosome', 'Start', 'End', 'Num_Probes', 'Segment_Mean')
	write.table(X, file = "$@", row.names = F, quote = F, sep = '\t')

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

absolute/reviewed/all.seq.ABSOLUTE.table.txt : absolute/review/all.PP-calls_tab.txt absolute/review/all.PP-modes.data.RData
	$(R_INIT)
	$(LIB_INIT)
	ExtractReviewedResults("$<", 'seq', "$(<<)", "absolute", "all", verbose = T, copy_num_type = "total")
