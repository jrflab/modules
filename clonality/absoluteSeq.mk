include ~/share/modules/Makefile.inc

LOGDIR = log/absoluteSeq.$(NOW)
MEM := 2G
PE := 1
SHELL = $(HOME)/share/scripts/Rshell
.SHELLFLAGS = -s -m $(MEM) -p $(PE) -n $(@F) -l $(LOGDIR) -e 

.ONESHELL:
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: all

PRIMARY_DISEASE ?= breast
PLATFORM ?= Illumina_WES

all : absolute/reviewed/all.segtab.txt

define LIB_INIT
library(ABSOLUTE)
endef

absolute/segment/%.seg.txt : varscan/segment/%.collapsed_seg.txt
	$(R_INIT)
	$(LIB_INIT)
	X <- read.table("$<", header = T, sep = '\t')
	colnames(X) <- c('Chromosome', 'Start', 'End', 'Num_Probes', 'Segment_Mean')
	write.table(X, file = "$@", row.names = F, quote = F, sep = '\t')

absolute/results/%.ABSOLUTE.RData : absolute/segment/%.seg.txt
	$(R_INIT)
	$(LIB_INIT)
	sigma.p <- 0
	max.sigma.h <- 0.07
	min.ploidy <- 0.95 
	max.ploidy <- 7
	primary.disease <- "$(PRIMARY_DISEASE)"
	sample.name <- "$*"
	platform <- "$(PLATFORM)"
	max.as.seg.count <- 1500
	copynum.type <- "total"
	max.neg.genome <- 0
	max.non.clonal <- 0
	min.mut.af <- 0
	seg.dat.fn <- "$<"
	results.dir <- "$(@D)"
	output.fn.base = "$*"
	RunAbsolute(seg.dat.fn, sigma.p, max.sigma.h, min.ploidy, max.ploidy, primary.disease, platform,sample.name, results.dir, max.as.seg.count, max.non.clonal, max.neg.genome, copynum.type, maf.fn = NULL, min.mut.af = NULL, output.fn.base = output.fn.base, verbose = T)

absolute/review/%.PP-calls_tab.txt absolute/review/%.PP-modes.data.RData : $(foreach pair,$(SAMPLE_PAIRS),absolute/results/$(pair).ABSOLUTE.RData)
	$(R_INIT)
	$(LIB_INIT)
	absolute.files <- qw("$^")
	indv.results.dir <- "$(@D)"
	copynum.type <- "total"
	CreateReviewObject(obj.name = "$*", absolute.files, indv.results.dir, copynum.type, plot.modes = T, verbose = T)

absolute/reviewed/all.segtab.txt : absolute/review/all.PP-calls_tab.txt absolute/review/all.PP-modes.data.RData
	$(R_INIT)
	$(LIB_INIT)
	ExtractReviewedResults("$<", 'absolute-workflow', "$(<<)", "absolute", "all", verbose = T, copy_num_type = "total")
