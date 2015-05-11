# run bicseq for cnvs

include modules/Makefile.inc

LOGDIR = log/bicseq.$(NOW)

SHELL = modules/scripts/Rshell
.SHELLFLAGS = -s -m $(MEM) -n $(@F) -l $(LOGDIR) -e 

.ONESHELL:
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: all

MEM = 9G

all : $(foreach pair,$(SAMPLE_PAIRS),expands/rdata/$(pair).cbs_snv.Rdata)



expands/rdata/%.cbs_snv.Rdata : mutect/tables/%.mutect.txt varscan/segment/%.varscan2copynumber.txt
	library(expands)
	snv <- read.table("$<", header = T, sep = "\t")
	snv <- subset(snv, judgement != "REJECT")
	colnames(snv)[1:2] <- c("chr", "startpos")
	cbs <- read.table("$(<<)", header = T, sep = "\t")
	cbs <- transform(cbs, CN_estimate = 2^Segmented)
	colnames(cbs)[c(2,3,4)] <- c("chr", "startpos", "endpos")
	dir.create("$(@D)", recursive = T)
	save(snv, cbs, file = "$@")



