# run expands for determining tumor ploidy

include modules/Makefile.inc

LOGDIR = log/expands.$(NOW)

SHELL = modules/scripts/Rshell
.SHELLFLAGS = -s -m $(MEM) -n $(@F) -l $(LOGDIR) -e 

.ONESHELL:
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: all

MEM = 20G

all : $(foreach pair,$(SAMPLE_PAIRS),expands/rdata/$(pair).cbs_snv.Rdata)

expands/rdata/%.cbs.Rdata : varscan/copycall/%.copycall
	cn <- read.table("$<", header=T, as.is=T)
	keep <- which(cn[,1] %in% c(1:22, "X"))
	if (length(rm) > 0) { cn <- cn[keep,]}
	cn[which(cn[,1]=="X"),1] <- 23
	cn[,1] <- as.numeric(cn[,1])
	cn <- cn[order(cn[,1], cn[,2]),]
	cn <- cbind(name = paste(cn[,1], cn[,2], sep="_"), cn[,c(1:3,7)])
	cgh <- make_cghRaw(cn)
	normalized <- normalize(cgh)
	segmented <- segmentData(normalized, relSDlong=2, undo.splits="sdundo", undo.SD=1.5)
	calls <- CGHcall(segmented, nclass=3)
	excalls <- ExpandCGHcall(calls, segmented)
	cbs <- with(fData(excalls), data.frame(chr = as.character(Chromosome[calls[[5]][,"wm"]]), startpos = Start[calls[[5]][,"wm"]], endpos = End[calls[[5]][,"wmend"]], CN_Estimate = 2^calls[[5]][, "smwh"], stringsAsFactors = F))
	cbs <- transform(cbs, segmentLength = endpos - startpos)


expands/rdata/%.snv.Rdata : mutect/tables/%.mutect.txt
	library(expands)
	snv <- read.table("$<", header = T, sep = "\t", stringsAsFactors = F)
	snv <- subset(snv, judgement != "REJECT")
	colnames(snv)[1:2] <- c("chr", "startpos")
	snv <- subset(snv, select = 'chr', 'startpos')
	snv$$chr <- as.integer(snv$$chr)
	snv <- as.matrix(snv[!is.na(snv$$chr), ])
	dir.create("$(@D)", recursive = T)
	dm <- assignQuantityToMutation(snv, cbs, "CN_Estimate")
	save(snv, cbs, file = "$@")



