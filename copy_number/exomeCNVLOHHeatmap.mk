include ~/share/modules/Makefile.inc

LOGDIR = log/exome_cnv_loh_heatmap.$(NOW)

SHELL = $(HOME)/share/scripts/Rshell
.SHELLFLAGS = -m $(MEM) -p $(PE) -n $(@F) -l $(LOGDIR) -e 

.ONESHELL:
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: all

MEM := 2G
PE := 1

all : exomecnv/lohheatmap.png

exomecnv/loh/lohmat.Rdata : $(foreach pair,$(SAMPLE_PAIRS),exomecnv/loh/$(pair).loh.txt)
	lohFiles <- unlist(strsplit("$^", " "))
	lohNames <- sub(".*/", "", sub("\\..*", "", lohFiles))
	suppressPackageStartupMessages(library("rtracklayer"));
	targets <- import('$(TARGETS_FILE)');
	for (i in 1:length(lohFiles)) {
		lohFile <- lohFiles[i]
		lohName <- lohNames[i]
		s <- read.delim(lohFile, header = T, as.is = T)
		lohGR <- GRanges(seqnames = sub('chr', '', s[, "chr"],), ranges = IRanges(start = s[, "position.start"], end = s[, "position.end"]), loh = s[, "LOH"])
		x <- suppressWarnings(findOverlaps(targets, lohGR))
		mcols(targets)[queryHits(x), lohName] <- lohGR[subjectHits(x)]$$loh
	}
	names(targets) <- paste(seqnames(targets), start(targets), sep="_")
	lohmat <- as.matrix(mcols(targets))
	rownames(lohmat) <- names(targets)
	lohmat[lohmat] <- 1
	lohmat[which(!lohmat | is.na(lohmat))] <- 0
	dir.create('$(@D)', showWarnings = F)
	save(lohmat, file = "$@")

exomecnv/lohheatmap.png : exomecnv/loh/lohmat.Rdata
	load("$<")
	suppressPackageStartupMessages(library("RColorBrewer"));
	suppressPackageStartupMessages(library("gplots"));
	cols <- c(brewer.pal(8, "Dark2"), brewer.pal(8, "Set1"), brewer.pal(8, "Set2"))
	chr <- unlist(lapply(rownames(lohmat), function(x) {strsplit(x, split="_", fixed=T)[[1]][1]}))
	dir.create('$(@D)', showWarnings = F)
	png("$@", height=600, width=1200, type="cairo")
	heatmap.2(t(lohmat), trace="none", scale = 'none', Colv = NA, col=c("white", "red"), margin=c(5,15), labCol="", ColSideColors=cols[as.integer(as.factor(chr))], cexCol=1.4, dendrogram = 'row', key = F)
	null <- dev.off()

