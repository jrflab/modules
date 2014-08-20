include ~/share/modules/Makefile.inc

LOGDIR = log/gistic.$(NOW)

SHELL = $(HOME)/share/scripts/Rshell
.SHELLFLAGS = -s -m $(MEM) -p $(PE) -n $(@F) -l $(LOGDIR) -e 

.ONESHELL:
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: all gistic_inputs

MEM := 2G
PE := 1

DGV_FILE = $(HOME)/share/reference/GRCh37_hg19_variants_2013-07-23.txt

CNV_SIZES = 100000 300000

all : gistic_inputs gistic/lohheatmap.png
gistic_inputs : gistic/markersfile.txt gistic/segmentationfile.txt $(foreach size,$(CNV_SIZES),gistic/cnv.$(size).txt)

gistic/varscanmat.Rdata : PE := 8
gistic/varscanmat.Rdata : MEM := 1G
gistic/varscanmat.Rdata : $(foreach pair,$(SAMPLE_PAIRS),varscan/segment/$(pair).seg.txt)
	suppressPackageStartupMessages(library("rtracklayer"));
	suppressPackageStartupMessages(library("foreach"));
	suppressPackageStartupMessages(library("doMC"));
	segFiles <- unlist(strsplit("$^", " "));
	segNames <- sub(".*/", "", sub("\\..*", "", segFiles))
	targets <- import('$(TARGETS_FILE)')
	names(targets) <- paste(seqnames(targets), start(targets), sep="_")
	registerDoMC(8)
	varscanmat <- foreach (i = 1:length(segFiles), .combine = 'cbind') %dopar% {
		segFile <- segFiles[i]
		segName <- segNames[i]
		s <- read.delim(segFile, header = T, as.is = T)
		s$$Chromosome[s$$Chromosome==23] <- "X"
		s$$Chromosome[s$$Chromosome==24] <- "Y"
		rr <- rle(paste(s$$Chromosome, s$$log2_ratio_seg, sep="_"))
		ends <- cumsum(rr$$lengths)
		starts <- c(1, ends[-length(ends)]+1)
		rr$$values <- unlist(lapply(rr$$values, function(x) { strsplit(x, split="_")[[1]][2]}))
		gr <- GRanges(seqnames = s$$Chromosome[starts], range = IRanges(start = s$$Start[starts], end = s$$End[ends]), segmented = as.numeric(rr$$values))
		x <- suppressWarnings(findOverlaps(targets, gr))
		#mcols(targets)[queryHits(x), segName] <-
		xx <- rep(NA, length(targets))
		xx[queryHits(x)] <- gr[subjectHits(x)]$$segmented
		xx
	}
	rownames(varscanmat) <- names(targets)
	colnames(varscanmat) <- segNames
	dir.create('$(@D)', showWarnings = F)
	save(varscanmat, file = "$@")


gistic/markersfile.txt : gistic/varscanmat.Rdata
	load('$<')
	markers <- matrix(unlist(lapply(rownames(varscanmat), function(x) { strsplit(x, split="_", fixed=T) })), ncol=2, byrow=2)
	rownames(markers) <- 1:nrow(markers)
	dir.create('$(@D)', showWarnings = F)
	write.table(markers, col.names=F, file="$@", sep="\t", quote=F, na="")

gistic/segmentationfile.txt : gistic/varscanmat.Rdata gistic/markersfile.txt
	load('$<')
	markers <- read.table("$(<<)", sep = "\t")
	seg <- matrix(nrow=0, ncol=6)
	colnames(seg) <- c("Sample", "Chromosome", "Start", "End", "Num.markers", "segmented")
	chr <- markers[,1]
	for (i in 1:ncol(varscanmat)) {
		notna <- which(!is.na(varscanmat[,i]))
		rr <- rle(paste(chr[notna], varscanmat[notna,i], sep="_"))
		ends <- cumsum(rr$$lengths)
		starts <- c(1, ends[-length(ends)]+1)
		rr$$values <- unlist(lapply(rr$$values, function(x) { strsplit(x, split="_")[[1]][2]}))
		tab <- cbind(colnames(varscanmat)[i], markers[starts,1], markers[starts,2], markers[ends,2], as.vector(rr$$lengths), as.vector(rr$$values))
		tab <- as.data.frame(tab, stringsAsFactors=F)
		colnames(tab) <- c("Sample", "Chromosome", "Start", "End", "Num.markers", "segmented")
		seg <- rbind(seg, tab)
	}
	dir.create('$(@D)', showWarnings = F)
	write.table(seg, file="$@", sep="\t", row.names=F, col.names=F, quote=F)


gistic/lohmat.Rdata : $(foreach pair,$(SAMPLE_PAIRS),exomecnv/loh/$(pair).loh.txt)
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

gistic/cnv.%.txt : gistic/markersfile.txt
	suppressPackageStartupMessages(library("GenomicRanges"));
	dgv <- read.delim("$(DGV_FILE)", as.is=T)
	dgv <- dgv[which(dgv[,5]=="CNV"),]
	dgv <- dgv[,1:4]
	dgv$$size = dgv[,4]-dgv[,3]+1
	dgv <-dgv[which(dgv$$size <= $*),]
	dgv <- dgv[which(dgv$$chr %in% 1:22),]
	markers <- read.delim("$<", as.is=T, header=F)
	dgvGR <- GRanges(seqnames = dgv$$chr, ranges = IRanges(start = dgv$$start, end = dgv$$end))
	markersGR <- GRanges(seqnames = markers[,2], ranges = IRanges(start = markers[,3], end = markers[,3]))
	markers <- cbind(markers, countOverlaps(markersGR, dgvGR))
	#markers <- cbind(markers, apply(markers, 1, function(x, dgv) {
	#	length(which(dgv$$chr==x[2] & dgv$$start <= x[3] & dgv$$end >= x[3])) }, dgv))
	cnv <- markers[which(markers[,4] > 0),]
	cnv <- cbind(cnv[,1], cnv[,1])
	dir.create('$(@D)', showWarnings = F)
	write.table(cnv, file = "$@", sep = "\t", row.names = F, col.names = F, quote = F, na = "")

gistic/lohheatmap.png : gistic/lohmat.Rdata
	load("$<")
	suppressPackageStartupMessages(library("RColorBrewer"));
	suppressPackageStartupMessages(library("gplots"));
	cols <- c(brewer.pal(8, "Dark2"), brewer.pal(8, "Set1"), brewer.pal(8, "Set2"))
	chr <- unlist(lapply(rownames(lohmat), function(x) {strsplit(x, split="_", fixed=T)[[1]][1]}))
	dir.create('$(@D)', showWarnings = F)
	png("$@", height=1200, width=600, type="cairo")
	heatmap.2(t(lohmat), trace="none", scale = 'none', Colv = NA, col=c("white", "red"), margin=c(5,15), labCol="", ColSideColors=cols[as.integer(as.factor(chr))], cexCol=1.4, dendrogram = 'row', key = F)
	null <- dev.off()
