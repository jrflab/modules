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

gistic/varscanmat.Rdata : MEM := 2G
gistic/varscanmat.Rdata : PE := 8
gistic/varscanmat.Rdata : $(foreach pair,$(SAMPLE_PAIRS),varscan/segment/$(pair).seg.txt)
	segFiles <- unlist(strsplit("$^", " "));
	segNames <- sub(".*/", "", sub("\\..*", "", segFiles))
	suppressPackageStartupMessages(library("rtracklayer"));
	suppressPackageStartupMessages(library("snow"));
	targets <- import('$(TARGETS_FILE)')
	cl <- makeSOCKcluster(8)
	results <- parLapply(cl, segFiles, function(file, targets) {
		s <- read.delim(file, header=T, as.is=T)
		s$$Chromosome[s$$Chromosome==23] <- "X"
		s$$Chromosome[s$$Chromosome==24] <- "Y"
		rr <- rle(paste(s$$Chromosome, s$$log2_ratio_seg, sep="_"))
		ends <- cumsum(rr$$lengths)
		starts <- c(1, ends[-length(ends)]+1)
		rr$$values <- unlist(lapply(rr$$values, function(x) { strsplit(x, split="_")[[1]][2]}))
		segName <- sub(".*/", "", sub("\\..*", "", file))
		tab <- cbind(segName, s$$Chromosome[starts], s$$Start[starts], s$$End[ends], rr$$lengths, rr$$values)
		tab <- as.data.frame(tab, stringsAsFactors=F)
		colnames(tab) <- c("Sample", "Chromosome", "Start", "End", "Num.markers", "segmented")
		thissample <- rep(NA, length(targets))
		for (j in 1:nrow(tab)) {
			thissample[which(seqnames(targets) == tab$$Chromosome[j] & as.numeric(start(targets)) >= as.numeric(tab$$Start[j]) & 
				as.numeric(start(targets)) <= as.numeric(tab$$End[j]))] <- as.numeric(tab$$segmented[j])
			if (j %% 100 == 0) { print (j); gc()}
		}
		as.numeric(thissample)
	}, targets)
	stopCluster(cl)
	varscanmat <- do.call("cbind", results)
	colnames(varscanmat) <- segNames
	rownames(varscanmat) <- paste(seqnames(targets), start(targets), sep="_")
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
		tab <- cbind(colnames(varscanmat)[i], 
		markers[starts,1], markers[starts,2], markers[ends,2], 
		as.vector(rr$$lengths), as.vector(rr$$values))
		tab <- as.data.frame(tab, stringsAsFactors=F)
		colnames(tab) <- c("Sample", "Chromosome", "Start", "End", "Num.markers", "segmented")
		seg <- rbind(seg, tab)
	}
	dir.create('$(@D)', showWarnings = F)
	write.table(seg, file="$@", sep="\t", row.names=F, col.names=F, quote=F)


gistic/lohmat.Rdata : MEM := 2G
gistic/lohmat.Rdata : PE := 8
gistic/lohmat.Rdata : $(foreach pair,$(SAMPLE_PAIRS),exomecnv/loh/$(pair).loh.txt)
	lohFiles <- unlist(strsplit("$^", " "))
	lohNames <- sub(".*/", "", sub("\\..*", "", lohFiles))
	suppressPackageStartupMessages(library("rtracklayer"));
	suppressPackageStartupMessages(library("snow"));
	targets <- import('$(TARGETS_FILE)');
	cl <- makeSOCKcluster(8)
	results <- parLapply(cl, lohFiles, function(file, targets) {
		s <- read.delim(file, header=T, as.is=T)
		s <- s[which(s$$LOH=="TRUE"),]
		s <- s[which(unlist(apply(s[,2:3],1,function(x){x[1]!=x[2]}))),]
		s$$chr <- gsub("chr", "", s$$chr)
		thissample <- rep(FALSE, length(targets))
		for (j in 1:nrow(s)) {
			x <- which(seqnames(targets) == s$$chr[j] & as.numeric(start(targets)) >= as.numeric(s$$position.start[j]) & as.numeric(start(targets)) <= as.numeric(s$$position.end[j]))
			if (length(x) > 0) {
				thissample[x] <- s$$LOH[j]
			}
		}
		as.numeric(thissample)
	}, targets)
	stopCluster(cl)
	lohmat <- do.call("cbind", results)
	colnames(lohmat) <- lohNames
	rownames(lohmat) <- paste(seqnames(targets), start(targets), sep="_")
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
	suppressPackageStartupMessages(library("RColorBrewer"));
	suppressPackageStartupMessages(library("gplots"));
	cols <- c(brewer.pal(8, "Dark2"), brewer.pal(8, "Set1"), brewer.pal(8, "Set2"))
	chr <- unlist(lapply(rownames(lohmat), function(x) {strsplit(x, split="_", fixed=T)[[1]][1]))
	dir.create('$(@D)', showWarnings = F)
	png("$@", height=1200, width=600, type="cairo")
	heatmap.2(lohmat, trace="none", Rowv=F, col=c("white", "red"), margin=c(12,5), labRow="", RowSideColors=col[chr], cexCol=1.4)
	null <- dev.off()
