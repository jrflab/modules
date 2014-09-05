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

export LD_LIBRARY_PATH = /home/limr/usr/MATLAB/v714/runtime/glnxa64:/home/limr/usr/MATLAB/v714/bin/glnxa64:/home/limr/usr/MATLAB/v714/sys/os/glnxa64:/home/limr/usr/MATLAB/v714/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:/home/limr/usr/MATLAB/v714/sys/java/jre/glnxa64/jre/lib/amd64/server:/home/limr/usr/MATLAB/v714/sys/java/jre/glnxa64/jre/lib/amd64
export XAPPLRESDIR = /home/limr/usr/MATLAB/v714/X11/app-defaults
MCR_DIR = $(HOME)/share/usr/MATLAB
GISTIC = $(HOME)/usr/gistic_2_0_21/gp_gistic2_from_seg
GISTIC_OPTS = -genegistic 0 -smallmem 1 -maxseg 5000 -savegene 1 -saveseg 1 -savedata 0 -v 30 -ta 0.4 -td 0.4 -js 15 -qvt 0.25 -conf 0.99 -broad 1 -brlen 0.5 -rx 0
DGV_FILE = $(HOME)/share/reference/GRCh37_hg19_variants_2013-07-23.txt

CNV_SIZES = 100000 300000

all : gistic_inputs $(foreach size,$(CNV_SIZES),gistic/gistic_cnv$(size).timestamp)
gistic_inputs : gistic/markersfile.txt gistic/segmentationfile.txt $(foreach size,$(CNV_SIZES),gistic/cnv.$(size).txt)

gistic/markersfile.txt :
	suppressPackageStartupMessages(library("rtracklayer"));
	targets <- import('$(TARGETS_FILE)')
	markers <- data.frame(chr = seqnames(targets), pos = start(targets))
	dir.create('$(@D)', showWarnings = F)
	write.table(markers, col.names = F, file = "$@", sep = "\t", quote = F, na = "")
	
gistic/segmentationfile.txt : PE := 8
gistic/segmentationfile.txt : MEM := 1G
gistic/segmentationfile.txt : $(foreach pair,$(SAMPLE_PAIRS),varscan/segment/$(pair).seg.txt)
	suppressPackageStartupMessages(library("rtracklayer"));
	suppressPackageStartupMessages(library("foreach"));
	suppressPackageStartupMessages(library("doMC"));
	segFiles <- unlist(strsplit("$^", " "));
	segNames <- sub(".*/", "", sub("\\..*", "", segFiles))
	targets <- import('$(TARGETS_FILE)')
	width(targets) <- 1
	registerDoMC(8)
	seg <- foreach (i = 1:length(segFiles), .combine = 'rbind') %dopar% {
		segFile <- segFiles[i]
		segName <- segNames[i]
		s <- read.delim(segFile, header = T, as.is = T, row.names = 1)
		s[['Chromosome']][s[['Chromosome']] == 23] <- "X"
		s[['Chromosome']][s[['Chromosome']] == 24] <- "Y"
		gr <- with(s, GRanges(seqnames = Chromosome, range = IRanges(start = Start, end = End), segmented = as.numeric(log2_ratio_seg)))
		redGr <- reduce(gr)
		x <- findOverlaps(redGr, gr, select = 'first')
		redGr$$segmented <- gr[x]$$segmented
		# reduced the genomic range, need to intersect with targets
		numMarkers <- countOverlaps(redGr, targets)
		Start <- start(targets)[findOverlaps(redGr, targets, select = 'first')]
		End <- start(targets)[findOverlaps(redGr, targets, select = 'last')]
		seg <- data.frame(segName, chrom = seqnames(redGr), start = Start, end = End, numMarkers, segmented = redGr$$segmented)
		seg <- subset(seg, numMarkers > 0)
		seg[!duplicated(seg), ]
	}
	splitSeg <- split(seg, list(as.factor(seg$$segName), as.factor(seg$$chrom)))
	seg <- do.call('rbind', lapply(splitSeg, function(x) {
		rx <- Rle(x$$segmented)
		nx <- x[start(rx), ]
		nx$$end <- x[end(rx), "end"]
		nx$$numMarkers <- aggregate(x$$numMarkers, rx, sum)
		nx
	}))
	dir.create('$(@D)', showWarnings = F)
	write.table(seg, file = "$@", sep = "\t", row.names = F, col.names = F, quote = F)

gistic/cnv.%.txt : gistic/markersfile.txt
	suppressPackageStartupMessages(library("GenomicRanges"));
	dgv <- read.delim("$(DGV_FILE)", as.is=T)
	dgv <- dgv[which(dgv[,5]=="CNV"), ]
	dgv <- dgv[, 1:4]
	dgv$$size = dgv[,4]-dgv[,3]+1
	dgv <- dgv[which(dgv$$size <= $*), ]
	dgv <- dgv[which(dgv$$chr %in% 1:22), ]
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


gistic/gistic_cnv%.timestamp : MEM := 8G
gistic/gistic_cnv%.timestamp : gistic/segmentationfile.txt gistic/markersfile.txt gistic/cnv.%.txt
	Sys.setenv(LD_LIBRARY_PATH = "/home/limr/usr/MATLAB/v714/runtime/glnxa64:/home/limr/usr/MATLAB/v714/bin/glnxa64:/home/limr/usr/MATLAB/v714/sys/os/glnxa64:/home/limr/usr/MATLAB/v714/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:/home/limr/usr/MATLAB/v714/sys/java/jre/glnxa64/jre/lib/amd64/server:/home/limr/usr/MATLAB/v714/sys/java/jre/glnxa64/jre/lib/amd64")
	Sys.setenv(XAPPLRESDIR = "/home/limr/usr/MATLAB/v714/X11/app-defaults")
	Sys.setenv(MCR_DIR = "$(HOME)/share/usr/MATLAB")
	dir.create('$(@D)/gistic_cnv$*', showWarnings = F, recursive = T)
	system("$(GISTIC) -b $(@D)/gistic_cnv$* -seg $< -mk $(<<) -refgene $(GISTIC_REF) -cnv $(<<<) $(GISTIC_OPTS) 2>&1 && touch $@")
