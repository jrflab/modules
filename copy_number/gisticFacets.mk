include modules/Makefile.inc

LOGDIR = log/gistic.$(NOW)

SHELL = modules/scripts/Rshell
.SHELLFLAGS = -m $(MEM) -s -p $(PE) -n $(@F) -l $(LOGDIR) -e 

.ONESHELL:
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: all gistic_inputs gistic_heatmaps

MEM := 2G
PE := 1

export LD_LIBRARY_PATH = /home/limr/usr/MATLAB/v714/runtime/glnxa64:/home/limr/usr/MATLAB/v714/bin/glnxa64:/home/limr/usr/MATLAB/v714/sys/os/glnxa64:/home/limr/usr/MATLAB/v714/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:/home/limr/usr/MATLAB/v714/sys/java/jre/glnxa64/jre/lib/amd64/server:/home/limr/usr/MATLAB/v714/sys/java/jre/glnxa64/jre/lib/amd64
export XAPPLRESDIR = /home/limr/usr/MATLAB/v714/X11/app-defaults
MCR_DIR = $(HOME)/share/usr/MATLAB
GISTIC = $(HOME)/share/usr/gistic_2_0_21/gp_gistic2_from_seg
GISTIC_THRESHOLD ?= 0.3
GISTIC_JS ?= 15
GISTIC_OPTS = -genegistic 0 -smallmem 1 -maxseg 5000 -savegene 1 -saveseg 1 -savedata 0 -v 30 -ta $(GISTIC_THRESHOLD) -td $(GISTIC_THRESHOLD) -js $(GISTIC_JS) -qvt 0.25 -conf 0.99 -broad 1 -brlen 0.5 -rx 0
DGV_FILE = $(HOME)/share/reference/GRCh37_hg19_variants_2013-07-23.txt
PLOT_GISTIC_HEATMAP = $(RSCRIPT) modules/copy_number/plotGisticHeatmap.R

CNV_SIZES = 100000 300000

all : gistic_inputs $(foreach size,$(CNV_SIZES),gistic/gistic_cnv$(size).timestamp) gistic_heatmaps
gistic_inputs : gistic/markersfile.txt gistic/segmentationfile.txt $(foreach size,$(CNV_SIZES),gistic/cnv.$(size).txt)
gistic_heatmaps : $(foreach size,$(CNV_SIZES),gistic/gistic_cnv$(size)/gistic_cnv_heatmap.pdf)

gistic/markersfile.txt : gistic/segmentationfile.txt
	suppressPackageStartupMessages(library("rtracklayer"));
	seg <- read.table('$<', sep = '\t', stringsAsFactors = F, col.names = c('samplePair', 'chr', 'start', 'end', 'numMarkers', 'logRatio'))
	targets <- import('$(TARGETS_FILE)')
	markers <- data.frame(chr = seqnames(targets), pos = start(targets))
	markers <- markers[markers$$chr %in% seg$$chr, ]
	dir.create('$(@D)', showWarnings = F)
	write.table(markers, col.names = F, file = "$@", sep = "\t", quote = F, na = "")
	
gistic/segmentationfile.txt : PE := 8
gistic/segmentationfile.txt : MEM := 1G
gistic/segmentationfile.txt : $(foreach pair,$(SAMPLE_PAIRS),facets/$(pair).cncf.txt)
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
		s <- read.delim(segFile, header = T, as.is = T, row.names = 1, col.names=c("Chromosome","Start","End","adjusted_log_ratio","nhet","log2_ratio_seg","mafR","segclust","cnlr.median.clust","mafR.clust","cf","tcn","lcn","cf.em","tcn.em","lcn.em"))
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


gistic/gistic_cnv%.timestamp : MEM := 12G
gistic/gistic_cnv%.timestamp : gistic/segmentationfile.txt gistic/markersfile.txt gistic/cnv.%.txt
	Sys.setenv(LD_LIBRARY_PATH = "/home/limr/usr/MATLAB/v714/runtime/glnxa64:/home/limr/usr/MATLAB/v714/bin/glnxa64:/home/limr/usr/MATLAB/v714/sys/os/glnxa64:/home/limr/usr/MATLAB/v714/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:/home/limr/usr/MATLAB/v714/sys/java/jre/glnxa64/jre/lib/amd64/server:/home/limr/usr/MATLAB/v714/sys/java/jre/glnxa64/jre/lib/amd64")
	Sys.setenv(XAPPLRESDIR = "/home/limr/usr/MATLAB/v714/X11/app-defaults")
	Sys.setenv(MCR_DIR = "$(HOME)/share/usr/MATLAB")
	dir.create('$(@D)/gistic_cnv$*', showWarnings = F, recursive = T)
	system("umask 002; $(GISTIC) -b $(@D)/gistic_cnv$* -seg $< -mk $(<<) -refgene $(GISTIC_REF) -cnv $(<<<) $(GISTIC_OPTS) 2>&1 && touch $@")

gistic/gistic_cnv%/gistic_cnv_heatmap.pdf : gistic/gistic_cnv%.timestamp
	system("$(PLOT_GISTIC_HEATMAP) --out $@ gistic/gistic_cnv$*/all_thresholded.by_genes.txt")



