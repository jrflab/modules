#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--file_in", default = NA, type = 'character', help = "input file name"),
				  make_option("--file_out", default = NA, type = 'character', help = "output file name"))

parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

chromStrToNum <- function(str) {
	suppressWarnings(cNum <- as.numeric(str))
	if (is.na(cNum) && str == "X" ) { 
		cNum <- 23
	} else if (is.na(cNum) && str == "Y") {
		cNum <- 24 
	}
	return(invisible(cNum))
}

GetChrominfo <- function() {
  f <- "modules/copy_number/hg19_chrominfo.txt"
  chrom <- read.table(file=f)
  chrom <- subset(chrom, grepl("^chr[0-9XY]{1,2}$", chrom[,1]))
  f <- "modules/copy_number/hg19_gaps.txt"
  gaps <- read.table(file=f)
  centro <- subset(gaps, gaps[,8] == "centromere")
  chrominfo <- merge(chrom[,1:2], centro[,2:4], by.x = 1, by.y = 1) 
  chrominfo$centromere <- rowMeans(chrominfo[,3:4]) 
  chrominfo <- chrominfo[,c(1,2,5,3,4)] 
  colnames(chrominfo) <- c("chr", "size", "centromere", "centstart", "centend")
  chrominfo[,1] <- as.character(chrominfo[,1])
  chrominfo$chr <- sub("chr", "", chrominfo$chr)
  chrominfo$chr <- sub("X", "23", chrominfo$chr)
  chrominfo$chr <- sub("Y", "24", chrominfo$chr)
  chrominfo[,1] <- as.numeric(chrominfo[,1])
  chrominfo <- chrominfo[order(chrominfo$chr), ]  
  rownames(chrominfo) <- as.character(chrominfo[,1])
  chrominfo <- as.matrix(chrominfo)
  return(invisible(chrominfo))
}

fix_facets_column_names <- function(dat) {
	colnames(dat)[which(colnames(dat)=="chrom")] <- "chromosome"
	colnames(dat)[which(colnames(dat)=="loc.start")] <- "startBP"
	colnames(dat)[which(colnames(dat)=="loc.end")] <- "endBP"
	colnames(dat)[which(colnames(dat)=="lcn.em")] <- "nB"
	sz <- dat[,"endBP"] - dat[,"startBP"]
	dat <- cbind(dat, size=sz)
    nA <- dat[,"tcn.em"] - dat[,"nB"]
    dat <- cbind(dat, nA=nA)
	return(invisible(dat))
}

join_adjacent_segments <- function(dat) {
	cur_segs <- dat
	something_changed <- 1
	while ( something_changed ) {
		new_segs <- c()
		something_changed <- 0
		x <- 2
		last_changed <- 0
		while (x <= nrow(cur_segs)) {
			last_changed <- 0
			if ( 	(cur_segs[x-1,"nB"] == cur_segs[x,"nB"]) && 
					(cur_segs[x-1,"nA"] == cur_segs[x,"nA"]) &&
					(cur_segs[x-1,"chromosome"] == cur_segs[x,"chromosome"])
			) {
				t <- cur_segs[x-1,]
				t["endBP"] <- cur_segs[x,"endBP"]
				t["end"] <- cur_segs[x,"end"]
				t["size"] <- t["endBP"] - t["startBP"]
				something_changed <- 1
				new_segs <- rbind(t, new_segs)
				x <- x+2
				last_changed <- 1
			} else {
				new_segs <- rbind(cur_segs[x-1,], new_segs)
				x<-x+1
			}
		}
		if (! last_changed ) {
			new_segs <- rbind(cur_segs[x-1,],new_segs)
		}
		n <- nrow(new_segs)
		new_segs <- new_segs[n:1,]
		cur_segs <- new_segs
	}	
	return(invisible(cur_segs))
}

fix_facet_segs <- function(dat) {
    i <- which(is.na(dat$nB))
    if ( length(i) > 0 )  {
        dat <- dat[-i, ]
    }
    dat <- join_adjacent_segments(dat)
    return(invisible(dat))
}

chrom_arm_LST_score <- function(dat) {
	score <- 0
	segs <- c()
	SIZE_THRESH <- 10e6
	SPACE_THRESH <- 3e6
	if ( nrow(dat) >= 2 ) {
		for (x in 2:nrow(dat)) {
			if ( 	(dat[x-1,"size"] >= SIZE_THRESH) && 
					(dat[x,"size"] >= SIZE_THRESH) &&
					( (dat[x,"startBP"] - dat[x-1,"endBP"]) <= SPACE_THRESH)
			) {
				score <- score +1
				segs <- rbind(dat[x-1,], segs)
			}
		}
	}
	tmp <- list()
	tmp$score <- score
	tmp$segs <- segs
	return(invisible(tmp))
}

lst_filter <- function(dat, size_thresh) {
	i <- which(dat[,"size"] < size_thresh)
	sz <- dat[i,"size"]
	i <- i[order(sz)]
	segs_removed <- 0
	while (length(i) > 0) {
		dat <- dat[-i[1], ]
		dat <- join_adjacent_segments(dat)
		i<- which(dat[,"size"] < size_thresh)
		sz <- dat[i,"size"]
		i <- i[order(sz)]	
		segs_removed <- segs_removed + 1
	}
	return(invisible(dat))
}

score_LST <- function(dat, chromInfo) {
	score <- 0
	segs <- c()
	dat <- lst_filter(dat, 3e6)
	for (c in unique(dat[,"chromosome"]) ) {
		i <- which(dat[,"chromosome"] == c)
		csegs <- dat[i,]
		cNum <- chromStrToNum(c)
		i <- which(csegs[,"startBP"] <= chromInfo[cNum,"centstart"])
		parm <- csegs[i,]
		tmp <- chrom_arm_LST_score(parm)
		score <- score + tmp$score
		segs <- rbind(tmp$segs, segs)
		i <- which(csegs[,"endBP"] >= chromInfo[cNum,"centend"])
		qarm <- csegs[i,]
		tmp <- chrom_arm_LST_score(qarm)
		score <- score + tmp$score
		segs <- rbind(tmp$segs, segs)
	}
	tmp <- list()
	tmp$score <- score
	tmp$segs <- segs
	return(invisible(tmp))
}

dat = read.table(opt$file_in, sep="\t", header=TRUE, stringsAsFactor=FALSE)
dat = fix_facets_column_names(dat)
segs = fix_facet_segs(dat)
chromInfo = GetChrominfo()
lst = score_LST(segs, chromInfo)
cat(paste0(gsub("facets/cncf/","", gsub(".cncf.txt", "", opt$file_in)), "\t", lst$score), file = opt$file_out, append=FALSE)
cat("\n", file = opt$file_out, append=TRUE)

warnings()

