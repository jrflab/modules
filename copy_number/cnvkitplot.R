#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("copynumber"))
suppressPackageStartupMessages(library("colorspace"))
suppressPackageStartupMessages(library("ASCAT"))
suppressPackageStartupMessages(library("GAP"))

'plot_log2_' <- function(x, title = "")
{
   	par(mar=c(5, 5, 4, 2)+.1)
   	data("CytoBand")
   	end = NULL
   	for (i in 1:23) {
   		end = c(end, max(CytoBand[CytoBand[,1]==i,"End"]))
   	}
   	end = cumsum(end)
   	start = c(1, end[1:22]+1)
   	CytoBand = cbind(start, end)
   	index = NULL
   	for (i in 1:23) {
   		index = c(index, seq(from = CytoBand[i, "start"], to=CytoBand[i, "end"], length=sum(x$chromosome==i)))
   	}
	plot(index, x$log2, type="p", pch=".", cex=1.95, col="grey80", axes=FALSE, frame=TRUE, xlab="", ylab="", main="", ylim=c(-4,5))
  	axis(2, at = c(-4, -2, 0, 2, 4), labels = c(-4, -2, 0, 2, 4), cex.axis = 1, las = 1)
	mtext(side = 2, text = expression(Log[2]~"Ratio"), line = 3.15, cex = 1.25)
	abline(v=1, col="goldenrod3", lty=3, lwd=.5)
	abline(h=0, col="red", lty=1, lwd=1)
	for (j in 1:23) {
		abline(v=CytoBand[j,"end"], col="goldenrod3", lty=3, lwd=.5)
	}
	axis(1, at = .5*(CytoBand[,"start"]+CytoBand[,"end"]), labels=c(1:22, "X"), cex.axis = 0.85, las = 1)
	rect(xleft=1-1e10, xright=CytoBand[23,"end"]+1e10, ybottom=4, ytop=6, col="lightgrey", border="black", lwd=1.5)
	title(main = title, line=-1, cex.main=.75, font.main=1)
    box(lwd=1.5)
}

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--in_file", default = NA, type = 'character', help = "input file name"))
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

outfile_on_target = gsub("cnr", "log2", gsub(".cnr", ".ontarget.pdf", opt$in_file, fixed=TRUE), fixed=TRUE)
outfile_off_target = gsub("cnr", "log2", gsub(".cnr", ".offtarget.pdf", opt$in_file, fixed=TRUE), fixed=TRUE)

data = read.table(file=opt$in_file, header=TRUE, sep="\t", comment.char="#", stringsAsFactors=FALSE)
data = subset(data, data[,"depth"]!=0)

if (nrow(data)==0) {
	system(paste0("touch ", outfile_on_target))
	system(paste0("touch ", outfile_off_target))
} else {
	data[,"chromosome"] = gsub(pattern="chr", replacement="", x=data[,"chromosome"], fixed=TRUE)
	data[data[,"chromosome"]=="X", "chromosome"] = 23
	data[data[,"chromosome"]=="Y", "chromosome"] = 24
	data[,"chromosome"] = as.numeric(data[,"chromosome"])
	data = subset(data, data[,"chromosome"]<=23)
	
	if (sum(data$gene=="-")>0) {
		flag = 1
	} else if (sum(data$gene=="Antitarget")>0) {
		flag = 2
	}
	
	if (flag==1) {
		ontarget = subset(data, data$gene=="-")
	} else if (flag==2) {
		ontarget = subset(data, data$gene!="Antitarget")
	}
	
	pdf(file=outfile_on_target, width=10, height=4.25)
	plot_log2_(x=ontarget, title=gsub("cnvkit/cnr/", "", gsub(".cnr", "", opt$in_file, fixed=TRUE), fixed=TRUE))
	dev.off()
	
	if (flag==1) {
		offtarget = subset(data, data$gene!="-")
	} else if (flag==2) {
		offtarget = subset(data, data$gene=="Antitarget")
	}
	
	tmp = offtarget[,c("chromosome", "start", "log2"),drop=FALSE]
	tmp = winsorize(data=tmp, tau=3.5, k=25, verbose=FALSE, return.outliers=TRUE)
	offtarget[tmp$wins.outliers[,3]!=0,"log2"] = NA
	pdf(file=outfile_off_target, width=10, height=4.25)
	plot_log2_(x=offtarget, title=gsub("cnvkit/cnr/", "", gsub(".cnr", "", opt$in_file, fixed=TRUE), fixed=TRUE))
	dev.off()
}
