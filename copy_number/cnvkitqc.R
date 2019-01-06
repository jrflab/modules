#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

args_list <- list(make_option("--normal_files", default = NA, type = 'character', help = "normal input files"),
				  make_option("--tumor_files", default = NA, type = 'character', help = "tumor input files"),
				  make_option("--out_file", default = NA, type = 'character', help = "output file"),
				  make_option("--option", default = NA, type = 'character', help = "1-0 for ontarget or offtarget"))
parser <- OptionParser(usage = "%prog", option_list = args_list)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

normal_files = unlist(strsplit(x=opt$normal_files, split=" ", fixed=TRUE))
normal_samples = gsub(pattern=".cnr", replacement="", x=gsub(pattern="cnvkit/cnr/", replacement="", x=normal_files, fixed=TRUE), fixed=TRUE)
tumor_files = unlist(strsplit(x=opt$tumor_files, split=" ", fixed=TRUE))
tumor_samples = gsub(pattern=".cnr", replacement="", x=gsub(pattern="cnvkit/cnr/", replacement="", x=tumor_files, fixed=TRUE), fixed=TRUE)
out_file = opt$out_file

'MAD' <- function(x)
{
	x = na.omit(x)
	q2 = mad(x)
	return(invisible(q2))
}

'MAPD' <- function(x)
{
	x = na.omit(x)
	q2 = median(abs(x[1:(length(x)-1)] - x[2:length(x)]))
	return(invisible(q2))
}

'MIQR' <- function(x)
{
	x = na.omit(x)
	iq = stats::IQR(abs(x[1:(length(x)-1)] - x[2:length(x)]))
	return(invisible(iq))
}

'scale.' <- function(x)
{
	y = (x-min(x))/(max(x)-min(x))
	return(invisible(y))
}

'transparentRgb' <- function (col = "black", alpha = 85) 
{
    tmp = c(col2rgb(col), alpha, 255)
    names(tmp) = c("red", "green", "blue", "alpha", "maxColorValue")
    out = do.call("rgb", as.list(tmp))
    return(invisible(out))
}


qc = matrix(NA, nrow=length(c(normal_samples, tumor_samples)), ncol=3, dimnames=list(c(normal_samples, tumor_samples), c("MAD", "MAPD", "IQR")))
for (i in 1:length(normal_files)) {
	print(i)
	data = read.csv(file=normal_files[i], header=TRUE, sep="\t", stringsAsFactors=FALSE)
	index = data[,"chromosome"] %in% 1:22 & data[,"gene"] == ifelse(opt$option==1, "-", "Antitarget")
	qc[normal_samples[i],1] = MAD(data[index,"log2"])
	qc[normal_samples[i],2] = MAPD(data[index,"log2"])
	qc[normal_samples[i],3] = MIQR(data[index,"log2"])
}
for (i in 1:length(tumor_files)) {
	print(i)
	data = read.csv(file=tumor_files[i], header=TRUE, sep="\t", stringsAsFactors=FALSE)
	index = data[,"chromosome"] %in% 1:22 & data[,"gene"] == ifelse(opt$option==1, "-", "Antitarget")
	qc[tumor_samples[i],1] = MAD(data[index,"log2"])
	qc[tumor_samples[i],2] = MAPD(data[index,"log2"])
	qc[tumor_samples[i],3] = MIQR(data[index,"log2"])
}
data = qc
colnames(data) = c("MAD", "MAPD", "IQR")
data = cbind("SAMPLE_NAME"=c(normal_samples, tumor_samples), "SAMPLE_TYPE"=c(rep("N", length(normal_samples)), rep("T", length(tumor_samples))), data)
write.table(data, file=out_file, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

# MAPD
file_name = paste0("cnvkit/qc/", ifelse(opt$option==1, "on", "off"), "target_mapd.pdf")
x = as.numeric(data[data[,"SAMPLE_TYPE"]=="T", "MAPD"])
y = as.numeric(data[data[,"SAMPLE_TYPE"]=="N", "MAPD"])
dx = density(x, from=0, to=max(x,y))
dx$y = scale.(dx$y)
dy = density(y, from=0, to=max(x,y))
dy$y = scale.(dy$y)
pdf(file=file_name, width=7, height=7)
par(mar = c(6.1, 6.5, 4.1, 1.1))
plot(0, 0, type="n", axes = FALSE, frame.plot = FALSE, main = "", xlab = "", ylab = "", xlim=c(0, max(max(x, y), 1.5)), ylim=c(0,1.2))
polygon(x=c(dx$x, rev(dx$x)), y=c(dx$y, rep(0, length(dx$y))), border="steelblue", col=transparentRgb("steelblue", 155), lwd=2)
polygon(x=c(dy$x, rev(dy$x)), y=c(dy$y, rep(0, length(dy$y))), border="grey50", col=transparentRgb("grey50", 155), lwd=2)
legend("topright", col=c("steelblue", "grey50"), pch=15, legend=c("Tumor", "Normal"), box.lwd=-1)
axis(1, at = NULL, cex.axis = 1.5, padj = 0.25, lwd=1.25, lwd.ticks=1.15)
axis(2, at = seq(0, 1, by=.2), labels = seq(0, 1, by=.2), cex.axis = 1.5, las = 1, lwd=1.25, lwd.ticks=1.15)
mtext(side = 1, text = "MAPD", line = 4, cex = 1.5)
mtext(side = 2, text = "Density", line = 4, cex = 1.5)
dev.off()

# MAD
file_name = paste0("cnvkit/qc/", ifelse(opt$option==1, "on", "off"), "target_mad.pdf")
x = as.numeric(data[data[,"SAMPLE_TYPE"]=="T", "MAD"])
y = as.numeric(data[data[,"SAMPLE_TYPE"]=="N", "MAD"])
dx = density(x, from=0, to=max(x,y))
dx$y = scale.(dx$y)
dy = density(y, from=0, to=max(x,y))
dy$y = scale.(dy$y)
pdf(file=file_name, width=7, height=7)
par(mar = c(6.1, 6.5, 4.1, 1.1))
plot(0, 0, type="n", axes = FALSE, frame.plot = FALSE, main = "", xlab = "", ylab = "", xlim=c(0, max(max(x, y), 1.5)), ylim=c(0,1.2))
polygon(x=c(dx$x, rev(dx$x)), y=c(dx$y, rep(0, length(dx$y))), border="steelblue", col=transparentRgb("steelblue", 155), lwd=2)
polygon(x=c(dy$x, rev(dy$x)), y=c(dy$y, rep(0, length(dy$y))), border="grey50", col=transparentRgb("grey50", 155), lwd=2)
legend("topright", col=c("steelblue", "grey50"), pch=15, legend=c("Tumor", "Normal"), box.lwd=-1)
axis(1, at = NULL, cex.axis = 1.5, padj = 0.25, lwd=1.25, lwd.ticks=1.15)
axis(2, at = seq(0, 1, by=.2), labels = seq(0, 1, by=.2), cex.axis = 1.5, las = 1, lwd=1.25, lwd.ticks=1.15)
mtext(side = 1, text = "MAD", line = 4, cex = 1.5)
mtext(side = 2, text = "Density", line = 4, cex = 1.5)
dev.off()

# IQR
file_name = paste0("cnvkit/qc/", ifelse(opt$option==1, "on", "off"), "target_iqr.pdf")
x = as.numeric(data[data[,"SAMPLE_TYPE"]=="T", "IQR"])
y = as.numeric(data[data[,"SAMPLE_TYPE"]=="N", "IQR"])
dx = density(x, from=0, to=max(x,y))
dx$y = scale.(dx$y)
dy = density(y, from=0, to=max(x,y))
dy$y = scale.(dy$y)
pdf(file=file_name, width=7, height=7)
par(mar = c(6.1, 6.5, 4.1, 1.1))
plot(0, 0, type="n", axes = FALSE, frame.plot = FALSE, main = "", xlab = "", ylab = "", xlim=c(0, max(max(x, y), 1.5)), ylim=c(0,1.2))
polygon(x=c(dx$x, rev(dx$x)), y=c(dx$y, rep(0, length(dx$y))), border="steelblue", col=transparentRgb("steelblue", 155), lwd=2)
polygon(x=c(dy$x, rev(dy$x)), y=c(dy$y, rep(0, length(dy$y))), border="grey50", col=transparentRgb("grey50", 155), lwd=2)
legend("topright", col=c("steelblue", "grey50"), pch=15, legend=c("Tumor", "Normal"), box.lwd=-1)
axis(1, at = NULL, cex.axis = 1.5, padj = 0.25, lwd=1.25, lwd.ticks=1.15)
axis(2, at = seq(0, 1, by=.2), labels = seq(0, 1, by=.2), cex.axis = 1.5, las = 1, lwd=1.25, lwd.ticks=1.15)
mtext(side = 1, text = "IQR", line = 4, cex = 1.5)
mtext(side = 2, text = "Density", line = 4, cex = 1.5)
dev.off()

