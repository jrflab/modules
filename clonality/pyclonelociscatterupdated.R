#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("mclust"))
suppressPackageStartupMessages(library("MASS"))

optList = list(make_option("--sample_name", default = NULL, help = "sample name"))

parser = OptionParser(usage = "%prog [options] mutation_file", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

'post_density' <- function (x)
{
	y = density(x=x, adjust=2)
	return(invisible(y))
}

'hex_cols' <- function(x)
{
	x = x%%8
	if (x==0) {
		x = 8
	}
	cols = c("#4865B1", "#FFA500", "#B22034", "#E9E0BA", "#D5D5D5", "#000000", "#DC0073", "#00A1E4")
	return(cols[x])
}

data = read.csv(file=paste0("pyclone/", opt$sample_name, "/summary.tsv"), header=TRUE, sep="\t", stringsAsFactors=FALSE)
sample_names = gsub("_ucf", "", x=colnames(data)[grep("_ucf", colnames(data), fixed=TRUE)], fixed=TRUE)
DP = data[,paste0("DP_", sample_names),drop=FALSE]
AD = round(data[,paste0("MAF_", sample_names),drop=FALSE] * DP)
MAF = data[,paste0("MAF_", sample_names),drop=FALSE]
CF = data[,paste0(sample_names, "_ucf"),drop=FALSE]

pdf(file=paste0("pyclone/", opt$sample_name, "/plots/all_loci_scatter_filtered.pdf"))
par(mar=c(6.1, 6.5, 4.1, 1.1))
for (i in 1:(length(sample_names)-1)) {
	for (j in (i+1):length(sample_names)) {
		plot(0, 0, type="n", axes=FALSE, frame.plot=FALSE, main="", xlab="", ylab="", xlim=c(0,1), ylim=c(0,1))
		dp_x = DP[,i]
		dp_y = DP[,j]
		ad_x = AD[,i]
		ad_y = AD[,j]
		maf_x = MAF[,i]
		maf_y = MAF[,j]
		cf_x = CF[,i]
		cf_y = CF[,j]
		
		ind = is.na(dp_x) | is.na(dp_y) | is.na(ad_x) | is.na(ad_y) | is.na(maf_x) | is.na(maf_y) | is.na(cf_x) | is.na(cf_y)
		dp_x = dp_x[!ind]
		dp_y = dp_y[!ind]
		ad_x = ad_x[!ind]
		ad_y = ad_y[!ind]
		maf_x = maf_x[!ind]
		maf_y = maf_y[!ind]
		cf_x = cf_x[!ind]
		cf_y = cf_y[!ind]
		
		ind = dp_x>=50 & dp_y>=50
		dp_x = dp_x[ind]
		dp_y = dp_y[ind]
		ad_x = ad_x[ind]
		ad_y = ad_y[ind]
		maf_x = maf_x[ind]
		maf_y = maf_y[ind]
		cf_x = cf_x[ind]
		cf_y = cf_y[ind]
		
		ind = ad_x<5 | maf_x<.03
		cf_x[ind] = 0
		ind = ad_y<5 | maf_y<.03
		cf_y[ind] = 0
		
		x = jitter(cf_x, amount=.01)
		x[x>1] = 1
		x[x<0] = 0
		y = jitter(cf_y, amount=.01)
		y[y>1] = 1
		y[y<0] = 0
		
		contour(kde2d(x, y, n=100, lims = c(c(-.01,1.01),c(-.01,1.01))), drawlabels=FALSE, nlevels=15, add=TRUE, lwd=.5, col=hex_cols(3))
		points(x, y, type="p", pch=1, col=hex_cols(3))
	    axis(1, at=seq(from=0, to=1, by=.2), labels=seq(from=0, to=1, by=.2), cex.axis=1.5, padj=0.25, lwd = 1.25, lwd.ticks = 1.15)
	    axis(2, at=seq(from=0, to=1, by=.2), labels=seq(from=0, to=1, by=.2), cex.axis=1.5, las=1, lwd = 1.25, lwd.ticks = 1.15)
	    points(c(.1,.1), c(-.1,1), type="l", col="orange", lty=3)
	    points(c(.9,.9), c(-.1,1), type="l", col="orange", lty=3)
	    points(c(-.1,1), c(.1,.1), type="l", col="orange", lty=3)
	    points(c(-.1,1), c(.9,.9), type="l", col="orange", lty=3)
	    mtext(side=1, text=sample_names[i], line=4, cex=1.5)
	    mtext(side=2, text=sample_names[j], line=4, cex=1.5)
	}	    
}
dev.off()


for (i in 1:length(sample_names)) {

	dp_x = data[,paste0("DP_", sample_names[i]),]
	ad_x = round(data[,paste0("MAF_", sample_names[i]),] * dp_x)
	maf_x = data[,paste0("MAF_", sample_names[i]),]
	cf_x = data[,paste0(sample_names[i], "_ucf"),]
		
	ind = is.na(dp_x) | is.na(ad_x) | is.na(maf_x) | is.na(cf_x)
	dp_x = dp_x[!ind]
	ad_x = ad_x[!ind]
	maf_x = maf_x[!ind]
	cf_x = cf_x[!ind]
	data = data[!ind,,drop=FALSE]
		
	ind = dp_x>=50
	dp_x = dp_x[ind]
	ad_x = ad_x[ind]
	maf_x = maf_x[ind]
	cf_x = cf_x[ind]
	data = data[ind,,drop=FALSE]
			
	ind = ad_x<5 | maf_x<.03
	maf_x[ind] = 0
	cf_x[ind] = 0
	
	data[,paste0("MAF_", sample_names[i])] = maf_x
	data[,paste0(sample_names[i], "_ucf")] = cf_x

}

MAF = data[,paste0("MAF_", sample_names),drop=FALSE]
CF = data[,paste0(sample_names, "_ucf"),drop=FALSE]

index = apply(MAF, 1, function(x) {sum(x==0)==length(x)}) | apply(CF, 1, function(x) {sum(x==0)==length(x)})
data = data[!index,,drop=FALSE]

write.table(data, file=paste0("pyclone/", opt$sample_name, "/summary_filtered.tsv"), sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
