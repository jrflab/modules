################## 
## an attempt to replicate the copy number procedures for MSK-IMPACT with modifications
#################

# the MSK-IMPACT paper specifies using GATK DoC but didn't say exactly how
# so split the BED file into approximate 100-bp bins and then DoC
# but leave a gap between intervals otherwise DoC merges them

bed <- read.delim("IMPACT410_TARGETS.bed", as.is=T, header=F)
bed$size <- bed$V3-bed$V2
bed$bins <- floor(bed$size/100)

newbed <- apply(bed, 1, function(x) {
	if (as.numeric(x[8])<=1 ) { x[1:6]
	} else {
		avgsize <- as.numeric(x[7])/as.numeric(x[8]); #print(avgsize);
		numlarger <- as.numeric(x[7])-(as.numeric(x[8])*floor(avgsize))
		ends <- as.numeric(cumsum(unlist(c(rep(ceiling(avgsize), numlarger), rep(floor(avgsize), as.numeric(x[8])-numlarger)))))
		starts <- c(0, ends[-length(ends)])
		starts <- as.numeric(x[2])+starts
		ends <- as.numeric(x[2])+ends

		cbind(x[1], starts+2, ends, x[4], x[5], x[6])
	}
})
newbed <- do.call("rbind", newbed)
newbed[,2] <- as.numeric(newbed[,2])
newbed[,3] <- as.numeric(newbed[,3])
newbed[,4] <- 1:nrow(newbed)

write.table(newbed[,1:6], file="IMPACT410.modified.bed", sep="\t", row.names=F, quote=F, col.names=F)

### DoC ###
#qsub ~/scripts2/qCMD.med.sh java -Xmx8g -jar /home/ngk1/share/usr/lib/java/GenomeAnalysisTK-3.3-0.jar -T DepthOfCoverage -R /home/ngk1/share/reference/GATK_bundle/2.3/human_g1k_v37.fa -I input_bams.list -o IMPACT410.modified.bed -L IMPACT410.modified.bed --omitDepthOutputAtEachBase --omitPerSampleStats --omitLocusTable;

########## normalisation of DoC output based on Loess on GC content

library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(limma)
library(DNAcopy)
library(snow)

setwd("~/Documents/Miscellaneous/20140905_sync_endometrioid/impact410/doc/")

files <- dir(pattern="interval_summary")
doc <- lapply(files, function(x) { tab <- read.delim(x, as.is=T, row.names=1); })
doc <- do.call("rbind", doc)

doc <- doc[,grep("total_cvg", colnames(doc))]
doc <- doc[,which(colnames(doc)!="FFPEPOOLEDNORMAL")]

# step 1: square-root transformed
doc <- sqrt(doc)

# step 2: Loess normalization of depth based on GC content
gcContent <- function(x) {
	alf <- alphabetFrequency(x, as.prob=TRUE)
	sum(alf[c("G", "C")])}

pos <- lapply(rownames(doc), function(x) { p <- strsplit(x, ":", fixed=T)[[1]]; p <- c(p[1], strsplit(p[2], "-", fixed=T)[[1]])})

cl <- makeCluster(4)
gc <- unlist(parLapply(cl, pos, function(x, gcContent) { 
	library(BSgenome)
	library(BSgenome.Hsapiens.UCSC.hg19)
	gcContent(Hsapiens[[paste("chr", x[1], sep="")]][x[2]:x[3]])}, gcContent))
stopCluster(cl)

doc_norm <- apply(doc, 2, function(x) { x+loessFit(x, gc)$residual })
rownames(doc_norm) <- rownames(doc)
colnames(doc_norm) <- gsub("_total_cvg", "", colnames(doc_norm))

# step 3: filter out regions with normalised depth in the top or bottom 5% in >=20% of normal samples
ft_low <- apply(doc_norm[,grep("N", colnames(doc_norm))], 2, function(x) { x<quantile(x, 0.05) })
toft_low <- which(apply(ft_low, 1, function(x) { length(which(x))>=(length(x)*0.2)}))

ft_high <- apply(doc_norm[,grep("N", colnames(doc_norm))], 2, function(x) { x>quantile(x, 0.95) })
toft_high <- which(apply(ft_high, 1, function(x) { length(which(x))>=(length(x)*0.2)}))

toft <- unique(c(toft_low, toft_high))
pos_ft <- matrix(unlist(pos), ncol=3, byrow=3)
pos_ft[which(pos_ft[,1]=="X"),1] <- 23
pos_ft[which(pos_ft[,1]=="Y"),1] <- 24
for (i in 1:ncol(pos_ft)) { pos_ft[,i] <- as.numeric(pos_ft[,i])}

if (length(toft)>0) { doc_norm <- doc_norm[-toft, , drop=F]; pos_ft <- pos_ft[-toft,] }

# step 4: compute log ratios
sample_sheet <- read.delim("sample_sheet.txt", header=F, as.is=T, sep=" ")
genomedat <- apply(sample_sheet, 1, function(x, cov_norm) { scale(log(cov_norm[,x[1]]/cov_norm[,x[2]], base=2), scale=F)}, doc_norm)
colnames(genomedat) <- sample_sheet[,1]

######################### THE REST SHOULD BE THE SAME (ALMOST THE SAME) AS VARSCAN COPYNUMBER #############
# step 5: segmentation
cna <- CNA(genomedat, as.numeric(pos_ft[,1]), round(as.numeric(pos_ft[,2])+as.numeric(pos_ft[,3])/2), sampleid= sample_sheet[,1])
smoothed.cna <- smooth.CNA(cna, outlier.SD.scale=3, trim=0.01)
seg <- segment(smoothed.cna, undo.SD=1, alpha=0.00000001, undo.splits="sdundo")


# step 6: plot (copied from existing script)

for (i in unique(seg$output$ID)) {

png(paste(i,".seg_plot.png", sep=""), type = 'cairo-png', height=400, width=2000)
obj <- subset(seg, sample=i)

objdat <- obj$data[which(!is.na(obj$data[,3])),]

plot(objdat[,3], pch=20, xlab='Position', ylab="Copy number", xaxt='n', ylim=c(min(objdat[,3]), max(objdat[,3])+0.5))
points(unlist(apply(obj$output, 1, function(x) {rep(x[6], x[5])})), pch = 20, col = 'blue')
abline(v=cumsum(rle(as.vector(objdat$chrom))$lengths), col="red", lty=3)

    cen <- read.table("centromere2.txt", sep = '\t')
    for (j in unique(cen[,1])) {
        pos <- cen[which(cen[,1]==j)[1],3]
        index <- which(objdat$chrom==j & objdat$maploc > pos)[1]
        if (!is.na(index)) {
            abline(v=index, col="darkgrey", lty=3)
        }
        text(cumsum(rle(as.vector(objdat$chrom))$lengths)-((rle(as.vector(objdat$chrom))$lengths)/2), max(objdat[,3])+0.5-0.25)
    }

dev.off()


}