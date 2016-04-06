#!/usr/bin/env Rscript
#### turn segmented copy number data to gene-based copy number with findOverlaps
## define HomDel as TCN=0, loss as TCN<ploidy, gain as TCN>ploidy, amp as TCN>=ploidy+4
## where ploidy= mode of TCN
### some variant of the below, also need one for the breast panel, IMPACT310 and exome

#---------------
# initialization
#---------------

# load base libraries
suppressMessages(pacman::p_load(optparse,RColorBrewer,GenomicRanges,plyr,dplyr,readr,stringr,tidyr,purrr,magrittr,rlist,crayon,foreach,Cairo,RMySQL,rtracklayer,colorspace,ggplot2,grid,gridExtra,RColorBrewer))
suppressPackageStartupMessages(library("facets",lib.loc="/home/bermans/R-dev/"));

#--------------
# parse options
#--------------

optList <- list(
				make_option("--outFile", default = NULL, help = "output file"),
				make_option("--mysqlHost", default = '10.0.200.48', help = "MySQL server hostname"),
				make_option("--mysqlPort", default = 38493, help = "MySQL server port"),
				make_option("--mysqlUser", default = 'embl', help = "MySQL server username"),
				make_option("--mysqlPassword", default = 'embl', help = "MySQL server password"),
				make_option("--mysqlDb", default = 'homo_sapiens_core_75_37', help = "MySQL server database"),
				make_option("--genesFile", default = NULL, help = "list of genes to include (hgnc symbols)"))
parser <- OptionParser(usage = "%prog [options] [facets files]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
	cat("Need facets output files\n")
	print_help(parser);
	stop();
} else if (is.null(opt$outFile)) {
	cat("Need output prefix\n")
	print_help(parser);
	stop();
} else {
	facetsFiles <- arguments$args
}

connect <- function() dbConnect(MySQL(), host = opt$mysqlHost, port = opt$mysqlPort, user = opt$mysqlUser, password = opt$mysqlPassword, dbname = opt$mysqlDb)
cat('Connecting to ensembl ... ')
mydb <- connect()
on.exit(dbDisconnect(mydb))

query <- "select r.name as chrom,
g.seq_region_start as start,
g.seq_region_end as end,
x.display_label as hgnc,
k.band as band
from gene as g
join seq_region as r on g.seq_region_id = r.seq_region_id
join xref as x on g.display_xref_id = x.xref_id
left join karyotype k on g.seq_region_id = k.seq_region_id
and ((g.seq_region_start >= k.seq_region_start and g.seq_region_start <= k.seq_region_end)
or (g.seq_region_end >= k.seq_region_start and g.seq_region_end <= k.seq_region_end))
where x.external_db_id = 1100;"
repeat {
	rs <- try(dbSendQuery(mydb, query), silent = T)
	if (is(rs, "try-error")) {
		cat("Lost connection to mysql db ... ")
		mydb <- connect()
		cat("reconnected\n")
	} else {
		break
	}
}
genes <- dbFetch(rs, -1)
cat(paste("Found", nrow(genes), "records\n"))

genes %<>% filter(chrom %in% as.character(c(1:22, "X", "Y"))) %>%
	filter(!duplicated(hgnc)) %>% 
	arrange(as.integer(chrom), start, end)

if (!is.null(opt$genesFile)) {
	g <- scan(opt$genesFile, what = 'character')
	genes %<>% filter(hgnc %in% g)
	absentGenes <- g[!g %in% genes$hgnc]
	if (length(absentGenes) > 0) {
		print("Unable to find", length(absentGenes), "in database\n");
		cat(absentGenes, sep = '\n');
	}
}

cat(paste("Filtering to", nrow(genes), "records\n"))

genesGR <- genes %$% GRanges(seqnames = chrom, ranges = IRanges(start, end), band = band, hgnc = hgnc)
			
mm <- lapply(facetsFiles, function(f) {
	tab <- read.delim(f, as.is=T)
	tab$chrom[which(tab$chrom==23)] <- "X"

	tabGR <- tab %$% GRanges(seqnames = chrom, ranges = IRanges(loc.start, loc.end))
	mcols(tabGR) <- tab %>% select(cnlr.median:lcn.em)

	fo <- findOverlaps(tabGR, genesGR)

	df <- as.data.frame(cbind(mcols(genesGR)[subjectHits(fo),], mcols(tabGR)[queryHits(fo),]))
	df %<>% group_by(hgnc) %>% top_n(1, abs(cnlr.median))

	ploidy <- table(df$tcn)
	ploidy <- as.numeric(names(ploidy)[which.max(ploidy)])

	df$GL <- 0
	df$GL[df$tcn < ploidy] <- -1
	df$GL[df$tcn == 0] <- -2
	df$GL[df$tcn > ploidy] <- 1
	df$GL[df$tcn >= ploidy + 4] <- 2

	load(gsub("cncf.txt", "Rdata", f, fixed=T))
	noise <- median(abs(out2$jointseg$cnlr-  unlist(apply(out2$out[,c("cnlr.median", "num.mark")], 1, function(x) {rep(x[1], each=x[2])}))))

	lrr <- sort(out2$jointseg$cnlr)
	if (noise <= 0.2) { lrr <- lrr[round(0.25*length(lrr)):round(0.75*length(lrr))]
	} else if ( noise <= 0.3 ) { lrr <- lrr[round(0.275*length(lrr)):round(0.725*length(lrr))]
	} else { lrr <- lrr[round(0.3*length(lrr)):round(0.7*length(lrr))]}

	df$GL2 <- 0
	df$GL2[df$cnlr.median < median(lrr)-(2.5*sd(lrr))] <- -1
	df$GL2[df$cnlr.median < median(lrr)-(7*sd(lrr))] <- -2
	df$GL2[df$cnlr.median > median(lrr)+(2*sd(lrr))] <- 1
	df$GL2[df$cnlr.median > median(lrr)+(6*sd(lrr))] <- 2

	df %>% select(hgnc, GL, GL2) %>% ungroup
})
names(mm) <- facetsFiles
for (f in facetsFiles) {
	n <- sub('\\..*', '', sub('.*/', '', f))
	colnames(mm[[f]])[2:3] <- paste(n, c("EM", "LRR_threshold"), sep="_")
}

mm <- left_join(genes, join_all(mm, type = 'full', by="hgnc")) %>% arrange(as.integer(chrom), start, end)
write.table(mm, file=opt$outFile, sep="\t", row.names=F, na="", quote=F)

#----------------
# fix facets file
#----------------

cncf.path <- 'facets/'
plot.file <- 'facets/geneCN.pdf'

cnv.matrix <-
	mm %>%
	mutate(chrom=
		ifelse(chrom=='X',23,
		ifelse(chrom=='Y',23,chrom))) %>%
	mutate(chrom=as.integer(chrom))

chrom <- cnv.matrix$chrom

cnv.matrix %<>% arrange(chrom,start,end)

cnv.matrix[is.na(cnv.matrix)] <- 3  # replace NA with numeric for use in rle function
breaks <- c(0,cumsum(table(chrom)[chrom %>% table %>% names %>% as.numeric %>% order])) # start-1 == end

for (sample.name in names(cnv.matrix) %>% list.filter(!. %in% c('chrom','start','end','hgnc','band'))){

	cat('\n * sample:',sample.name)

	#loop over chromosomes in junction table
	for(chromosome in 1:length(cnv.matrix$chrom %>% unique)){

		#cat('\n   chr',formatC(chromosome, width=2, flag='0'),': ',sep='')
		start <- breaks[chromosome]+1
		end <- breaks[chromosome+1]
		chbit <- cnv.matrix[start:end,sample.name] %>% unlist %>% rle

		# replace an entirely empty chromosome with calls of 0
		if((chbit$lengths %>% length)==1){
			cnv.matrix[start:end,sample.name] <- 0
		}else{
			# extend chromosome start calls from nearest integer on same chromosome
			if(cnv.matrix[start,sample.name]==3){
				#cat('start..')
				NAbottom <- min(which(chbit$values==3))
				Ibottom <- start + cumsum(chbit$lengths)[NAbottom]
				cnv.matrix[start:(Ibottom-1),sample.name] <- cnv.matrix[Ibottom,sample.name]
			}
			# extend chromosome end calls from nearest integer on same chromosome
			if(cnv.matrix[end,sample.name]==3){
				#cat('end..')
				NAtop <- min(which(rev(chbit$values==3)))
				Itop <- end - cumsum(rev(chbit$lengths))[NAtop]
				cnv.matrix[(Itop+1):end,sample.name] <- cnv.matrix[Itop,sample.name]
			}
		}

		nav <- which(chbit$values[-c(1,length(chbit$values))]==3)+1
		if(length(nav>0)){
			#cat('gaps..')
			# find nearest neighbors of inner-chromosome gaps
			gstarts <- cumsum(chbit$lengths)[nav-1]+start
			glengths <- chbit$lengths[nav]

			nrows <- 
				gstarts %>%
				map2(glengths, ~ .x:(.x+.y-1)) %>%
				unlist %>%
				unname

			nfills <-
					cnv.matrix %>%
					select(chrom,start,end,hgnc,get(sample.name)) %>%
					tbl_df %>%
					slice(nrows) %>%
					mutate(mid=rowMeans(.[,c('start','end')])) %>%
					rowwise %>%
					inner_join(
						cnv.matrix %>%
						filter(chrom==chromosome) %>%
						mutate(qmid=rowMeans(.[,c('start','end')])) %>% 
						select(chrom,qmid,fill=get(sample.name)),
						by='chrom'
						) %>%
					ungroup %>%
					group_by(chrom,start,end,hgnc) %>%
					filter(fill!=3) %>%
					slice(which.min(abs(mid-qmid))) %>%
					.$fill

			# write neighbor calls to master table
			if(length(nfills)>0){
				cnv.matrix[nrows,sample.name] <- nfills
			}else{
				cnv.matrix[nrows,sample.name] <- NA
			}
		}
	} # loop over chromosomes

	#cat('\n')
}

cnv.matrix$chrom <- chrom
cnv.matrix %<>% mutate(mid=(start+end)/2) %>% select(chrom,start,mid,end,everything()) %>% ungroup



OverCall <- function(sample.name){

	sample.rle <- cnv.matrix[,sample.name] %>% unlist %>% unname %>% rle

	mark.amp <- intersect(which(sample.rle$lengths > 100),which(sample.rle$values == 2))

	if(length(mark.amp)!=0){

		loc.amp <- mapply(function(x,y) x:y ,c(0,cumsum(sample.rle$lengths))[mark.amp]+1 , cumsum(sample.rle$lengths)[mark.amp]) %>% unlist %>% as.vector

		cncf <-
			read.delim(str_c(cncf.path,'/',sample.name %>% strsplit('_') %>% unlist %>% head(2) %>% paste(collapse='_'),'.cncf.txt'),sep='\t',stringsAsFactors=FALSE) %>%
			tbl_df %>%
			mutate(loc.mid=(loc.start+loc.end)/2) # add genome-wide stats

		logic.amp <-
			cnv.matrix[loc.amp,c('chrom','mid')] %>%
			add_rownames(var='index') %>%
			mutate(index=as.numeric(index)) %>%
			group_by(index) %>%
			left_join(cncf,by='chrom') %>%
			slice(which.min(abs(mid-loc.mid))) %>%
			ungroup %>%
			filter(cnlr.median>1.1) %>%
			.$index

		loc.amp <- loc.amp[logic.amp]
		cnv.matrix[loc.amp,sample.name] <<- 1
	}

	mark.del <- intersect(which(sample.rle$lengths > 100),which(sample.rle$values == -2))

	if(length(mark.del)!=0){

		loc.del <- mapply(function(x,y) x:y ,c(0,cumsum(sample.rle$lengths))[mark.del]+1 , cumsum(sample.rle$lengths)[mark.del]) %>% unlist %>% as.vector

		logic.del <-
			cnv.matrix[loc.del,c('chrom','mid')] %>%
			add_rownames(var='index') %>%
			mutate(index=as.numeric(index)) %>%
			group_by(index) %>%
			ungroup %>%
			left_join(cncf,by='chrom') %>%
			slice(which.min(abs(mid-loc.mid))) %>%
			filter(cnlr.median<0.5) %>%
			.$index

		loc.del <- loc.del[logic.del]
		cnv.matrix[loc.del,sample.name] <<- -1
	}

}


names(cnv.matrix) %>%
list.filter(!. %in% c('chrom','start','mid','end','hgnc','band')) %>%
list.map(.,.) %>%
map(~ OverCall(.x))

write_tsv(cnv.matrix,'facets/geneCN_fill.txt')

cat(green("\n  [done]\n\n"))




