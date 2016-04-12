#!/usr/bin/env Rscript
# fill facets gene CN file
# load base libraries
suppressMessages(pacman::p_load(optparse,
                                RColorBrewer,
                                GenomicRanges,
                                plyr,
                                dplyr,
                                tibble,
                                readr,
                                stringr,
                                tidyr,
                                purrr,
                                magrittr,
                                rlist,
                                crayon,
                                foreach,
                                Cairo,
                                RMySQL,
                                rtracklayer,
                                colorspace,
                                ggplot2,
                                grid,
                                gridExtra,
                                RColorBrewer))
suppressPackageStartupMessages(library("facets", lib.loc="/home/bermans/R-dev/"));

#--------------
# parse options
#--------------

optList <- list(
				make_option("--geneCNFile", default = NULL, help = "gene copy-number file"),
				make_option("--outFile", default = NULL, help = "output file"))

parser <- OptionParser(usage = "%prog [options] [facets files]", option_list = optList);

arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (length(arguments$args) < 1) {
	cat("Need cncf files\n")
	print_help(parser);
	stop();
} else if (is.null(opt$geneCNFile)) {
	cat("Need gene copynumber file\n")
	print_help(parser);
	stop();
} else if (is.null(opt$outFile)) {
	cat("Need output file\n")
	print_help(parser);
	stop();
} else {
	cncf.files <- arguments$args
}

cnv.matrix <- read.table(opt$geneCNFile, sep='\t', header=T, check.names=F)

cncfs <- lapply(cncf.files, function(x) {
                read.delim(x, sep='\t', check.names=F, stringsAsFactors=FALSE) %>%
                tbl_df %>%
                mutate(loc.mid=(loc.start+loc.end)/2) # add genome-wide stats
                })
names(cncfs) <- sub('\\.cncf\\.txt', '', basename(cncf.files))

cnv.matrix %<>%
	mutate(chrom=
		ifelse(chrom=='X', 23,
		ifelse(chrom=='Y', 23, chrom))) %>%
	mutate(chrom=as.integer(chrom))

chrom <- cnv.matrix$chrom

cnv.matrix %<>% arrange(chrom, start, end)

cnv.matrix[is.na(cnv.matrix)] <- 3  # replace NA with numeric for use in rle function
breaks <- c(0, cumsum(table(chrom)[chrom %>% table %>% names %>% as.numeric %>% order])) # start-1 == end

for (sample.name in names(cnv.matrix) %>% list.filter(!. %in% c('chrom', 'start', 'end', 'hgnc', 'band'))) {

	cat('\n * sample:', sample.name)

	#loop over chromosomes in junction table
	for (chromosome in 1:length(cnv.matrix$chrom %>% unique)) {

		#cat('\n   chr', formatC(chromosome, width=2, flag='0'), ': ', sep='')
		start <- breaks[chromosome]+1
		end <- breaks[chromosome+1]
		chbit <- cnv.matrix[start:end, sample.name] %>% unlist %>% rle

		# replace an entirely empty chromosome with calls of 0
		if ((chbit$lengths %>% length)==1) {
			cnv.matrix[start:end, sample.name] <- 0
		}else{
			# extend chromosome start calls from nearest integer on same chromosome
			if (cnv.matrix[start, sample.name]==3) {
				#cat('start..')
				NAbottom <- min(which(chbit$values==3))
				Ibottom <- start + cumsum(chbit$lengths)[NAbottom]
				cnv.matrix[start:(Ibottom-1), sample.name] <- cnv.matrix[Ibottom, sample.name]
			}
			# extend chromosome end calls from nearest integer on same chromosome
			if (cnv.matrix[end, sample.name]==3) {
				#cat('end..')
				NAtop <- min(which(rev(chbit$values==3)))
				Itop <- end - cumsum(rev(chbit$lengths))[NAtop]
				cnv.matrix[(Itop+1):end, sample.name] <- cnv.matrix[Itop, sample.name]
			}
		}

		nav <- which(chbit$values[-c(1, length(chbit$values))]==3)+1
		if (length(nav>0)) {
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
					select(chrom, start, end, hgnc, get(sample.name)) %>%
					tbl_df %>%
					slice(nrows) %>%
					mutate(mid=rowMeans(.[, c('start', 'end')])) %>%
					rowwise %>%
					inner_join(
						cnv.matrix %>%
						filter(chrom==chromosome) %>%
						mutate(qmid=rowMeans(.[, c('start', 'end')])) %>% 
						select(chrom, qmid, fill=get(sample.name)),
						by='chrom'
						) %>%
					ungroup %>%
					group_by(chrom, start, end, hgnc) %>%
					filter(fill!=3) %>%
					slice(which.min(abs(mid-qmid))) %>%
					.$fill

			# write neighbor calls to master table
			if (length(nfills)>0) {
				cnv.matrix[nrows, sample.name] <- nfills
			}else{
				cnv.matrix[nrows, sample.name] <- NA
			}
		}
	} # loop over chromosomes

	#cat('\n')
}

cnv.matrix$chrom <- chrom
cnv.matrix %<>% mutate(mid=(start+end)/2) %>% select(chrom, start, mid, end, everything()) %>% ungroup



OverCall <- function(column.name) {
    sample.name <- sub('_EM$', '', column.name)
    sample.name <- sub('_LRR_threshold$', '', sample.name)

	sample.rle <- cnv.matrix[, column.name] %>% unlist %>% unname %>% rle

	mark.amp <- intersect(which(sample.rle$lengths > 100), which(sample.rle$values == 2))

	if (length(mark.amp)!=0) {

		loc.amp <- mapply(function(x, y) x:y , c(0, cumsum(sample.rle$lengths))[mark.amp]+1 , cumsum(sample.rle$lengths)[mark.amp]) %>% unlist %>% as.vector


		logic.amp <-
			cnv.matrix[loc.amp, c('chrom', 'mid')] %>%
			rownames_to_column(var='index') %>%
			mutate(index=as.numeric(index)) %>%
			group_by(index) %>%
			left_join(cncfs[[sample.name]], by='chrom') %>%
			slice(which.min(abs(mid-loc.mid))) %>%
			ungroup %>%
			filter(cnlr.median>1.1) %>%
			.$index

		cnv.matrix[logic.amp, column.name] <<- 1
	}

	mark.del <- intersect(which(sample.rle$lengths > 100), which(sample.rle$values == -2))

	if (length(mark.del)!=0) {

		loc.del <- mapply(function(x, y) x:y , c(0, cumsum(sample.rle$lengths))[mark.del]+1 , cumsum(sample.rle$lengths)[mark.del]) %>% unlist %>% as.vector

		logic.del <-
			cnv.matrix[loc.del, c('chrom', 'mid')] %>%
			rownames_to_column(var='index') %>%
			mutate(index=as.numeric(index)) %>%
			group_by(index) %>%
			ungroup %>%
			left_join(cncfs[[sample.name]], by='chrom') %>%
			slice(which.min(abs(mid-loc.mid))) %>%
			filter(cnlr.median<0.5) %>%
			.$index

		cnv.matrix[logic.del, column.name] <<- -1
	}

}


names(cnv.matrix) %>%
list.filter(!. %in% c('chrom', 'start', 'mid', 'end', 'hgnc', 'band')) %>%
list.map(., .) %>%
map(~ OverCall(.x))

write_tsv(cnv.matrix, opt$outFile)

cat(green("\n  [done]\n\n"))



