#!/usr/bin/env Rscript
# fill facets gene CN file

# load base libraries
pacman::p_load( optparse,
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
                RColorBrewer )

suppressPackageStartupMessages(library("facets", lib.loc="/home/bermans/R-dev/"));

#--------------
# parse options
#--------------

if(!interactive()) {

    optList <- list( make_option("--geneCNFile", default = NULL, help = "gene copy-number file"),
                     make_option("--outFile", default = NULL, help = "output file") )

    parser <- OptionParser(usage = "%prog [options] [facets files]", option_list = optList)
    arguments <- parse_args(parser, positional_arguments = T)
    opt <- arguments$options

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

} else {
    opt <- list( geneCNFile = 'facets/geneCN.txt',
                 outFile    = 'facets/geneCN.fill.txt' )
    cncf.files <- list.files('facets/cncf', pattern='*cncf.txt', full.names=TRUE)
}


#---------------------
# CNA heatmap function
#---------------------

PlotCNHeatmap <- function(gene.cn, file.name, sample.names=NULL, threshold=FALSE) {

    if(is.null(sample.names) & threshold==TRUE) {
        sample.names <- gene.cn %>% select(matches('threshold')) %>% names %>% sort
    } else if(is.null(sample.names)) {
        sample.names <- gene.cn %>% names %>% list.filter(! . %in% c('hgnc','gene','chrom','start','mid','end','band')) %>% sort
    }

    chr.rle <- gene.cn$chrom %>% rle
    chr.sep <- chr.rle$lengths %>% cumsum
    chr.mid <- c(0, chr.sep[-length(chr.sep)]) + chr.rle$lengths/2

    pdf(file.name, width=12, height=3 + 0.25*length(sample.names))

        g.cn <- gene.cn %>% select(one_of(rev(sample.names)))

        layout(matrix(c(0,1),2,1,byrow=TRUE), c(length(sample.names),1), TRUE)  

        par(mfrow=c(1,1), mar=c(8,5,1,1))
        image(as.matrix(g.cn), col=c('#CF3A3D', '#DC9493', '#FFFFFF', '#7996BA', '#2A4B94'), xaxt='n', yaxt='n', zlim=c(-2, 2))

        for (i in seq(-1, max(((2*(ncol(g.cn)-1))+1),1), 2)) {
            abline(h=i/(2*(ncol(g.cn)-1)), col="white", lwd=2)
        }

        for (i in (chr.sep*2)-1) {
            abline(v=i/((max(chr.sep)-1)*2), lwd=1.5, col='grey', lty="dashed")
        }

        box()

         axis( 1,
               at       = chr.mid/(max(chr.sep)-1),
               label    = chr.rle$values,
               cex.axis = 0.8,
               tick     = FALSE )

         axis( 2,
               at       = if(ncol(g.cn)==1){0.5}else{seq(0, 1, 1/max((ncol(g.cn)-1),1))},
               label    = sub('T_.*', '', colnames(g.cn)),
               las      = 2,
               cex.axis = 1,
               tick     = FALSE )

         legend( 'bottom',
                 inset  = c(0, -0.28),
                 legend = c('Homozygous deletion', 'Loss', 'Gain', 'Amplification'),
                 fill   = c('#CF3A3D', '#DC9493', '#7996BA', '#2A4B94'),
                 xpd    = TRUE,
                 ncol   = 2 )
     dev.off()

}


#--------------------
# CNA overcalling fix
#--------------------

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
            dplyr::slice(which.min(abs(mid-loc.mid))) %>%
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
            dplyr::slice(which.min(abs(mid-loc.mid))) %>%
            filter(cnlr.median<0.5) %>%
            .$index

        cnv.matrix[logic.del, column.name] <<- -1
    }

}


#------------
# facets fill
#------------

cnv.matrix <- read.delim(opt$geneCNFile, sep='\t', check.names=FALSE, stringsAsFactors=FALSE)

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

PlotCNHeatmap(cnv.matrix, file.name='facets/gene.cn.pdf', threshold=TRUE)

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
                    select(chrom, start, end, hgnc, one_of(sample.name)) %>%
                    tbl_df %>%
                    filter(row_number() %in% nrows) %>%
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
                    dplyr::slice(which.min(abs(mid-qmid))) %>%
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


# plot final heatmap
PlotCNHeatmap(cnv.matrix, file.name='facets/gene.cn.fill.pdf', threshold=TRUE)

# call overcall fix function
names(cnv.matrix) %>%
list.filter(!. %in% c('chrom', 'start', 'mid', 'end', 'hgnc', 'band')) %>%
list.map(., .) %>%
map(~ OverCall(.x))

write_tsv(cnv.matrix, opt$outFile)

cat(green("\n  [done]\n\n"))
