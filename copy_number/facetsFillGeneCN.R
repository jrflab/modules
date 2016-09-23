#!/usr/bin/env Rscript

# fill facets gene CN file

#--------------------
# load base libraries
#--------------------

for (lib in c("optparse","RColorBrewer","GenomicRanges","plyr","dplyr","stringr","tidyr","magrittr","foreach","Cairo","rtracklayer","grid","rlist","purrr","tibble","readr","crayon")) {
    suppressPackageStartupMessages(library(lib, character.only=TRUE))
}

# use showtext if installed for monospaced system font, useful for TCGA
# barcodes on y axis
if ('showtext' %in% rownames(installed.packages())) {
  suppressMessages(library(showtext))
  font.add('DejaVuSansMono', 'DejaVuSansMono.ttf')
  showtext.auto()
  fontfamily <- 'DejaVuSansMono'
} else {
  fontfamily <- 'serif'
}

#--------------
# parse options
#--------------

if(!interactive()) {

    optList <- list( make_option("--geneCN_txt", default = NULL, help = "gene copy number input file"),
                     make_option("--geneCN_fill_txt", default = NULL, help = "gene copy number fill output file"),
                     make_option("--geneCN_pdf", default = NULL, help = "gene copy number heatmap plot"),
                     make_option("--geneCN_fill_pdf", default = NULL, help = "gene copy number fill heatmap plot") )

    parser <- OptionParser(usage = "%prog [options] [facets files]", option_list = optList)
    arguments <- parse_args(parser, positional_arguments = T)
    opt <- arguments$options

    if (length(arguments$args) < 1) {
        message('Need cncf files')
        print_help(parser)
        stop()
    } else if (is.null(opt$geneCN_txt)) {
        message('Need gene copy number input file location')
        print_help(parser)
        stop()
    } else if (is.null(opt$geneCN_fill_txt)) {
        message('Need gene copy number fill output file location')
        print_help(parser)
        stop()
    } else if (is.null(opt$geneCN_pdf)) {
        message('Need gene copy number heatmap plot location')
        print_help(parser)
        stop()
    } else if (is.null(opt$geneCN_fill_pdf)) {
        message('Need gene copy number fill heatmap plot location')
        print_help(parser)
        stop()
    } else {
        cncf.files <- arguments$args
    }

} else {
    opt <- list( geneCN_txt = 'facets/geneCN.txt',
                 geneCN_fill_txt = 'facets/geneCN.fill.txt',
                 geneCN_pdf = 'facets/geneCN.pdf',
                 geneCN_fill_pdf = 'facets/geneCN.fill.pdf' )
    cncf.files <- list.files('facets/cncf', pattern='*cncf.txt', full.names=TRUE)
}


#---------------------
# CNA heatmap function
#---------------------

PlotCNHeatmap <- function(gene.cn, file.name, sample.names=NULL, threshold=FALSE) {

    if(all(class(gene.cn) == 'character')) {
      gene.cn <- read.delim(gene.cn, sep='\t', stringsAsFactors=FALSE) %>% tbl_df
    }

    if(is.null(sample.names) & threshold==TRUE) {
        sample.names <- gene.cn %>% select(matches('threshold')) %>% names %>% sort
    } else if(is.null(sample.names)) {
        sample.names <- gene.cn %>% names %>% list.filter(! . %in% c('hgnc','gene','chrom','start','mid','end','band'))
    }

    # convert X & Y to 23
    gene.cn %<>% mutate(chrom=as.numeric(ifelse(chrom == 'X' | chrom == 'Y', '23', chrom)))

    # segment spans
    chr.rle <- gene.cn$chrom %>% rle
    chr.sep <- chr.rle$lengths %>% cumsum

    # midpoints
    chr.mid <- c(0, chr.sep[-length(chr.sep)]) + chr.rle$lengths/2

    # remove annotation cols
     gene.cn %<>% select(one_of(rev(sample.names)))

    pdf(file.name, width=24, height=4+length(sample.names)/2)

        par(mar=c(14, 14, 1, 1), oma=c(1, 1, 1, 1))  # bottom, left, top, right

        image(as.matrix(gene.cn), col=c('#cf3a3d', '#dc9493', '#FFFFFF', '#7996ba', '#2a4b94'), xaxt='n', yaxt='n', zlim=c(-2, 2))

        for (i in seq(-1, max(((2*(ncol(gene.cn)-1))+1),1), 2)) {
            abline(h=i/(2*(ncol(gene.cn)-1)), col="white", lwd=2)
        }

        for (i in (chr.sep*2)-1) {
            abline(v=i/((max(chr.sep)-1)*2), lwd=1.5, col='grey', lty="dashed")
        }

         axis( 1,
               at       = chr.mid/(max(chr.sep)-1),
               label    = chr.rle$values,
               cex.axis = 1.3,
               tick     = FALSE )

         axis( 2,
               at       = if(ncol(gene.cn)==1){ 0.5 }else{ seq(0, 1, 1/max((ncol(gene.cn)-1),1)) },
               label    = if(threshold==TRUE){ sub('_LRR_threshold$', '', colnames(gene.cn)) }else{ colnames(gene.cn) },
               las      = 2,
               cex.axis = 1.1,
               tick     = FALSE )

        box()

         legend( 'bottom',
                 inset  = c(0, -0.5),
                 legend = c('Amplification', 'Gain', 'Loss', 'Homozygous deletion'),
                 fill   = c('#2a4b94', '#7996ba', '#dc9493', '#cf3a3d'),
                 xpd    = TRUE,
                 ncol   = 2,
                 cex    = 1.1 )
     dev.off() %>% invisible

}


#--------------------
# CNA overcalling fix
#--------------------

OverCall <- function(column.name) {
    sample.name <- sub('_EM$', '', column.name)
    sample.name <- sub('_LRR_threshold$', '', sample.name)

    sample.rle <- cnv.matrix[, column.name] %>% unlist %>% unname %>% rle

    mark.amp <- intersect(which(sample.rle$lengths > 50), which(sample.rle$values == 2))

    if (length(mark.amp)!=0) {

        loc.amp <- mapply(function(x, y) x:y , c(0, cumsum(sample.rle$lengths))[mark.amp]+1 , cumsum(sample.rle$lengths)[mark.amp]) %>% unlist %>% as.vector


        logic.amp <-
            cnv.matrix[loc.amp, c('chrom', 'mid')] %>%
            rownames_to_column(var='index') %>%
            mutate(index=as.numeric(index)) %>%
            group_by(index) %>%
            left_join(cncfs[[sample.name]] %>% mutate(chrom=ifelse(chrom=='Y' | chrom=='X', 23, chrom)), by='chrom', copy=TRUE) %>%
            dplyr::slice(which.min(abs(mid-loc.mid))) %>%
            ungroup %>%
            filter(cnlr.median>0.8) %>%
            .$index

        cnv.matrix[logic.amp, column.name] <<- 1
    }

    mark.del <- intersect(which(sample.rle$lengths > 50), which(sample.rle$values == -2))

    if (length(mark.del)!=0) {

        loc.del <- mapply(function(x, y) x:y , c(0, cumsum(sample.rle$lengths))[mark.del]+1 , cumsum(sample.rle$lengths)[mark.del]) %>% unlist %>% as.vector

        logic.del <-
            cnv.matrix[loc.del, c('chrom', 'mid')] %>%
            rownames_to_column(var='index') %>%
            mutate(index=as.numeric(index)) %>%
            group_by(index) %>%
            ungroup %>%
            left_join(cncfs[[sample.name]] %>% mutate(chrom=ifelse(chrom=='Y' | chrom=='X', 23, chrom)), by='chrom', copy=TRUE) %>%
            dplyr::slice(which.min(abs(mid-loc.mid))) %>%
            filter(cnlr.median<0.4) %>%
            .$index

        cnv.matrix[logic.del, column.name] <<- -1
    }

}


#------------
# facets fill
#------------

cnv.matrix <- read.delim(opt$geneCN_txt, sep='\t', check.names=FALSE, stringsAsFactors=FALSE)

cncfs <- lapply(cncf.files, function(x) {

                cncf <- read.delim(x, sep='\t', check.names=F, stringsAsFactors=FALSE) %>%
                tbl_df
                # rename column names if they exist
                names.key <- c('loc.start'='start', 'loc.end'='end')
                existing <- match(names(names.key), names(cncf))
                names(cncf)[na.omit(existing)] <- names.key[which(!is.na(existing))]
                # add midpoints
                cncf %<>% mutate(loc.mid=(start+end)/2) # add genome-wide stats
                return(cncf)
                })
names(cncfs) <- sub('\\.cncf\\.txt', '', basename(cncf.files))

cnv.matrix %<>%
    mutate(chrom=
        ifelse(chrom=='X', 23,
        ifelse(chrom=='Y', 23, chrom))) %>%
    mutate(chrom=as.integer(chrom))

PlotCNHeatmap(gene.cn=cnv.matrix, file.name=opt$geneCN_pdf, threshold=TRUE)

cnv.matrix[is.na(cnv.matrix)] <- 3  # replace NA with numeric for use in rle function
breaks <- c(0, cumsum(table(cnv.matrix$chrom)[cnv.matrix$chrom %>% table %>% names %>% as.numeric %>% order])) # start-1 == end

# main loop
for (sample.name in names(cnv.matrix) %>% list.filter(!. %in% c('chrom', 'start', 'end', 'hgnc', 'band'))) {

    message(' * sample: ', sample.name)

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
                        by='chrom', copy=TRUE) %>%
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
    }  # loop over chromosomes
}


# add midpoints
cnv.matrix %<>% mutate(mid=(start+end)/2) %>% select(chrom, start, mid, end, band, hgnc, everything())

# call overcall fix function
names(cnv.matrix) %>%
list.filter(!. %in% c('chrom', 'start', 'mid', 'end', 'hgnc', 'band')) %>%
list.map(., .) %>%
map(~ OverCall(.x)) %>% invisible

# write fill file
write_tsv(cnv.matrix, opt$geneCN_fill_txt)

# plot fill heatmap
PlotCNHeatmap(cnv.matrix, file.name=opt$geneCN_fill_pdf, threshold=TRUE) %>% invisible


message(green(' [done]'))
