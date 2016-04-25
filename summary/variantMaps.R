#!/usr/bin/env Rscript

# script for build mutation heatmaps, trees & other standard plots

#-----------#
#           #
# LIBRARIES #
#           #
#-----------#

pacman::p_load( dplyr,readr,tidyr,magrittr,purrr,stringr,rlist,openxlsx, # base
                crayon,colorspace,RColorBrewer, # coloring
                ggplot2,grid,gridExtra,gplots, # plot layout
                dendextend,dendextendRcpp,dynamicTreeCut,gclus,phangorn, # dentrogram
                digest ) # hashing

source('modules/summary/variantCN.R')

#----------------#
#                #
# OPTION PARSING #
#                #
#----------------#

optList <- list( make_option("--mutationSummary", default = NULL, help = "mutation summary table"),
                 make_option("--geneCN", default = NULL, help = "filled geneCN file"),
                 make_option("--mutOutFile", default = NULL, help = "mutation output file"),
                 make_option("--mutRecurrentOutFile", default = NULL, help = "recurrent mutation output file"),
                 make_option("--cnOutFile", default = NULL, help = "copy number output file"),
                 make_option("--cnRecurrentOutFile", default = NULL, help = "recurrent copy number output file"),
                 make_option("--cnAmpDelOutFile", default = NULL, help = "amp del table") )

parser <- OptionParser(usage = "%prog [options] [mutation summary file]", option_list = optList)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

if (is.null(opt$geneCN)) {
    cat("Need geneCN file\n")
    print_help(parser);
    stop();
else if (is.null(opt$mutationSummary)) {
    cat("Need mutation summary file\n")
    print_help(parser);
    stop();
} else if (is.null(opt$mutOutFile)) {
    cat("Need mut output file\n")
    print_help(parser);
    stop();
} else if (is.null(opt$mutRecurrentOutFile)) {
    cat("Need mut recurrent output file\n")
    print_help(parser);
    stop();
} else if (is.null(opt$cnOutFile)) {
    cat("Need cn output file\n")
    print_help(parser);
    stop();
} else if (is.null(opt$cnRecurrentOutFile)) {
    cat("Need cn recurrent output file\n")
    print_help(parser);
    stop();
} else if (is.null(opt$cnAmpDeltOutFile)) {
    cat("Need cn amp del output file\n")
    print_help(parser);
    stop();
} else {
    muts.file <- opt$mutationSummary
    cn.file <- opt$geneCN
    muts.out.file <- opt$mutOutFile
    muts.recurrent.out.file <- opt$mutRecurrentOutFile
    cn.out.file <- opt$cnOutFile
    cn.recurrent.out.file <- opt$cnRecurrentOutFile
}

# set graphics device
options(device=pdf)
#options(encoding='ISOLatin2.enc')
#pdf.options(encoding='ISOLatin2.enc')


#--------------------#
#                    #
# INPUT & PARAMETERS #
#                    #
#--------------------#

# run options
include.silent       = FALSE
dist.method          = 'hamming'   # 'hamming', euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'     | see ?dist for details [hamming method implemented manually]
clust.method         = 'complete' # 'complete', ward', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid' | see ?hclust for details
exclude.values       = c('', '.', 'Normal', 'Not Performed', 'Performed but Not Available', FALSE)
tree.sort            = 'distance'  # ladder 'ladderize'
tree.labels          = TRUE
pheno                = NULL
color.seed           = 3
color.clusters       = TRUE
min.cluster          = 3
random.pheno.color   = TRUE
pheno.palette        = c('#2d4be0', '#20e6ae', '#ccb625', '#969696')
cn.cols              = 'threshold'  # 'threshold' = reserved string for selecting columns with threshold sufix, 'all' = reserved string for using all columns, else a string specifing columns to use
use.keys             = FALSE

# read keys file
if('keys.yaml' %>% file.exists){
    keys <- list.load('keys.yaml') %>% unlist
    values <- names(keys) %>% setNames(keys)
}
if(use.keys==FALSE){keys=NULL}

# read subsets file
if('subsets.yaml' %>% file.exists){
    if(use.keys !=FALSE){
        subsets <- list.load('subsets.yaml') %>% map(~ { keys[.x]})
    } else {
        subsets <- list.load('subsets.yaml')
    }
}

# read subset_groups.yaml file
if('subset_groups.yaml' %>% file.exists){
    subset.groups <- list.load('subset_groups.yaml')
}

# set detaults if not supplied (for interactive use)
if(!exists('muts.file')){muts.file = 'summary/mutation_summary.xlsx'}
if(!exists('muts.out.file')){muts.out.file = 'summary/mutation_heatmap.pdf'}
if(!exists('muts.recurrent.out.file')){muts.out.file = 'summary/mutation_heatmap_recurrent.pdf'}
if(!exists('cn.file')){cn.file = 'facets/geneCN_fill.txt'}
if(!exists('cn.out.file')){cn.out.file = 'summary/cn_heatmap.pdf'}
if(!exists('cn.recurrent.out.file')){cn.out.file = 'summary/cn_heatmap_recurrent.pdf'}
if(!exists('cn.amp.del.out')){cn.amp.del.out = 'summary/cn_amp_del.tsv'}


#-----------#
#           #
# FUNCTIONS #
#           #
#-----------#

# Hamming distance function
Hamming <- function(event.matrix) {
    D <- (1 - event.matrix) %*% t(event.matrix)
    D + t(D)
}


# wrap Hamming function, for easy specification
DistExtra <- function(event.matrix,dist.method){
    if(dist.method=='hamming'){
        dist.matrix <-
            event.matrix %>%
            Hamming %>%
            as.dist
    }else{
        dist.matrix <-
            event.matrix %>%
            dist(method=dist.method)
    }
    return(dist.matrix)
}


# character hashing
hash <- function(x) {
    hstr <- digest(x, algo='xxhash32')
    as.numeric(paste0('0x', hstr)) %% 80
}


#------------------------------------
# rename columns & clean nomenclature
#------------------------------------

FormatMuts <- function(muts, col.names=NULL, drop=FALSE, allosome='merge', keys=NULL) {

# col.names:
#
#   NULL = use only the columns available in input table 
#   character string = dumy column names to be added if absent
#   'all' = reserved character string to specify addition of all maintained columns
#
#   drop = TRUE / FALSE, should additional columns be dropped
#
# allosome:
#
#   'none'      = exclude X & Y chromosome
#   'merge'     = treat X & Y coordinates as homologous pair
#   'distinct'  = treat X & Y as seperate, sequential chromosomes

    # lookup table
    col.keys <- c( 'alt'='alt', 'ALT'='alt', 'Alternate.allele'='alt',
                   'band'='band', 'Band'='band',
                   'cancer.gene'='cancer.gene', 'Cancer.Gene.Census'='cancer.gene', 'Cancer Gene Census'='cancer.gene','Cancer5000.S.genes..Lawrence.et.al.'='cancer.gene',
                   'ccf'='ccf', 'CCF'='ccf', 'cancer_cell_frac'='ccf', 'Cancer cell fraction'='ccf', 'Cancer.cell.fraction'='ccf',
                   'chasm'='chasm', 'CHASM'='chasm',
                   'chrom'='chrom','Chrom'='chrom','Chromosome'='chrom','CHROM'='chrom',
                   'clonality'='clonality', 'Clonality'='clonality',
                   'cn'='cn', 'CN'='cn',
                   'ci95.low'='ci95.low', 'ccf_CI95_low'='ci95.low',
                   'effect'='effect','Hugo_Symbol'=='gene', 'Effect'='effect', 'Variant_Classification'='effect', 'ANN....EFFECT'='effect', 'ANN[*].EFFECT'='effect',
                   'end'='end', 'end'='stop',
                   'fathmm'='fathmm', 'FATHMM'='fathmm',
                   'gene'='gene', 'Gene'='gene', 'Hugo_Symbol'='gene', 'GENE'='gene', 'hgnc'='gene', 'ANN....GENE'='gene', 'ANN[*].GENE'='gene', 'Gene.symbol'='gene',
                   'haploinsufficient'='haploinsufficient',
                   'kandoth'='kandoth', '127 significantly mutated genes (Kandoth et al)'='kandoth', '127 significantly.mutated genes.(Kandoth et al)'='kandoth','X127.significantly.mutated.genes..Kandoth.et.al.'='kandoth',
                   'lawrence'='lawrence', 'Cancer5000-S genes (Lawrence et al)'='lawrence', 'Cancer5000-S genes.(Lawrence et al)'='lawrence',
                   'loh'='loh', 'LOH'='loh', 'Loss.of.heterozygocity.(LOH)'='loh', 'Loss of heterozygocity (LOH)'='loh','Loss.of.heterozygocity..LOH.'='loh',
                   'maf'='maf','Mutant allele fraction','maf','Mutant.allele.fraction'='maf',
                   'mut.taster'='mut.taster',
                   'pathogenic'='pathogenic', 'Pathogenic'='pathogenic',
                   'pheno'='pheno', 'pheno.bar'='pheno',
                   'pos'='pos','POS'='pos', 'Position'='pos', 'position'='position',
                   'pr.somatic.clonal'='pr.somatic.clonal', 'Pr_somatic_clonal'='pr.somatic.clonal',
                   'provean'='provean', 'Provean'='provean',
                   'purity'='purity',
                   'ref'='ref', 'REF'='ref', 'Reference.allele'='ref',
                   'sample'='sample', 'Sample'='sample', 'Tumor_Sample_Barcode'='sample', 'TUMOR_SAMPLE'='sample', 'Sample.ID'='sample',
                   'start'='start', 'Start'='start', 'Start_position'='start' )

    # function to dummy column if absent
    DummyCols <- function(muts, col.names) {
        for(col.name in col.names) {
            if(!col.name %in% colnames(muts)) {
                message(green(str_c('adding dummy column: ', col.name)))
                muts.names <- colnames(muts) %>% replace(which(is.na(.)),'empty') %>% list.filter(.!='add.col')
                muts$add.col <- NA
                muts %<>% setNames(c(muts.names,col.name))
            }
        }
        return(muts)
    }

    # convert chrom vector to desired format
    ChromMod <- function(muts, allosome){

        if('X' %in% muts$chrom) { message(yellow('X chromosome labelling found')) }
        if('Y' %in% muts$chrom) { message(yellow('Y chromosome labelling found')) }

        if(allosome=='distinct') {
                muts %>%
                mutate(.,chrom=as.numeric(ifelse(chrom=='X', 23,chrom))) %>%
                mutate(chrom=as.numeric(ifelse(chrom=='Y', 24,chrom)))
        } else if(allosome =='merge') {
                muts %>% 
                mutate(chrom=as.numeric(ifelse(chrom %in% c('X','Y'), 23, chrom)))
        } else {
                muts %>%
                filter(!chrom %in% c('X','Y'))
        }
    }

    # save column names to fill unhandled
    fallback.cols <- colnames(muts)

    # rename columns
    names(muts) <- col.keys[names(muts)]

    # specify output columns
    if(is.null(col.names)) {
        col.names <- colnames(muts)
        cols.match <- col.names %>% list.filter(!is.na(.))
        col.names[which(is.na(col.names))] <- fallback.cols[which(is.na(col.names))]
        names(muts) <- col.names
    } else if(col.names=='all') {
        col.names <- c( 'alt','band','cancer.gene','ccf','chasm','chrom','clonality','cn','ci95.low','effect',
                        'end','fathmm','gene','haploinsufficient','kandoth','lawrence','loh','maf','mut.taster',
                        'pathogenic','pheno','pos','pr.somatic.clonal','provean','purity','ref','sample','start' )
    }

    TypeCol <- function(muts, columns, types) {
        for (col in 1:length(columns)) {
            if(columns[col] %in% colnames(muts)) {
                if(types[col] == 'char') {
                    muts[,columns[col]] <- as.character(muts[,columns[col]])
                } else if(type == 'num') {
                    muts[,columns[col]] <- as.numeric(muts[,columns[col]])
                } else if(type == 'logical') {
                    muts[,columns[col]] <- as.logical(muts[,columns[col]])
                }
            }
        }
        return(muts)
    }

    # df typing to avoid row name conflicts
    muts %<>% as.data.frame

    # add columns if absent
    muts %<>% DummyCols(col.names)

    # select columns
    muts <- muts[colnames(muts) %>% list.filter(.!='empty')]

    # chromosome type conversion
    muts %<>% ChromMod(allosome)

    # column typing
    muts %<>% TypeCol(c('sample','gene','effect','ref','alt'), c('char','char','char','char','char'))

    if('effect' %in% colnames(muts)){
        # rename variant classifications
        muts %<>%
            tbl_df %>%
            { muts <- .
                if(!all(is.na(muts$effect))){ filter(muts,!is.na(effect)) }
                muts
            } %>%
            unique %>%
            mutate(effect=
                ifelse(effect%in%c('STOP_GAINED','Nonsense_Mutation','stop_gained&splice_region_variant','stop_gained','Nonsense_Mutation','Stop_Codon_Ins','nonsense','truncating snv','Truncating snv','Truncating snv','Truncating SNV'),'Truncating SNV',
                ifelse(effect%in%c('FRAME_SHIFT','FRAME_SHIFT','Frame_Shift_Del','Frame_Shift_Ins','frameshift_variant','frameshift_variant&stop_gained','frameshift_variant&splice_region_variant','frameshift_variant&splice_acceptor_variant&splice_region_variant&splice_region_variant&intron_variant','Frame_Shift_Del','Frame_Shift_Ins','frame_shift_del','frame_shift_ins','frameshift indel','Frameshift indel','Frameshift In-Del','frameshift_variant'),'Frameshift In-Del',
                ifelse(effect%in%c('NON_SYNONYMOUS_CODING','STOP_LOST','Missense_Mutation','missense_variant','missense_variant&splice_region_variant','missense_variant|missense_variant','Missense_Mutation','missense','missense snv','Missense snv','Missense SNV'),'Missense SNV',
                ifelse(effect%in%c('CODON_CHANGE_PLUS_CODON_DELETION','CODON_DELETION','CODON_INSERTION','In_Frame_Ins','In_Frame_Del','disruptive_inframe_deletion','disruptive_inframe_insertion','inframe_deletion','inframe_insertion','disruptive_inframe_deletion&splice_region_variant','inframe_deletion&splice_region_variant','In_Frame_Del','In_Frame_Ins','in_frame_del','in_frame_ins','inframe indel','Inframe indel','Inframe In-Del'),'Inframe In-Del',
                ifelse(effect%in%c('SPLICE_SITE_DONOR','SPLICE_SITE_ACCEPTOR','SPLICE_SITE_REGION','Splice_Site','splice_donor_variant&intron_variant','splice_acceptor_variant&intron_variant','splicing','splice_donor_variant&splice_region_variant&intron_variant','splice_donor_variant&disruptive_inframe_deletion&splice_region_variant&splice_region_variant&intron_variant','Splice_Site','splice','splice site variant','Splice site variant','missense_variant & splice_region_variant'),'Splice site variant',
                ifelse(effect%in%c('STOP_LOST','START_LOST','START_GAINED','UTR_5_PRIME','start_lost','stop_lost',"5'UTR","5'Flank",'De_novo_Start_InFrame','De_novo_Start_OutOfFrame','Stop_Codon_Del','Start_Codon_SNP','Start_Codon_Ins','Start_Codon_Del','Nonstop_Mutation','nonstop','upstream, start/stop, or de novo modification','Upstream, start/stop, or de novo modification'),'Upstream, start/stop, or de novo modification',
                ifelse(effect%in%c('synonymous_variant','splice_region_variant&synonymous_variant','splice_region_variant&synonymous_variant','non_coding_exon_variant','upstream_gene_variant','downstream_gene_variant','intron_variant','frameshift_variant&splice_donor_variant&splice_region_variant&splice_region_variant&intron_variant','non_coding_exon_variant|synonymous_variant','SYNONYMOUS_CODING','synonymous_variant|synonymous_variant','splice_region_variant&synonymous_variant|splice_region_variant&non_coding_exon_variant','splice_acceptor_variant & intron_variant','intragenic_variant',"3'UTR",'IGR','lincRNA','RNA','Intron','silent','intron_exon','silent','Silent','intron_variant & missense_variant'),'Silent',
                ifelse(effect%in%c('Amplification','amplification','amp','2'),'Amplification',
                ifelse(effect%in%c('Deletion','deletion','del','-2'),'Deletion',
                ifelse(is.na(effect),NA,
            NA)))))))))))

        # warn on unknown effects
        if(all(is.na(muts$effect))){
            message(yellow('warning: empty effect column, left as NA'))
        } else if(muts %>% filter(is.na(effect)) %>% nrow > 0) {
            message(yellow('warning: some variants not accounted for'))
        }
    }

    # remove unmutated LOH if present & format
    if('loh' %in% colnames(muts)) {
        muts %<>% mutate(loh=ifelse(loh%in%c('LOH','loh','Loss of heterozygosity') & !is.na(effect),'LOH',NA))
    }

    if(drop!=FALSE) {
        muts %<>% select(one_of(cols.match))
    }

    if(!is.null(keys)) {
       muts %<>% mutate(sample=keys[sample])
    }

    return(muts)
}


#----------------------------------
# prepare melted table for plotting
#----------------------------------

OrgMuts <- function(muts, sample.order=NULL, recurrence=1, allosome='merge', cn=FALSE) {

# sample.order:
#   NULL = arrange alphabetically
#   character string = specify specific order
#
# recurrence:
#   integer cutoff for recurrent mutations
#
# allosome:
#
#   'none'      = exclude X & Y chromosome
#   'merge'     = treat X & Y coordinates as homologous pair
#   'distinct'  = treat X & Y as seperate, sequential chromosomes

    # function to dummy column if absent
    DummyCols <- function(muts, col.names) {
        for(col.name in col.names) {
            if(!col.name %in% colnames(muts)) {
                message(green(str_c('adding dummy column: ',col.name)))
                muts.names <- colnames(muts) %>% replace(which(is.na(.)),'empty') %>% list.filter(.!='add.col')
                muts$add.col <- NA
                muts %<>% setNames(c(muts.names,col.name))
            }
        }
        return(muts)
    }

    # convert chrom vector to desired format
    ChromMod <- function(muts, allosome){

        if('X' %in% muts$chrom) { message(yellow('X chromosome labelling found'))}
        if('Y' %in% muts$chrom) { message(yellow('Y chromosome labelling found'))}

        if(allosome=='distinct') {
                muts %>%
                mutate(.,chrom=as.numeric(ifelse(chrom=='X', 23,chrom))) %>%
                mutate(chrom=as.numeric(ifelse(chrom=='Y', 24,chrom)))
        } else if(allosome =='merge') {
                muts %>% 
                mutate(chrom=as.numeric(ifelse(chrom %in% c('X','Y'), 23, chrom)))
        } else {
                muts %>%
                filter(!chrom %in% c('X','Y'))
        }
    }

    # specify output columns
    col.names <- c('sample','gene','chrom','effect','pheno','band','pos','maf','ccf','loh','cn','clonality','pathogenic')

    # add dummy column names
    muts %<>% DummyCols(col.names)

    # effect prescedence
    muts %<>%
        select(one_of(col.names)) %>%
        unique %>%
        filter(!is.na(effect)) %>%
        # rank variant importance for plot overlay
        mutate(precedence=
            ifelse(effect=='Deletion',1,
            ifelse(effect=='Amplification',2,
            ifelse(effect=='Truncating SNV',3,
            ifelse(effect=='Frameshift In-Del',4,
            ifelse(effect=='Missense SNV',5,
            ifelse(effect=='Inframe In-Del',6,
            ifelse(effect=='Splice site variant',7,
            ifelse(effect=='Upstreop, or de novo modification',8,
            ifelse(effect=='Silent',9,
            NA)))))))))) %>%
        # remove genes with lower prescedence
        group_by(sample,gene) %>%
        arrange(gene,precedence) %>%
        slice(rank(precedence, ties.method="first")==1) %>%
        # count number of variants per gene
        unique %>%
        group_by(gene) %>%
        mutate(n.gene=n()) %>%
        ungroup %>%
        { muts <- .
            if(recurrence>0) {muts %<>% filter(n.gene>recurrence)}
            return(muts)
        } %>%
        # count number of variants per band
        unique %>%
        group_by(band) %>%
        mutate(n.band=n()) %>%
        ungroup %>%
        # define plot gene order
        arrange(desc(n.gene),sample,desc(ccf),precedence,gene) %>%
        mutate(order=row_number()) %>%
        ungroup %>%
        select(-precedence) %>%
        unique

    # fill empty tiles
    if(cn!=TRUE) {
        muts.fill <-
            muts %>%
            select(sample,gene,effect,pheno) %>%
            group_by(pheno) %>%
            spread(gene,effect) %>%
            gather(gene,effect,-pheno,-sample) %>%
            filter(is.na(effect)) %>%
            ungroup
    } else {
        muts.fill <-
            muts %>%
            select(sample,band,effect,pheno) %>%
            group_by(pheno) %>%
            spread(band,effect) %>%
            gather(band,effect,-pheno,-sample) %>%
            filter(is.na(band)) %>%
            ungroup
    }

    # plot aesthetics
    muts %<>%
        full_join(muts.fill) %>%
        # push NAs to bottom of stack
        mutate(order=ifelse(is.na(effect),-Inf,order)) %>%
        unique %>%
        # fix plot ordering & assign gene factor levels
        arrange(!is.na(effect),desc(order)) %>%
        mutate(gene=factor(gene,levels=filter(.,!is.na(effect)) %>% .$gene %>% unique)) %>%
        # ccf binning
        mutate(ccf=
            ifelse(ccf==0,               'CCF = 0%',
            ifelse(ccf>0.00 & ccf<=0.05, '0% < CCF ≤ 5%',
            ifelse(ccf>0.05 & ccf<=0.20, '5% < CCF ≤ 20%',
            ifelse(ccf>0.20 & ccf<=0.40, '20% < CCF ≤ 40%',
            ifelse(ccf>0.40 & ccf<=0.60, '40% < CCF ≤ 60%',
            ifelse(ccf>0.60 & ccf<=0.80, '60% < CCF ≤ 80%',
            ifelse(ccf>0.80,             '80% < CCF ≤ 100%',
            NA)))))))) %>%
        mutate(ccf=ifelse(is.na(effect),NA,ccf)) %>%
        mutate(ccf=factor(ccf, levels=c('CCF = 0%','0% < CCF ≤ 5%','5% < CCF ≤ 20%','20% < CCF ≤ 40%','40% < CCF ≤ 60%','60% < CCF ≤ 80%','80% < CCF ≤ 100%'))) %>%
        arrange(!is.na(effect), !is.na(band), n.band, sample, desc(chrom), desc(pos)) %>%
        rowwise %>%
        mutate(band=str_c('chr',chrom,': ',band)) %>%
        ungroup %>%
        mutate(band=factor(band, levels=filter(.,!is.na(effect)) %>% .$band %>% unique)) %>%
        select(-order) %>%
        mutate(pathogenic=as.factor(ifelse(pathogenic=='Pathogenic','darkgoldenrod2', NA))) %>%
        mutate(clonality=ifelse(clonality=='Clonal',clonality,NA))

    # fix sample order
    if(is.null(sample.order)) {
        muts %<>% mutate(sample=factor(sample, levels=unique(sort(muts$sample))))
    } else {
        muts %<>% mutate(sample=factor(sample, levels=sample.order))
    }

    return(muts)
}


#-----------------------
# main plotting function
#-----------------------

PlotVariants <- function(muts, output.file, color.table=NULL, clonality=FALSE, ccf=FALSE, loh=TRUE, cn=FALSE, width=20, height=20, text.size=18){

    # rename for plot output
    muts %<>% rename(Sample=sample, Gene=gene, Band=band, Effect=effect, Pathogenic=pathogenic, Clonality=clonality, CCF=ccf, CN=cn)

    # plot aesthetic definitions
    palette  <- c(
        'Truncating SNV'='#C84DDD',
        'Frameshift In-Del'='#C17200',
        'Missense SNV'='#00A5A8',
        'Inframe indel'='#E44988',
        'Splice site variant'='#008AE9',
        'Upstream, start/stop, or de novo modification'='#749000',
        'Silent'='#666666',
        'Amplification'='#333399',
        'Deletion'='#e60000',
        'CCF = 0%'='#e5e5e5',
        `0% < CCF ≤ 5%`='#c7dbee',
        `5% < CCF ≤ 20%`='#a0cae0',
        `20% < CCF ≤ 40%`='#6eaed4',
        `40% < CCF ≤ 60%`='#2772b3',
        `60% < CCF ≤ 80%`='#10539a',
        `80% < CCF ≤ 100%`='#0b3269')

    geometry <- c(`LOH`=3)
    surround <- c(`Pathogenic`='darkgoldenrod2')

    # main plot params
    if(cn!=TRUE) {
        hp <- ggplot(muts,aes(Sample,Gene))
    } else {
        hp <- ggplot(muts,aes(Sample,Band))
    }

    if(ccf==TRUE & cn!=TRUE ){  # CCF coloring
        hp <- hp +
        geom_tile(data=muts, aes(fill=CCF), colour='white') +
        scale_fill_manual(breaks=names(palette), values=palette, na.value='gray90', drop=FALSE)
    } else {  # draw tiles and color
        hp <- hp +
        geom_tile(data=muts, aes(fill=Effect), colour='white') +
        scale_fill_manual(breaks=names(palette), values=palette, na.value='gray90', drop=FALSE)
    }

    # add loh as +
    if(loh==TRUE & cn!=TRUE) {
        hp <- hp +
        geom_point(data=muts, aes(shape=loh, stroke=1.5), size=2) +
        scale_shape_manual(values=geometry, guide=guide_legend(colour = 'white'))
    }

    if(clonality==TRUE) {
        hp <- hp +
        geom_tile(data=muts %>% filter(!is.na(Effect) & !is.na(Clonality)), aes(colour=Clonality), size=1.5, fill=NA)
    }

    hp <- hp +
    # specify legend
    guides(colour='white') +
    guides(colour = guide_legend(override.aes=list(alpha = 1,fill=NA))) + 

    # tile groups
    facet_wrap(~pheno, nrow=1, scales='free_x') +
    scale_x_discrete(expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +

    # theme params
    theme(  legend.title        = element_blank(),
            panel.grid.major    = element_blank(),
            panel.grid.minor    = element_blank(),
            text                = element_text(size=18),
            axis.title.x        = element_blank(),
            axis.title.y        = element_blank(),
            axis.text.x         = element_text(angle=90, vjust=0.5, hjust=1),
            axis.text.y         = element_text(face='italic', hjust=1),
            axis.text           = element_text(size=text.size),
            axis.ticks.x        = element_blank(),
            axis.ticks.y        = element_blank(),
            legend.key          = element_rect(colour='black', fill=NULL, size=1),
            legend.key.size     = unit(1.8, 'lines'),
            legend.text         = element_text(size=text.size),
            strip.background    = element_rect(fill='white'),
            strip.text.x        = element_text(colour='white', size=text.size*1.2),
            panel.background    = element_rect(fill=NA, color=NA),
            panel.border        = element_rect(fill=NA, colour='black', size=2),
            plot.margin         = unit(c(0,0,0,0), 'in'))

    # build Grob object
    hpg <- ggplotGrob(hp)

    # count number of samples in each group
    plot.lengths <- muts %>% split(.$pheno) %>% map(~ .x$Sample %>% unique %>% length) %>% unlist

    # get the column indexcorresponding to the panels.
    panelI <- hpg$layout$l[grepl('panel', hpg$layout$name)]

    # replace the default panel widths with relative heights.
    hpg$widths <- grid:::unit.list(hpg$widths)
    hpg$widths[panelI] <- lapply(plot.lengths, unit, 'null')

    # add extra width between panels
    for(gap in 1:(length(panelI)-1)){
        hpg$widths[panelI[gap]+1]=list(unit(0.3, 'cm'))
    }

    # facet coloring
    names.grobs <- grid.ls(grid.force(hpg),print=FALSE)$name
    strips <- names.grobs[which(grepl('strip.background',names.grobs))]
    if(is.null(color.table)) {
        cols <- colorRampPalette(brewer.pal(8,'Dark2'))(length(strips))  # default coloring
    } else { 
        cols <- color.table %>% select(pheno,pheno.color) %>% unique %>% arrange(pheno) %>% .$pheno.color
    }

    # modify grob object
    for(strip in 1:length(strips)){
        hpg=editGrob(grid.force(hpg), gPath(strips[strip]), gp=gpar(fill=cols[strip]))
    }

    # draw plot
    pdf(output.file, width, height, bg='white')
        grid.draw(hpg)
    dev.off()
}


#-----------------#
#                 #
# DATA PROCESSING #
#                 #
#-----------------#

# random color selection
set.seed(color.seed)
if(random.pheno.color==TRUE){
    color.palette <- grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
    color.palette %<>% sample(.,length(.))
} else {
    color.palette <- colorRampPalette(pheno.palette)(length(subsets))  # static palette
}


# read mutations file
if(length(grep('xlsx',muts.file))){
    muts <- tryCatch({
        read.xlsx(muts.file, 'MUTATION_SUMMARY') %>% tbl_df},
        error=function(e){read.xlsx(muts.file) %>% tbl_df
    })
} else {
    muts <- read.delim(muts.file, sep='\t', stringsAsFactors=FALSE) %>% tbl_df
}


# format muts with function above
muts %<>% FormatMuts(keys=keys)


# remove silent mutations if present
if(include.silent!=TRUE){
    muts %<>% filter(effect!='silent' & effect!='Silent')
}


# build pheno bar lookup tables
color.table <- muts <- muts %>% FormatMuts(col.names='pheno')
pheno.table <- map2(subsets, names(subsets), ~ { data_frame(a=.x, sample=.x) %>% setNames(c(.y, 'sample')) %>% mutate(pheno=.y) })

for(table in 1:length(pheno.table)){
    color.table <<- left_join(color.table, pheno.table[[table]] %>% select(-pheno), by='sample', copy=TRUE)
}


# recode table values as colors
color.table %<>%
    #select(-sample, -gene, -effect, -ccf, -loh, -chrom, -pos, -start, -end, -band, -pathogenic, -clonality, -ci95.low, -pr.somatic.clonal, -ref, -alt, -purity, -cn, -fathmm, -chasm, -cancer.gene, -mut.taster, -provean, -haploinsufficient) %>%
    select(pheno) %>%
    rename(pheno.color=pheno) %>%
    mutate(pheno.color=ifelse(all(is.na(pheno.color)), 'Mutations', pheno.color)) %>%
    mutate_each(funs(as.character)) %>%
    map(~ { seed = hash(.x)
            values = .x %>% ifelse(. %in% c('.', NA), '', .)
            colors = color.palette[(1:length(unique(values)))+seed]
            colors[unique(values) == ''] <- 'grey'
            setNames(colors, unique(values))[values]
        }) %>%
    bind_cols(muts[c('gene', 'sample', 'pheno')], .) %>%
    mutate(exists=1) %>%
    select(sample,gene,exists,pheno,pheno.color,everything())  # table expansion & rearrange for column naming


# reshape melted data into wide-form for distance calculation
event.matrix <-
    color.table %>%
    select(sample,gene) %>%
    unique %>%
    mutate(exists=1) %>%
    spread(sample,exists,fill=0) %>%
    data.frame(.,row.names=1,stringsAsFactors=FALSE,check.names=FALSE) %>%
    as.matrix %>%
    t


# compute distance matrix using method specified above
dist <- DistExtra(event.matrix, dist.method)


# cluster using method specified above
hc <- dist %>% hclust(method=clust.method)


# construct dendrogram
dend <- 
    hc %>%
    as.dendrogram %>% 
    set('branches_lwd', 10)


# rotate tree to order using method specified above
if(tree.sort == 'distance') {
    dend %<>% rotate_DendSer(ser_weight=dist(x))
} else {
    dend %<>% ladderize
}


# reorder rows using distance matrix min
dend <- tree.dend <- reorder(dend, dist)


#-------------------------------
# FORMAT COPY NUMBER INFORMATION
#-------------------------------

# select threshold columns
gene.cn <-
    read.delim(cn.file, sep='\t',stringsAsFactors=FALSE) %>%
    tbl_df

if(cn.cols=='threshold') {
    gene.cn %<>% select(gene=hgnc,chrom,start,end,band,matches('threshold'))
    cn.sample.keys <- grep('threshold',colnames(gene.cn), value=TRUE) %>% list.map(strsplit(.,'_') %>% unlist %>% head(1) %>% keys[.] %>% unname) %>% unlist
    names(gene.cn)[which(names(gene.cn) %in% names(cn.sample.keys))] <- cn.sample.keys[which(names(cn.sample.keys) %in% names(gene.cn))]
} else if(cn.cols!='all') {
    gene.cn %<>% select(gene=hgnc,chrom,start,end,band,one_of(cn.cols))
    cn.sample.keys <- colnames(gene.cn) %>% list.filter(!. %in% c('gene','hgnc','chrom','start','end','band'))
} else {
    cn.sample.keys <- colnames(gene.cn) %>% list.filter(!. %in% c('gene','hgnc','chrom','start','end','band'))
}


# subset amp / del rows
cnas <-
    gene.cn %>%
    mutate(amp=rowSums(.[cn.sample.keys]==2)) %>%
    mutate(del=rowSums(.[cn.sample.keys]==-2)) %>%
    filter(amp==TRUE|del==TRUE) %>%
    mutate_each(funs(ifelse(.==1,0,.)), one_of(cn.sample.keys)) %>%
    mutate_each(funs(ifelse(.==-1,0,.)), one_of(cn.sample.keys)) %>%
    distinct_(c('chrom','band',cn.sample.names), .keep_all=TRUE) %>%
    select(-gene, -amp, -del) %>%
    gather(sample, effect, -band, -chrom, -start, -end) %>%
    filter(effect==2|effect==-2) %>%
    rowwise %>%
    mutate(sample=strsplit(sample,'_') %>% unlist %>% head(1)) %>%
    ungroup %>%
    unique %>%
    FormatMuts

# write CNA file
cnas %>%
select(sample,band,effect) %>%
spread(sample,effect) %>%
write_tsv(cn.amp.del.out)

# build color tables
cna.color.table <-
    color.table %>%
    select(sample,pheno,pheno.color) %>%
    unique %>%
    right_join(cnas) %>%
    mutate(exists=1) %>%
    select(sample, band, exists, pheno, pheno.color)


#-------------------
# ABSOLUTE CCF CALLS
#-------------------

abs.maf <-
    list.files('absolute/reviewed/SEG_MAF', pattern='_ABS_MAF.txt', full.names=TRUE) %>%
    map(~ { read.delim(.x, sep='\t',stringsAsFactors=FALSE) %>% tbl_df }) %>%
    bind_rows %>%
    separate(sample, into=c('sample','normal'), sep='_') %>%
    rename(pos=Start_position, alt.dept=alt) %>%
    FormatMuts %>%
    mutate(clonality=ifelse(pr.somatic.clonal>=0.5 | ci95.low >=0.9, 'Clonal', 'Subclonal')) %>%
    select(sample,gene,pos,ccf,clonality,purity)

abs.pos <-
    list.files('absolute/tables', pattern='.somatic.txt', full.names=TRUE) %>%
    map(~ { read.delim(.x, sep='\t',stringsAsFactors=FALSE) %>% tbl_df %>% select(sample=Sample, gene=Gene, pos=Position, alt=Alt)}) %>%
    bind_rows %>%
    separate(sample, into=c('sample','normal'), sep='_')

abs.ccf <- left_join(abs.maf, abs.pos)

#----------
# LOH CALLS
#----------

# read cnf & make calls
cncf.loh <-
    list.files('facets',pattern='.cncf.txt',full.names=TRUE) %>%
    map(~ {
        cncf <-
            read.delim(.x,stringsAsFactors=FALSE,sep="\t") %>%
            tbl_df %>%
            separate(ID, into=c('sample','normal'), sep='_') %>%
            # loh calclulation
            mutate(loh=
                ifelse(lcn.em==0,'LOH',
                ifelse(lcn.em>0,'.',
                ifelse(lcn==0,'LOH',
                ifelse(lcn>0,'.',
                ifelse(!is.na(lcn.em) & !is.na(lcn), '.',
                NA)))))) %>%
            mutate(loh=ifelse(lcn==0 & cnlr.median.clust < 0 & mafR > 0.2, 'LOH', '.')) %>%
            mutate(loc.mid=(loc.start+loc.end)/2) %>%
            mutate(id=row_number())

        cncf %<>%
            group_by(id) %>%
            full_join(
                cncf %>% filter(!is.na(loh)) %>%
                select(chrom, loc.mid.adj=loc.mid, cnlr.median.clust.adj=cnlr.median.clust, adj.lcn.em=lcn.em, adj.lcn=lcn, adj.loh=loh) %>%
                bind_rows(data_frame(chrom=1:max(cncf$chrom), loc.mid.adj=-Inf, cnlr.median.clust.adj=1, adj.lcn.em=1, adj.lcn=1, adj.loh='.')),
            by='chrom') %>%
            slice(which.min(abs(loc.mid-loc.mid.adj))) %>%
            ungroup %>%
            mutate(loh=ifelse(is.na(loh), adj.loh, loh)) %>%
            mutate(loh=ifelse(cnlr.median.clust < 0 & mafR > 0.2, 'LOH', '.')) %>%
            mutate(loh=ifelse(loh=='.',NA,loh)) %>%
            select(sample,chrom,loc.start,loc.mid,loc.end,loh,lcn.em,lcn,tcn.em,tcn,mafR,cnlr.median,cnlr.median.clust,cf,seg,num.mark,cf)
        return(cncf)
    }) %>%
    bind_rows %>%
    select(sample, chrom, loc.start, loc.mid, loc.end, loh, cf)


#-----------------------
# COMBINE MUTS, CCF, LOH
#-----------------------

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

geno.pos <-
    tcga_maf[, c('Hugo_Symbol', 'Chromosome', 'POS')] %>%
        tbl_df %>%
        select(gene=Hugo_Symbol, chrom=Chromosome, pos=POS) %>%
        FormatMuts %>%
        filter(!is.na(gene) & !is.na(chrom) & !is.na(pos)) %>%
        bind_rows(select(msk.sufam.muts,gene,chrom, pos) %>% unique) %>%
        filter(!is.na(gene) & !is.na(chrom) & !is.na(pos))%>%
        tbl_df %>%
        unique

muts <-
all.variants %>%
mutate(id=row_number()) %>%
group_by(id) %>%
full_join(geno.pos) %>%
arrange(!is.na(chrom), !is.na(pos)) %>%
top_n(1) %>%
ungroup %>%
select(sample,gene,effect,loh=LOH,chrom,pos) %>%
full_join(foo %>% select(-loh), by=c('sample','gene','effect','chrom','pos')) %>%
unique



muts %<>% mutate(sample=values[sample]) %>%
    filter(!is.na(sample))

muts <- all.variants %>% select(sample,gene,effect,loh=LOH)
muts0 <- muts

muts %<>%
    #select(-ccf,-loh,-clonality, -purity) %>%
    #select(-ccf,-clonality, -purity) %>%
    left_join(abs.ccf, by=c('sample','gene','pos'))

muts <-
muts %>%
    mutate(id=row_number()) %>%
    group_by(id) %>%
#   full_join(cncf.loh, by=c('sample','chrom')) %>%
    full_join(cncf.loh, by=c('sample','chrom')) %>%
    mutate(pos=as.integer(pos)) %>%
    mutate(in.seg=(loc.start<=pos & pos<=loc.end)) %>%
    #filter(in.seg==TRUE | sum(in.seg)==0) %>%
    #top_n(1) %>%
    #slice(which.min(abs(pos-loc.mid))) %>%
    slice(which.min(pos)) %>%
    ungroup %>%
    mutate(loh=loh.x)
    #mutate(loh=ifelse(in.seg==FALSE,NA,loh)) %>%
    #select(-in.seg)

muts %<>% mutate(ccf=ifelse(is.na(ccf),cf,ccf))

muts %<>% full_join(
            all.variants %>%
            filter(pheno.bar=='TCGA') %>%
            select(sample,gene,effect) %>%
            mutate(sample=values[sample]) %>%
            mutate(pheno='TCGA') %>%
            FormatMuts %>%
            select(sample,gene,effect,pheno) ) %>%

muts %<>% mutate(sample=keys[sample]) %>% filter(!is.na(sample))

#----------------------------
# MUTS PATHOGENIC CALCULATION
#----------------------------

muts %<>%
    mutate(pathogenic=
        ifelse(effect=='Missense SNV' & fathmm=='CANCER' & cancer.gene=='TRUE', 'Pathogenic',   # 'Potentially Pathogenic' should be added, mutation taster should be added
        ifelse(effect=='silent' & haploinsufficient == 'TRUE', 'Pathogenic',    # in-frame needs to be broken out above, haploinsufficient needs to come from somewhere
        ifelse((effect=='Truncating SNV' | effect=='Splice site variant' | effect=='Frameshift In-Del') & (loh=='LOH' | cancer.gene=='TRUE'), 'Pathogenic',
        'Passenger')))) %>%
    mutate(pathogenic=ifelse(pathogenic=='Pathogenic',pathogenic,NA))


#--------------------------
# BUILD VARIANT & CCF PLOTS
#--------------------------

muts <- all.variants %>% left_join(abs.ccf %>% mutate(sample=keys[sample]), by=c('sample','gene')) %>% filter(!is.na(ccf)) %>% FormatMuts %>% mutate(pheno='all') %>% filter(!is.na(effect))



muts %<>% filter(effect!='silent') %>% filter(effect!='Silent') %>% filter(!is.na(effect))

# variant plots
muts %>%
#OrgMuts(labels(dend)) %>%
OrgMuts(foo) %>%
mutate(pheno=ifelse(all(is.na(pheno)),'Variants by type',pheno)) %>%
PlotVariants('summary/mutation_heatmap_variant.pdf', color.table, clonality=FALSE, ccf=FALSE, loh=FALSE, cn=FALSE, width=length(unique(.$sample))+10, height=(length(unique(.$gene))/3)+10, text.size=30)

muts %>%
#OrgMuts(labels(dend), recurrence=2) %>%
OrgMuts(foo, recurrence=2) %>%
mutate(pheno=ifelse(all(is.na(pheno)),'Recurrent variants by type',pheno)) %>%
PlotVariants('summary/mutation_heatmap_variant_recurrent.pdf', color.table, clonality=FALSE, ccf=FALSE, loh=FALSE, cn=FALSE, width=length(unique(.$sample))+10, height=(length(unique(.$gene))/3)+10, text.size=30)

# ccf plots
muts %>%
#OrgMuts(labels(dend)) %>%
OrgMuts(foo) %>%
mutate(pheno=ifelse(all(is.na(pheno)),'Variant CCF',pheno)) %>%
PlotVariants('summary/mutation_heatmap_ccf.pdf', color.table, clonality=TRUE, ccf=TRUE, loh=TRUE, cn=FALSE, width=length(unique(.$sample))+10, height=(length(unique(.$gene))/3)+10, text.size=30)

muts %>%
#OrgMuts(labels(dend), recurrence=2) %>%
OrgMuts(foo, recurrence=2) %>%
mutate(pheno=ifelse(all(is.na(pheno)),'Recurrent variant CCF',pheno)) %>%
PlotVariants('summary/mutation_heatmap_ccf_recurrent.pdf', color.table, clonality=TRUE, ccf=TRUE, loh=TRUE, cn=FALSE, width=length(unique(.$sample))+10, height=(length(unique(.$gene))/3)+10, text.size=30)


# subset comparison plots
system("mkdir summary/subsets &>/dev/null")

pheno.table.strip <- pheno.table %>% map(~ {
        .x %>%
        #setNames(c('sample.id','sample','pheno')) })
        setNames(c('sample','sample.id','pheno')) })

map2(subset.groups, names(subset.groups), ~ {

    # for file naming
    subgroup.ext <- gsub(' ', '_', .y)


    #------------------
    # mutation heatmaps
    #------------------

    # subset mutation table
    message(green('[subset mutation heatmaps]'))
    muts <- 
        pheno.table.strip[unlist(.x)] %>%
        bind_rows %>%
        inner_join(muts %>% select(-pheno), .) %>%
        unique %>%
        filter(!is.na(pheno))

    # variant plots
    muts %>%
    OrgMuts(labels(dend), recurrence=0) %>%
    PlotVariants(str_c('summary/subsets/mutation_heatmap_variant.',subgroup.ext,'.pdf'), color.table, clonality=FALSE, ccf=FALSE, loh=TRUE, cn=FALSE, width=length(unique(.$sample))+10, height=(length(unique(.$gene))/3)+10, text.size=30)

    muts %>%
    OrgMuts(labels(dend), recurrence=1) %>%
    PlotVariants(str_c('summary/subsets/mutation_heatmap_variant_recurrent.',subgroup.ext,'.pdf'), color.table, clonality=FALSE, ccf=FALSE, loh=TRUE, cn=FALSE, width=length(unique(.$sample))+10, height=(length(unique(.$gene))/3)+10, text.size=30)

    # ccf plots
    muts %>%
    OrgMuts(labels(dend), recurrence=0) %>%
    PlotVariants(str_c('summary/subsets/mutation_heatmap_ccf.',subgroup.ext,'.pdf'), color.table, clonality=TRUE, ccf=TRUE, loh=FALSE, cn=FALSE, width=length(unique(.$sample))+10, height=(length(unique(.$gene))/3)+10, text.size=30)

    muts %>%
    OrgMuts(labels(dend), recurrence=1) %>%
    PlotVariants(str_c('summary/subsets/mutation_heatmap_ccf_recurrent.',subgroup.ext,'.pdf'), color.table, clonality=TRUE, ccf=TRUE, loh=FALSE, cn=FALSE, width=length(unique(.$sample))+10, height=(length(unique(.$gene))/3)+10, text.size=30)


    #------------------------------------
    # copy number Fisher's exact plotting
    #------------------------------------

    samples.a <- subsets[.x[1]] %>% unlist
    gene.cn.a <- gene.cn[c('gene', 'start', 'end', samples.a)]

    samples.b <- subsets[.x[2]] %>% unlist
    gene.cn.b <- gene.cn[c('gene', 'start', 'end', samples.b)]

    GenCN(gene.cn.a, gene.cn.b, plot.title.main=.y, plot.title.a=.x[1], plot.title.b=.x[2], allosome, targets.file=NULL, suffix='', threshold.a=FALSE, threshold.b=FALSE, gene.names=FALSE)

})


# subsets plot dir
combn.subsets <-
    expand.grid(a=names(subsets), b=names(subsets), stringsAsFactors=0) %>%
    rowwise %>%
    filter(intersect(unlist(subsets[a]), unlist(subsets[b])) %>% length == 0)


#-----------
# CLUSTERING
#-----------

# adaptive branch pruning, order by ladderized tree layout
if(color.clusters == TRUE) {
    clusters <-
        cutreeDynamic(
            hc,
            minClusterSize = min.cluster,
            distM = as.matrix(dist),
            method = 'hybrid',
            deepSplit = 4,              # clustering sensitivity [1-4]
            maxCoreScatter = NULL,      # max scatter of the core for a branch to be a cluster given as absolute heights [0-1]
            minGap = NULL,              # min cluster gap given as fraction of the difference between ‘cutHeight’ and the 5th percentile of joining heights [0-1]
            maxAbsCoreScatter = NULL,   # max scatter of the core for a branch to be a cluster given as absolute heights
            minAbsGap = NULL            # min cluster gap given as absolute height difference
        ) %>%
        .[order.dendrogram(tree.dend)]

    # cluster coloring
    cluster.v <- unique(clusters) %>% list.filter(.!=0)
    n.cluster <- length(cluster.v)

    #cluster palette
    cluster.palette <- colorRampPalette(brewer.pal(8,'Dark2'))(n.cluster)

    # color branches according to cluster
    tree.dend %<>% branches_attr_by_clusters(clusters, values=cluster.palette)
}


# remove labels from tree
if(tree.labels == FALSE){
    dend.tree %<>% set('labels', '')
}


#----------------
# TREE GENERATION
#----------------

# distance tree
pdf('summary/genomic_distance_tree_ladder.pdf',140,38)
    par(mar = c(10,10,10,10))  # bottom, left, top, right
    layout(matrix(c(1,2),nrow=1), widths=c(10,1))
    plot(dend, cex.axis=6)
    colored_bars( colors = color.table %>% select(-sample, -gene, -exists, -pheno),
                  dend = dend,
                  sort_by_labels_order = FALSE,
                  add = TRUE,
                  rowLabels = color.table$sample,
                  y_scale = 10,
                  cex.rowLabels = 5.8 )
    color.table %>%
    select(-sample, -gene, -exists) %>%
    map(~ accumulate(., ~ legend(x=3, y=3, legend=unique(names(.x)), fill=unique(.x), cex=4) ))

dev.off()

# weighted tree (charlotte's)
pdf('summary/genomic_distance_tree_weighted.pdf',60,25)
    heatmap.2(
        event.matrix,
        trace='none',
        Rowv=FALSE,
        hclustfun=function(x){hclust(x, 'ward.D2')},
        col=c('white','black'),
        reorderfun=function(d,w) rev(reorder(d,w)),
        distfun=function(x) as.dist(Hamming(x)),
        key=FALSE
            )
dev.off()


pdf('summary/genomic_distance_tree_bw.pdf',60,25)
    heatmap.2(
        t(event.matrix),
        trace='none',
        dendrogram='column',
        Colv=rev(dend),
        col=c('white','black'),
        key=FALSE
    )
dev.off()


#-----------------------
# LINEAGE MAP GENERATION
#-----------------------

plot_file          = 'summary/tree.pdf'
#colnames(muts)[1]  = 'Sample.ID'
colnames(muts)[7]  = 'Gene'
colnames(muts)[8]  = 'AA'
colnames(muts)[10] = 'Effect'
colnames(muts)[26] = 'CCF'
muts               = muts[!is.na(muts$ccf),]
tumor              = unique(muts$sample)


muts %<>% as.data.frame

sample_names   = as.list(sort(unique(muts$sample)))
mutation_genes = unique(muts$Gene)
rownames(muts) = 1:nrow(muts)
TCGA=FALSE


# Make a matrix of "blank" values the with nrow= #mutations and ncol=#samples
mutation_heatmap <- matrix(0, nrow=sum(unlist(lapply(sample_names, length))), ncol=sum(unlist(lapply(mutation_genes, length))))
rownames(mutation_heatmap) <- unlist(sample_names)
colnames(mutation_heatmap) <- mutation_genes


# Make sure the sample and mutations are both in the list of gene mutations and gene samples
if (!TCGA) { smallmaf <- muts[which(muts$Gene %in% unlist(mutation_genes) & muts$Sample.ID %in% unlist(sample_names)),]
}else {
        muts$id <- unlist(lapply(muts$Sample.ID, function(x){substr(x, 1, 12)}))
        print(head(muts$id))
        print(head(unlist(sample_names)))
        smallmaf <- muts[which(muts$Hugo_Symbol %in% unlist(mutation_genes) & muts$id %in% unlist(sample_names)),]
}

# for each row read the Effect and create the type based on which category it fits in
for (i in 1:nrow(smallmaf)) {
        if(!TCGA) { type = smallmaf$CCF[i] } else { type = smallmaf$Variant_Classification[i] }
        print(paste(i,type,sep="_"))

        if (!TCGA) {
        if (mutation_heatmap[which(rownames(mutation_heatmap)==smallmaf$Sample.ID[i],), which(colnames(mutation_heatmap)==smallmaf$Gene[i])] < type) {
            mutation_heatmap[which(rownames(mutation_heatmap)==smallmaf$Sample.ID[i],), which(colnames(mutation_heatmap)==smallmaf$Gene[i])] <- type}
        } else { mutation_heatmap[which(rownames(mutation_heatmap)==smallmaf$id[i]), which(colnames(mutation_heatmap)==smallmaf$Hugo[i])] <- type }
}

smalltab = t(mutation_heatmap)
smalltab = smalltab>0

smalltab = cbind(smalltab, FALSE)
colnames(smalltab)[ncol(smalltab)] = "Parental"

# smalltab = cbind(smalltab, FALSE)
# colnames(smalltab)[ncol(smalltab)] = "Test"

pd <- phyDat(t(smalltab), type="USER", levels=c(FALSE, TRUE))
dm <- dist.hagene.cning(pd)
tree <- njs(dm)
treeRatchet <- pratchet(pd, start=tree)
treeRatchet <- acctran(treeRatchet, pd)
lineage.map <- root(treeRatchet, "Parental")

pdf("summary/lineage_map.pdf")
    plot(lineage.map)
dev.off()



