#!/usr/bin/env Rscript

# script for build mutation heatmaps, trees & other standard plots

#-----------#
#           #
# LIBRARIES #
#           #
#-----------#

pacman::p_load( dplyr,lazyeval,readr,tidyr,magrittr,purrr,stringr,rlist,openxlsx,  # base
                crayon,colorspace,RColorBrewer,  # coloring
                ggplot2,grid,gridExtra,gplots,  # plot layout
                dendextend,dendextendRcpp,dynamicTreeCut,gclus,phangorn,  # dentrogram
                gtools,  # permutation
                optparse )  # option parsing

loadNamespace('plyr')

# set graphics device
options(device=pdf)  # options(encoding='ISOLatin2.enc') | pdf.options(encoding='ISOLatin2.enc')


#-------------------#
#                   #
# OPTIONS / PARSING #
#                   #
#-------------------#

if(!interactive()) {

    # define input options
    opts.list <- list( # input
                       make_option('--project_config',      default='', help='project configuration file'),
                       make_option('--samples_config',      default='', help='sample configuration file'),
                       make_option('--mutation_table_in',   default='', help='mutation summary file'),
                       make_option('--gene_cn_in',          default='', help='per-gene copy number file'),
                       make_option('--cncf_path',           default='', help='location of facets cncf files'),
                       make_option('--seg_maf_path',        default='', help='location of absolute segmented mafs'),
                       # run options
                       make_option('--include_silent',      default='', help='should silent variants be plotted'),
                       # output
                       make_option('--muts_table_out',      default='', help='mutation output summary'),
                       make_option('--cnas_table_out',      default='', help='copy number aberration output summary'),
                       make_option('--variants_table_out',  default='', help='all-variant output summary') )

    # parse input
    parser <- OptionParser(usage="%prog [options] [project_config] [samples_config] [mutation_table_in] [gene_cn_in] [cncf_path] [seg_maf_path] [muts_table_out] [cnas_table_out] [variants_table_out]", option_list=options.list)

    # build options table
    opts <-
        parse_args(parser, positional_arguments=TRUE)$options %>%
        head(-1) %>%
        stack %>%
        mutate(ind=as.character(ind)) %>%
        rowwise %>%
        mutate(pass=(if(values == ''){      # throw error if argument is missing
                print_help(parser)
                stop(str_c('missing ', ind, ' file'))
            } else {TRUE}))

    opts <- opts$values %>% set_names(opts$ind) %>% as.list

} else {
    # set detaults for interactive use
    opts = list( # input
                 project_config         = 'project_config.yaml',
                 samples_config         = 'samples.yaml',
                 muts_table_in          = 'summary/mutation_summary.xlsx',
                 gene_cn_in             = 'facets/geneCN.fill.txt',
                 seg_maf_path           = 'absolute/reviewed/SEG_MAF',
                 cncf_path              = 'facets/cncf',
                 # run options
                 include_silent         = 'True',
                 # output
                 muts_table_out         = 'summary/mutations.tsv',
                 cnas_table_out         = 'summary/copy_number_alterations.tsv',
                 variants_table_out     = 'summary/variants.tsv' )
}

# make subdirectories if absent
system("mkdir summary/trees &>/dev/null")
system("mkdir summary/sub_groups &>/dev/null")

#------------
# run options
#------------

include.silent       = if(opts$include_silent == 'true') { TRUE } else { FALSE }
dist.method          = 'hamming'   # 'hamming', euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'    | see ?dist for details [hamming method implemented manually]
clust.method         = 'complete' # 'complete', ward', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid' | see ?hclust for details
exclude.values       = c('', '.', 'Normal', 'Not Performed', 'Performed but Not Available', 'FALSE', FALSE)
sort.method          = 'distance'  # for tree sorting: ladder 'ladderize'
tree.labels          = TRUE
pheno                = NULL
color.seed           = 3
color.clusters       = TRUE
min.cluster          = 3
random.pheno.color   = TRUE
pheno.palette        = c('#2d4be0', '#20e6ae', '#ccb625', '#969696')
cn.cols              = 'threshold'  # 'threshold' = reserved string for selecting columns with threshold sufix, 'all' = reserved string for using all columns, else a string specifing columns to use
use.keys             = FALSE
loh.closest          = TRUE  # should copy number / loh assignments be made according to closest segment if variant does not fall within segment: bool
call.loh             = TRUE
call.abs             = TRUE
gene.names           = FALSE


#-----------#
#           #
# FUNCTIONS #
#           #
#-----------#

# Fishers functions
source('modules/summary/variantFishers.R')

#--------------------------
# Hamming distance function
#--------------------------

Hamming <- function(event.matrix) {
    D <- (1 - event.matrix) %*% t(event.matrix)
    D + t(D)
}


#---------------------------
# add dummy column if absent
#---------------------------

DummyCols <- function(events, col.names) {

    for(col.name in col.names) {
        if(!col.name %in% colnames(events)) {
            message(green(str_c('adding dummy column: ', col.name)))
            events.names <- colnames(events) %>% replace(which(is.na(.)),'empty') %>% list.filter(.!='add.col')
            events$add.col <- NA
            events %<>% setNames(c(events.names,col.name))
        }
    }
    return(events)

}


#--------------------------
# convert chromosome format
#--------------------------

ChromMod <- function(events, allosome){

    if('X' %in% events$chrom) { message(yellow('X chromosome labelling converted to numeric')) }

    if(allosome=='distinct') {
            if('Y' %in% events$chrom) { message(yellow('Y chromosome labelling converted to 24-based numeric')) }
            events %>%
            mutate(.,chrom=as.numeric(ifelse(chrom=='X', 23,chrom))) %>%
            mutate(chrom=as.numeric(ifelse(chrom=='Y', 24,chrom)))
    } else if(allosome =='merge') {
            if('Y' %in% events$chrom) { message(yellow('Y chromosome labelling converted to 23-based numeric')) }
            events %>% 
            mutate(chrom=as.numeric(ifelse(chrom %in% c('X','Y'), 23, chrom)))
    } else {
            events %>%
            filter(!chrom %in% c('X','Y'))
    }
}


#-----------------------------------
# soft sample key:value modification
#-----------------------------------

KeyMod <- function(sample.object, keys, force.keys=FALSE) {

    if(is.vector(sample.object)) {
        keys.index <- which(sample.object %in% names(keys))
        if(length(keys.index) == 0 & force.keys==FALSE) {
            message(green('no keys to convert'))
        } else if(force.keys==FALSE) {
            conversions <- data.frame(pre=as.character(sample.object[keys.index]), post=keys[as.character(sample.object[keys.index])], row.names=NULL) %>% unique
            message(green(str_c('key conversions: ',length(conversions$pre),'/',length(unique(sample.object)))))
            print(conversions)
            sample.object <- keys[sample.object]
        } else {
            conversions <- data.frame(pre=as.character(unique(sample.object)), post=keys[unique(as.character(sample.object))], in.keys=unique(as.character(events)) %in% names(keys), row.names=NULL)
            message(yellow('forcing key conversions'))
            print(conversions)
            sample.object[keys.index] <- keys[sample.object[keys.index]]
        }
    } else {
        keys.index <- which(sample.object$sample %in% names(keys))
        if(length(keys.index) == 0 & force.keys==FALSE) {
            message(green('no keys to convert'))
        } else if(force.keys==FALSE) {
            conversions <- data.frame(pre=as.character(sample.object$sample[keys.index]), post=keys[as.character(sample.object$sample[keys.index])], row.names=NULL) %>% unique
            message(green(str_c('key conversions: ',length(conversions$pre),'/',length(unique(sample.object$sample)))))
            print(conversions)
            sample.object %<>% mutate(sample=ifelse(row_number() %in% keys.index, keys[sample], sample))
        } else {
            conversions <- data.frame(pre=as.character(unique(sample.object$sample)), post=keys[unique(as.character(sample.object$sample))], in.keys=unique(as.character(events$sample)) %in% names(keys), row.names=NULL)
            message(yellow('forcing key conversions'))
            print(conversions)
            sample.object %<>% mutate(sample=keys[sample])
        }
    }
    return(sample.object)

}


#-----------------------
# coerce columns to type
#-----------------------

TypeCol <- function(events, columns, types) {
    for (col in 1:length(columns)) {
        if(columns[col] %in% colnames(events)) {
            if(types[col] == 'char') {
                events[,columns[col]] <- events[,columns[col]] %>% unlist %>% as.character
            } else if(types[col] == 'num') {
                events[,columns[col]] <- events[,columns[col]] %>% unlist %>%  ifelse(.=='.',NA,.) %>% as.numeric
            } else if(types[col] == 'logic') {
                events[,columns[col]] <- events[,columns[col]] == 'TRUE' %>% as.vector
            }
        }
    }
    return(events)
}


#-----------------------
# melted table to matrix
#-----------------------

# function to convert event table into sample-gene occurance matrix
MeltToMatrix <- function(event.melt) {
    # reshape melted events to wide matrix for distance calculation
    event.melt %>%
    select(sample, span) %>%
    unique %>%
    mutate(exists=1) %>%
    spread(sample,exists, fill=0) %>%
    data.frame(., row.names=1, stringsAsFactors=FALSE, check.names=FALSE) %>%
    as.matrix %>%
    t
}


#------------------------
# distance method wrapper
#------------------------

# wrap Hamming function for easy specification
DistExtra <- function(event.matrix, dist.method){
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


#---------------
# fix file names
#---------------

FixName <- function(file.name) {
    file.name %>% gsub('\\s', '_', .) %>% gsub('/', '-', .)
}


#-------------------------------------
# melted event table into ordered tree
#-------------------------------------

MeltToTree <- function(events, dist.method='hamming', clust.method='complete', sort.method='distance', span='gene', tree.samples) {

    events %<>% rename_(span=span)

    # build master event matrix
    event.matrix <- events %>% MeltToMatrix

    missing.samples <- tree.samples[!tree.samples %in% row.names(event.matrix)]

    if(length(missing.samples) > 0) {
        message(yellow('found samples with zero events, adding common root'))
        missing.matrix <- matrix(0, ncol=length(missing.samples), nrow=ncol(event.matrix)) %>% as_data_frame %>% set_names(missing.samples) %>% as.matrix %>% t
        event.matrix %<>% rbind(missing.matrix)
        event.matrix %<>% cbind( data_frame(project=rep(1,length(tree.samples))) %>% as.matrix)
    }

    # compute distance matrix using method specified
    dist <- DistExtra(event.matrix, dist.method)

    # cluster using method specified
    hc <- dist %>% hclust(method=clust.method)

    # construct dendrogram
    dend <- 
        hc %>%
        as.dendrogram %>% 
        set('branches_lwd', 10)

    # rotate tree, order using method specified
    if(sort.method == 'distance') {
        dend %<>% rotate_DendSer(ser_weight=dist(x))
    } else {
        dend %<>% ladderize
    }

    # reorder rows using distance matrix min
    dist <- reorder(dend, dist)

    tree <- list( event.matrix=event.matrix,
                 dist=dist,
                 hc=hc,
                 dend=dend,
                 dist=dist )

    return(tree)
}


#---------------------------------
# multiple pheno columns to single
#---------------------------------

MergePheno <- function(events, merge.cols, col.name='pheno') {
    events %>%
    do({ df <- .
        pheno <- select(df, one_of(merge.cols)) %>%
        apply(1, function(pheno){
            pheno %<>% na.omit
            if(length(pheno)==0){
                pheno=NA
            } else {
                pheno
            }
        })
        df <- mutate(df, col.name=pheno)
        colnames(df)[colnames(df) == 'col.name'] <- col.name
        return(df)
    })
}

#------------------------------------
# rename columns & clean nomenclature
#------------------------------------

FormatEvents <- function(events, col.names=NULL, drop=FALSE, allosome='merge', keys=NULL, force.keys=FALSE) {

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
                   'cancer.gene'='cancer.gene', 'Cancer.Gene.Census'='cancer.gene', 'Cancer Gene Census'='cancer.gene','Cancer5000.S.genes..Lawrence.et.al.'='cancer.gene', 'cancer_gene_census'='cancer.gene', 'cancer.gene'='cancer.gene',
                   'ccf'='ccf', 'CCF'='ccf', 'cancer_cell_frac'='ccf', 'Cancer cell fraction'='ccf', 'Cancer.cell.fraction'='ccf',
                   'chasm'='chasm', 'CHASM'='chasm', 'Breast_chasm_score'='chasm',
                   'chrom'='chrom','Chrom'='chrom','Chromosome'='chrom','CHROM'='chrom',
                   'Clonal'='clonal', 'clonal'='clonal',
                   'cn'='cn', 'CN'='cn',
                   'ci95.low'='ci95.low', 'ccf_CI95_low'='ci95.low',
                   'NORMAL.DP'='depth.n', 'depth.n'='depth.n',
                   'TUMOR.DP'='depth.t', 'depth.t'='depth.t',
                   'effect'='effect','Hugo_Symbol'=='gene', 'Effect'='effect', 'Variant_Classification'='effect', 'ANN....EFFECT'='effect', 'ANN[*].EFFECT'='effect',
                   'ExAC_AF'='ex.af', 'ex.af'='ex.af',
                   'end'='end', 'end'='stop',
                   'fathmm'='fathmm', 'FATHMM'='fathmm', 'fathmm_pred'='fathmm',
                   'gene'='gene', 'Gene'='gene', 'Hugo_Symbol'='gene', 'GENE'='gene', 'hgnc'='gene', 'ANN....GENE'='gene', 'ANN[*].GENE'='gene', 'Gene.symbol'='gene',
                   'haploinsufficient'='haplo', 'hap_insuf'='haplo', 'haplo'='haplo',
                   'ANN[*].HGVS_C'='hgvs.c',
                   'ANN[*].HGVS_P'='hgvs.p',
                   'ANN[*].IMPACT'='impact', 'impact'='impact',
                   'kandoth'='kandoth', '127 significantly mutated genes (Kandoth et al)'='kandoth', '127 significantly.mutated genes.(Kandoth et al)'='kandoth','X127.significantly.mutated.genes..Kandoth.et.al.'='kandoth',
                   'lawrence'='lawrence', 'Cancer5000-S genes (Lawrence et al)'='lawrence', 'Cancer5000-S genes.(Lawrence et al)'='lawrence',
                   'loh'='loh', 'LOH'='loh', 'Loss.of.heterozygocity.(LOH)'='loh', 'Loss of heterozygocity (LOH)'='loh','Loss.of.heterozygocity..LOH.'='loh',
                   'NORMAL_SAMPLE'='normal', 'normal'='normal',
                   'NORMAL_MAF'='maf.n', 'maf.n'='maf.n',
                   'maf.t'='maf.t', 'Mutant allele fraction', 'maf.t', 'Mutant.allele.fraction'='maf.t', 'TUMOR_MAF'='maf.t',
                   'mut.taster'='mut.taster', 'dbNSFP_MutationTaster_pred'='mut.taster',
                   'pathogenic'='pathogenic', 'Pathogenic'='pathogenic', 'pathogenicity'='pathogenic', 'Pathogenicity'='pathogenic',
                   'pheno'='pheno', 'pheno.bar'='pheno',
                   'pos'='pos','POS'='pos', 'Position'='pos', 'position'='position',
                   'pr.somatic.clonal'='pr.somatic.clonal', 'Pr_somatic_clonal'='pr.somatic.clonal',
                   'provean'='provean', 'Provean'='provean',
                   'purity'='purity',
                   'ref'='ref', 'REF'='ref', 'Reference.allele'='ref',
                   'sample'='sample', 'Sample'='sample', 'Tumor_Sample_Barcode'='sample', 'TUMOR_SAMPLE'='sample', 'Sample.ID'='sample',
                   'start'='start', 'Start'='start', 'Start_position'='start' )

    # save column names to fill unhandled
    fallback.cols <- colnames(events)

    # rename columns
    names(events) <- col.keys[names(events)]

    # specify output columns
    if(is.null(col.names)) {
        col.names <- colnames(events)
        cols.match <- col.names %>% list.filter(!is.na(.))
        col.names[which(is.na(col.names))] <- fallback.cols[which(is.na(col.names))]
        names(events) <- col.names
    } else if('all' %in% col.names) {
        col.names <- c( 'alt','band','cancer.gene','ccf','chasm','chrom','clonal','cn','ci95.low','effect',
                        'end','fathmm','gene','haploinsufficient','kandoth','lawrence','loh','maf','mut.taster',
                        'pathogenic','pheno','pos','pr.somatic.clonal','provean','purity','ref','sample','start' )
    }

    # df typing to avoid row name conflicts
    events %<>% as.data.frame

    # add columns if absent
    events %<>% DummyCols(col.names)

    # select columns
    events <- events[colnames(events) %>% list.filter(.!='empty')]

    # chromosome type conversion
    events %<>% ChromMod(allosome)

    # column typing
    events %<>% TypeCol( c('sample', 'gene', 'effect', 'ref',  'alt',  'cancer.gene', 'kandoth', 'lawrence', 'haplo', 'fathmm', 'chasm'),
                       c('char',   'char', 'char',   'char', 'char', 'logic',       'logic',   'logic',    'logic', 'char',   'num') )

    if('effect' %in% colnames(events)){
        # rename variant classifications
        events %<>%
            tbl_df %>%
            { events <- .
                if(!all(is.na(events$effect))){ filter(events,!is.na(effect)) }
                events
            } %>%
            unique %>%
            mutate(effect=
                ifelse(effect%in%c('STOP_GAINED','Nonsense_Mutation','stop_gained&splice_region_variant','stop_gained','Nonsense_Mutation','Stop_Codon_Ins','nonsense','truncating snv','Truncating snv','Truncating snv','Truncating SNV'),'Truncating SNV',
                ifelse(effect%in%c('FRAME_SHIFT','FRAME_SHIFT','Frame_Shift_Del','Frame_Shift_Ins','frameshift_variant&splice_acceptor_variant&splice_region_variant&intron_variant','frameshift_variant','frameshift_variant&stop_gained','frameshift_variant&splice_region_variant','frameshift_variant&splice_acceptor_variant&splice_region_variant&splice_region_variant&intron_variant','Frame_Shift_Del','Frame_Shift_Ins','frame_shift_del','frame_shift_ins','frameshift indel','Frameshift indel','Frameshift In-Del','frameshift_variant'),'Frameshift In-Del',
                ifelse(effect%in%c('NON_SYNONYMOUS_CODING','STOP_LOST','Missense_Mutation','missense_variant','missense_variant&splice_region_variant','missense_variant|missense_variant','Missense_Mutation','missense','missense snv','Missense snv','Missense SNV'),'Missense SNV',
                ifelse(effect%in%c('CODON_CHANGE_PLUS_CODON_DELETION','CODON_DELETION','CODON_INSERTION','In_Frame_Ins','In_Frame_Del','disruptive_inframe_deletion','disruptive_inframe_insertion','inframe_deletion','inframe_insertion','disruptive_inframe_deletion&splice_region_variant','inframe_deletion&splice_region_variant','In_Frame_Del','In_Frame_Ins','in_frame_del','in_frame_ins','inframe indel','Inframe indel','Inframe In-Del'),'Inframe In-Del',
                ifelse(effect%in%c('SPLICE_SITE_DONOR','SPLICE_SITE_ACCEPTOR','SPLICE_SITE_REGION','Splice_Site','splice_donor_variant&intron_variant','splice_acceptor_variant&intron_variant','splicing','splice_donor_variant&splice_region_variant&intron_variant','splice_donor_variant&disruptive_inframe_deletion&splice_region_variant&splice_region_variant&intron_variant','Splice_Site','splice','splice site variant','Splice site variant','missense_variant & splice_region_variant'),'Splice site variant',
                ifelse(effect%in%c('STOP_LOST','START_LOST','START_GAINED','UTR_5_PRIME','start_lost','stop_lost',"5'UTR","5'Flank",'De_novo_Start_InFrame','De_novo_Start_OutOfFrame','Stop_Codon_Del','Start_Codon_SNP','Start_Codon_Ins','Start_Codon_Del','Nonstop_Mutation','nonstop','upstream, start/stop, or de novo modification','Upstream, start/stop, or de novo modification'),'Upstream, start/stop, or de novo modification',
                ifelse(effect%in%c('synonymous_variant','splice_region_variant&synonymous_variant','splice_region_variant&synonymous_variant','non_coding_exon_variant','upstream_gene_variant','downstream_gene_variant','intron_variant','frameshift_variant&splice_donor_variant&splice_region_variant&splice_region_variant&intron_variant','non_coding_exon_variant|synonymous_variant','SYNONYMOUS_CODING','synonymous_variant|synonymous_variant','splice_region_variant&synonymous_variant|splice_region_variant&non_coding_exon_variant','splice_acceptor_variant & intron_variant','intragenic_variant',"3'UTR",'IGR','lincRNA','RNA','Intron','silent','intron_exon','silent','Silent','intron_variant & missense_variant'),'Silent',
                ifelse(effect%in%c('Amplification','amplification','amp','2'),'Amplification',
                ifelse(effect%in%c('Deletion','deletion','del','-2'),'Deletion',
                ifelse(is.na(effect),NA,str_c('UNACCOUNTED: ', effect))))))))))))

        # warn on unknown effects
        if(all(is.na(events$effect))){
            message(yellow('warning: empty effect column, left as NA'))
        } else if(events %>% filter(is.na(effect)) %>% nrow > 0) {
            message(yellow('warning: some variants not accounted for'))
        }
    }

    # remove unmutated LOH if present & format
    if('loh' %in% colnames(events)) {
        events %<>% mutate(loh=ifelse(loh%in%c('LOH','loh','Loss of heterozygosity') & !is.na(effect),'Loss of heterozygosity',NA))
    }

    if(drop!=FALSE) {
        events %<>% select(one_of(cols.match))
    }

    if(!is.null(keys)) {
        events %<>% KeyMod(keys)
    }

    return(events)
}


#----------------------------------
# prepare melted table for plotting
#----------------------------------

OrgEvents <- function(events, sample.order, pheno.order, sub.sets.pheno, recurrence=1, allosome='merge', event.type='gene') {

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

    missing.samples <- sample.order[!sample.order %in% events$sample]

    missing.pheno <-
        sub.sets.pheno %>%
        MergePheno(merge.cols=pheno.order, col.name='pheno') %>%
        select(sample, pheno) %>%
        filter(sample %in% missing.samples)

    if(length(missing.samples) > 0){

        message(yellow(str_c('missing specified samples: ',missing.samples,'\n')))

        missing.fill <-
            expand.grid(missing.samples,unique(unlist(events[,event.type])), stringsAsFactors=FALSE) %>%
            set_names(c('sample', event.type)) %>%
            tbl_df %>%
            left_join(missing.pheno)
    }

    # specify output columns
    col.names <- c('sample','gene','chrom','effect','pheno','band','span','pos','maf','ccf','loh','clonal','pathogenic')

    events %<>% filter(!is.na(pheno)) %>% filter(!is.na(effect) | !is.na(ccf))

    # add dummy column names
    events %<>% DummyCols(col.names)

    # effect prescedence
    events %<>%
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
            NA))))))))))

    if(event.type=='gene') {
        events %<>%
            # remove genes with lower prescedence
            group_by(sample,gene) %>%
            arrange(gene,precedence) %>%
            slice(rank(precedence, ties.method="first")==1) %>%
            # count number of variants per gene
            unique %>%
            ungroup %>%
            group_by(gene) %>%
            mutate(n.gene=n()) %>%
            ungroup %>%
            { events <- .
                if(recurrence>0) {events %<>% filter(n.gene>=recurrence)}
                return(events)
            } %>%
            unique %>%
            arrange(desc(n.gene),sample,desc(ccf),precedence,gene)
    }

    if(event.type=='band') {
        events %<>%
            group_by(band) %>%
            mutate(n.band=n()) %>%
            ungroup %>%
            rowwise %>%
            mutate(band=str_c('chr',chrom,': ',band))
    }

    if(event.type=='span') {
        events %<>%
            # remove spans with lower prescedence
            group_by(sample,span) %>%
            arrange(span,precedence) %>%
            slice(rank(precedence, ties.method="first")==1) %>%
            # count number of variants per gene
            unique %>%
            ungroup %>%
            group_by(span) %>%
            mutate(n.span=n()) %>%
            ungroup %>%
            { events <- .
                if(recurrence>0) {events %<>% filter(n.span>=recurrence)}
                return(events)
            } %>%
            unique %>%
            arrange(desc(n.span),sample,desc(ccf),precedence,span)
    }

    # define plot gene order
    events %<>%
        mutate(order=row_number()) %>%
        ungroup %>%
        select(-precedence) %>%
        unique

    if(event.type=='gene') {
        events.fill <-
            events %>%
            select(sample,effect,pheno,gene) %>%
            unique %>%
            group_by(pheno) %>%
            spread(gene,effect) %>%
            gather(gene,effect,-pheno,-sample) %>%
            filter(is.na(effect)) %>%
            ungroup
    } else if(event.type=='band') {
        events.fill <-
            events %>%
            select(sample,chrom,effect,pheno,band) %>%
            unique %>%
            group_by(pheno) %>%
            spread(band,effect) %>%
            gather(band,effect,-chrom,-pheno,-sample) %>%
            filter(is.na(effect)) %>%
            ungroup
    } else {
        events.fill <-
            events %>%
            select(sample,chrom,effect,pheno,span) %>%
            unique %>%
            group_by(pheno) %>%
            spread(span,effect) %>%
            gather(span,effect,-chrom,-pheno,-sample) %>%
            filter(is.na(effect)) %>%
            ungroup
    }

    if(length(missing.samples) > 0){
        events %<>% bind_rows(missing.fill)
    }

    # plot aesthetics
    events %<>%
        full_join(events.fill) %>%
        filter(!is.na(sample)) %>%
        filter_(paste('!is.na(', event.type, ')' )) %>%
        filter(sample %in% sample.order) %>%
        filter(!is.na(pheno)) %>%
        # push NAs to bottom of stack
        mutate(order=ifelse(is.na(effect),-Inf,order)) %>%
        unique %>%
        # fix plot ordering & assign gene factor levels
        arrange(!is.na(effect),desc(order)) %>%
        mutate(gene=factor(gene,levels=filter(.,!is.na(effect)) %>% .$gene %>% unique)) %>%
        mutate(span=factor(span,levels=filter(.,!is.na(effect)) %>% .$span %>% unique)) %>%
        # ccf binning
        mutate(ccf=
            ifelse(ccf==0,               'CCF = 0%',
            ifelse(ccf>0.00 & ccf<=0.05, '0% < CCF <= 5%',
            ifelse(ccf>0.05 & ccf<=0.20, '5% < CCF <= 20%',
            ifelse(ccf>0.20 & ccf<=0.40, '20% < CCF <= 40%',
            ifelse(ccf>0.40 & ccf<=0.60, '40% < CCF <= 60%',
            ifelse(ccf>0.60 & ccf<=0.80, '60% < CCF <= 80%',
            ifelse(ccf>0.80,             '80% < CCF <= 100%',
            NA)))))))) %>%
        mutate(ccf=ifelse(is.na(effect),NA,ccf)) %>%
        mutate(ccf=factor(ccf, levels=c('CCF = 0%','0% < CCF <= 5%','5% < CCF <= 20%','20% < CCF <= 40%','40% < CCF <= 60%','60% < CCF <= 80%','80% < CCF <= 100%'))) %>%
        mutate(pathogenic=as.factor(ifelse(pathogenic=='Pathogenic', 'darkgoldenrod2', NA))) %>%
        mutate(clonal=ifelse(clonal=='Clonal',clonal,NA)) %>%
        select(-order) %>%
        mutate(sample=factor(sample, levels=sample.order))

    if(event.type=='gene') {
        events %<>% filter(!is.na(gene))
    } else if(event.type=='band') {
        events %<>%
            arrange(!is.na(effect), desc(n.band), desc(chrom), desc(band)) %>%
            mutate(band=factor(band, levels=filter(.,!is.na(effect)) %>% .$band %>% unique)) %>%
            filter(!is.na(band))
    } else if( event.type=='span') {
        events %<>% filter(!is.na(span))
    }

    return(events)
}


#-----------------------
# main plotting function
#-----------------------

PlotVariants <- function(events, output.file, clonal=FALSE, pathogenic=FALSE, ccf=FALSE, loh=TRUE, width=20, height=20, text.size=18, event.type='gene'){

    # rename for plot output
    events %<>% plyr::rename(replace=c(sample='Sample', gene='Gene', band='Band', span='Span', variant='Variant', effect='Effect', pathogenic='Pathogenic', clonal='clonal', ccf='CCF', cn='CN'), warn_duplicated=FALSE)

    # plot aesthetic definitions
    palette  <- c( 'Truncating SNV'='#C84DDD',
                    'Frameshift In-Del'='#C17200',
                    'Missense SNV'='#00A5A8',
                    'Inframe In-Del'='#E44988',
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
                    `80% < CCF ≤ 100%`='#0b3269' )

    geometry <- c(`LOH`=3)
    surround.clonal <- c(Clonal='darkgoldenrod2')
    surround.pathogenic <- c(Pathogenic='purple')

    # main plot params
    if(event.type=='gene') {
        hp <- ggplot(events, aes(Sample, Gene, drop=FALSE))
    } else if(event.type=='band') {
        hp <- ggplot(events, aes(Sample, Band, drop=FALSE))
    } else {
        hp <- ggplot(events, aes(Sample, Span, drop=FALSE))
    }

    if(ccf==TRUE){  # CCF coloring
        hp <- hp +
        geom_tile(data=events, aes(fill=CCF, drop=FALSE), colour='white') +
        scale_fill_manual(breaks=names(palette), values=palette, na.value='gray90', drop=FALSE)
    } else {  # draw tiles and color
        hp <- hp +
        geom_tile(data=events, aes(fill=Effect, drop=FALSE), colour='white') +
        scale_fill_manual(breaks=names(palette), values=palette, na.value='gray90', drop=FALSE)
    }

    # add loh as +
    if(loh==TRUE) {
        hp <- hp +
        geom_point(data=events, aes(shape=loh, stroke=1.5), size=2) +
        scale_shape_manual(values=geometry, guide=guide_legend(colour = 'white'))
    }

    if(clonal==TRUE) {
        hp <- hp +
        geom_tile(data=events %>% filter(!is.na(Effect) & !is.na(clonal)), aes(colour=clonal), size=1.5, fill=NA)
    }

    if(pathogenic==TRUE) {
        hp <- hp +
        geom_tile(data=events %>% filter(!is.na(Effect) & !is.na(clonal)), aes(colour=clonal), size=1.5, fill=NA)
    }

    hp <- hp +
    # specify legend
    guides(colour='white') +
    guides(colour = guide_legend(override.aes=list(alpha = 1,fill=NA))) + 

    # tile groups
    facet_wrap(~pheno, nrow=1, scales='free_x') +
    scale_x_discrete(expand=c(0, 0), drop=FALSE) + scale_y_discrete(expand=c(0, 0)) +

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
    plot.lengths <- events %>% split(.$pheno) %>% map(~ .x$Sample %>% unique %>% length) %>% unlist

    # get the column indexcorresponding to the panels.
    panelI <- hpg$layout$l[grepl('panel', hpg$layout$name)]

    # replace the default panel widths with relative heights.
    hpg$widths <- grid:::unit.list(hpg$widths)
    hpg$widths[panelI] <- lapply(plot.lengths, unit, 'null')

    # add extra width between panels
    for(gap in 1:(length(panelI)-1)){
        hpg$widths[panelI[gap]+1]=list(unit(0.3, 'cm'))
    }

    # # facet coloring
    # names.grobs <- grid.ls(grid.force(hpg),print=FALSE)$name
    # strips <- names.grobs[which(grepl('strip.background',names.grobs))]
    # #cols <- colorRampPalette(brewer.pal(8,"Dark2"))(length(strips))
    # if(length(events$pheno %>% unique) == 1) {
    #     cols <- '#666666'
    # } else {
    #     cols <- events$pheno %>% unique
    # }

    # # modify grob object
    # for(strip in 1:length(strips)){
    #     hpg=editGrob(grid.force(hpg), gPath(strips[strip]), gp=gpar(fill=cols[strip]))
    # }

    # draw plot
    pdf(output.file, width, height, bg='white')
        grid.draw(hpg)
    dev.off()
}


#------------
# CNA heatmap
#------------

PlotCNHeatmap <- function(gene.cn, file.name, sample.names=NULL, threshold=FALSE) {

    if(is.null(sample.names) & threshold==TRUE) {
        sample.names <- gene.cn %>% select(matches('threshold')) %>% names %>% sort
    } else if(is.null(sample.names)) {
        sample.names <- gene.cn %>% names %>% list.filter(! . %in% c('hgnc','gene','chrom','start','mid','end','band'))
    }

    chr.rle <- gene.cn$chrom %>% rle
    chr.sep <- chr.rle$lengths %>% cumsum
    chr.mid <- c(0, chr.sep[-length(chr.sep)]) + chr.rle$lengths/2

    pdf(file.name, width=12, height=5*length(sample.names))

        g.cn <- gene.cn %>% select(one_of(rev(sample.names)))

        par(mfrow=c(length(sample.names),1), mar=c(8,5,1,2))
        image(as.matrix(g.cn), col=c("red", "darksalmon", "white", "lightblue", "blue"), xaxt='n', yaxt='n', zlim=c(-2, 2))
        box()

        for (i in (chr.sep*2)-1) {
            abline(v=i/((max(chr.sep)-1)*2), col='grey')
        }

        for (i in seq(-1, max(((2*(ncol(g.cn)-1))+1),1), 2)) {
            abline(h=i/(2*(ncol(g.cn)-1)), col="white", lwd=2)
        }

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
                 inset  = c(0, -0.4),
                 legend = c('Homozygous deletion', 'Loss', 'Gain', 'Amplification'),
                 fill   = c('red', 'darksalmon', 'lightblue', 'blue'),
                 xpd    = TRUE,
                 ncol   = 2 )
     dev.off()

}


#------------------------------#
#                              #
# INPUT & PARAMETER PROCESSING #
#                              #
#------------------------------#

#---------------------
# assign config params
#---------------------

# load config yaml file
config <- list.load(opts$project_config)

# load sample list
samples <- list.load(opts$samples_config) %>% list.stack %>% tbl_df

# sample key values
keys <- config$keys %>% unlist
if(is.null(keys)) {
    message(yellow('no keys found in config file, using samples file'))
    keys <- samples$name
    names(keys) <- samples$tumor
}

# define sub.sets
if(use.keys !=FALSE){
    sub.sets <- config$subsets %>% map(~ { keys[.x]})
} else {
    sub.sets <- config$subsets
}

# add all samples as sub.set
sub.sets <- c(sub.sets, list(all=samples$tumor), samples$tumor %>% as.list %>% set_names(samples$tumor %>% KeyMod(keys)))

# sub.set groups for pairwise comparisons
if(!is.null(config$subset_groups)) {
    sub.groups <-
        config$subset_groups %>%
        list.map(data_frame(.) %>% t %>% as_data_frame) %>%
        plyr::rbind.fill(.) %>%
        set_names(letters[1:(config$subset_groups %>% list.mapv(length(.)) %>% max)]) %>%
        tbl_df %>%
        gather(col,sub.set) %>%
        group_by(col) %>%
        mutate(group.id=row_number()) %>%
        ungroup %>%
        group_by(group.id) %>%
        mutate(subv=sub.sets[sub.set]) %>%
        mutate(overlap=anyDuplicated(unlist(subv))) %>%
        ungroup %>%
        select(-subv) %>%
        spread(col,sub.set) %>%
        mutate(comparison=names(config$subset_groups))

        if(any(sub.groups$overlap>0)){ message(red('warning: sub.set groups contain overlapping samples')) }

} else if(length(sub.sets) > 1) {
    sub.groups <-
        permutations(n=length(sub.sets), r=2, v=names(sub.sets)) %>%
        as_data_frame %>%
        set_names(c('a', 'b')) %>%
        mutate(group.id=row_number()) %>%
        select(group.id, a, b) %>%
        rowwise %>%
        mutate(overlap=length(intersect(unlist(sub.sets[a]), unlist(sub.sets[b])))) %>%
        ungroup %>%
        mutate(comparison=str_c(a, ' x ', b))
} else if(length(sub.sets) == 1) {
    sub.groups <- data_frame(overlap=NA, group.id=1, extension=names(config$subsets), a=names(config$subsets), b=NA)
} else {
    sub.groups <- data_frame(overlap=1, group.id=1, extension=NA, comparison=NA, a=NA, b=NA)
}

# stop if sub.set specifications absent
if(sub.groups %>% select(a, b) %>% unlist %in% names(sub.sets) %>% all == FALSE & length(sub.sets) > 1) {
    print((sub.groups %>% select(a, b) %>% unlist)[!sub.groups %>% select(a, b) %>% unlist %in% names(sub.sets)] %>% unname %>% unique)
    stop('missing sub.sets specified in sub.set groups')
}

# define sub groups
sub.groups %<>%
    filter(overlap==0 | is.na(overlap)) %>%
    select(-overlap) %>%
    bind_rows(data_frame(group.id=0, extension='all', comparison=NA, a='all'), .) %>%
    mutate(extension=ifelse(is.na(extension), str_c(comparison, group.id, sep='_') %>% FixName, extension)) %>%
    select(group.id, extension, comparison, everything())


#-----------------#
#                 #
# DATA PROCESSING #
#                 #
#-----------------#

#--------------------------------
# choose colors for pheno palette
#--------------------------------

# pheno color selection
if(random.pheno.color==TRUE){
    pheno.palette <-
        grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)] %>%
        sample(.,length(sub.sets)) %>%
        setNames(names(sub.sets))
} else {
    pheno.palette <-
        colorRampPalette(pheno.palette)(length(sub.sets)) %>%
        setNames(names(sub.sets))  # static palette
}

# build pheno bar lookup tables
sub.sets.pheno <-
    map2(sub.sets, pheno.palette, ~{ data_frame(sample=.x, .y=.y) %>% set_names(c('sample', .y)) }) %>%
    plyr::join_all(by='sample', type='full', match='all') %>%
    plyr::rename(replace=names(pheno.palette) %>% set_names(pheno.palette)) %>%
    tbl_df %>%
    KeyMod(keys)


#---------------------
# mutations processing
#---------------------

# read mutations file
if(length(grep('xlsx',opts$muts_table_in))){
    muts <- tryCatch({
        read.xlsx(opts$muts_table_in, 'MUTATION_SUMMARY') %>% tbl_df},
        error=function(e){read.xlsx(opts$muts_table_in) %>% tbl_df
    })
} else {
    muts <- read.delim(opts$muts_table_in, sep='\t', stringsAsFactors=FALSE) %>% tbl_df
}


# format muts with cleaning function
muts <-
    muts %>%
    FormatEvents(keys=keys) %>%
    mutate(mut.id=str_c(sample, ':', row_number())) %>%
    select(sample, chrom, pos, gene, effect, everything())


# remove silent mutations if present
if(include.silent!=TRUE){
    muts %<>% filter(effect!='Silent')
}


#-------------------
# ABSOLUTE CCF calls
#-------------------

if(call.abs == TRUE) {
    # retreive absolute ccf calls
    abs.ccf <-
        list.files(opts$seg_maf_path, pattern='_ABS_MAF.txt', full.names=TRUE) %>%
        map(~ { read.delim(.x, sep='\t',stringsAsFactors=FALSE) %>% tbl_df }) %>%
        bind_rows %>%
        separate(sample, into=c('sample','normal'), sep='_') %>%
        rename(pos=Start_position, alt.dept=alt) %>%
        FormatEvents %>%
        mutate(clonal=ifelse(pr.somatic.clonal>=0.5 | ci95.low >=0.9, 'Clonal', 'Subclonal')) %>%
        select(sample,gene,chrom,pos,ccf,clonal,purity) %>%
        mutate(sample=keys[sample])

    # add absolute calls to mutation table
    muts %<>% left_join(abs.ccf)
}


#----------
# LOH calls
#----------

if(call.loh == TRUE) {
    # read cnf & make calls
    segs.loh <-
        list.files(opts$cncf_path, pattern='.cncf.txt',full.names=TRUE) %>%
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
                mutate(seg.id=str_c(sample, ':', row_number())) %>%
                select(sample,chrom,loc.start,loc.mid,loc.end,loh,lcn.em,lcn,tcn.em,tcn,mafR,cnlr.median,cnlr.median.clust,cf,seg,num.mark,seg.id)

            return(cncf)
        }) %>%
        bind_rows %>%
        mutate(sample=keys[sample])


    # assign loh calls for each mutation
    muts.loh <-
        muts %>%
        select(mut.id, sample, chrom, pos) %>%
        group_by(mut.id) %>%
        full_join(segs.loh, by=c('sample', 'chrom')) %>%
        slice(which.min(pos-loc.mid)) %>%
        mutate(in.seg=(loc.start<=pos & pos<=loc.end)) %>%
        ungroup


    # filter loh calls if specified
    if(loh.closest!=TRUE) {
        muts.loh %<>% filter(in.seg==TRUE)
    }

    if('loh' %in% colnames(muts)) {
        message(yellow('warning: removing existing loh column'))
        muts %<>% select(-loh)
    }

    # add loh to mutation tables
    muts %<>% left_join(muts.loh)
}

# write muts file
muts %>%
arrange(sample, desc(ccf)) %>%
write_tsv(opts$muts_table_out)


#----------------------------------------
# format per-gene copy number information
#----------------------------------------

# select threshold columns
gene.cn <- read.delim(opts$gene_cn_in, sep='\t',stringsAsFactors=FALSE) %>% tbl_df

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


# sub.set amp / del rows
cnas <-
    gene.cn %>%
    mutate(amp=rowSums(.[cn.sample.keys]==2)) %>%
    mutate(del=rowSums(.[cn.sample.keys]==-2)) %>%
    filter(amp==TRUE|del==TRUE) %>%
    mutate_each(funs(ifelse(.==1,0,.)), one_of(cn.sample.keys)) %>%
    mutate_each(funs(ifelse(.==-1,0,.)), one_of(cn.sample.keys)) %>%
    unique %>%
    arrange(chrom, start) %>%
    group_by(chrom) %>%
    mutate(max.end=max(end)) %>%
    mutate(lag.end=lag(end)) %>%
    (function(df=., col.v=c("chrom", "band", unname(cn.sample.keys))) {
        DistinctCalls <- sapply(col.v, function(col) {
            lazyeval::interp(~ col.name, col.name=as.name(col))
        })
        df %>% distinct_(.dots=DistinctCalls, .keep_all=TRUE)
    }) %>%
    mutate(lead.lag.end=lead(lag.end)) %>%
    mutate(end=ifelse(!is.na(lead.lag.end), lead.lag.end, max.end)) %>%
    select(chrom, start, end, band, one_of(cn.sample.keys)) %>%
    gather(sample, effect, -band, -chrom, -start, -end) %>%
    filter(effect==2|effect==-2) %>%
    rowwise %>%
    mutate(sample=strsplit(sample,'_') %>% unlist %>% head(1)) %>%
    ungroup %>%
    unique %>%
    FormatEvents %>%
    select(sample, band, chrom, start, end, effect)


# write CNA file
cnas %>%
spread(sample, effect) %>%
arrange(chrom, start) %>%
write_tsv(opts$cnas_table_out)


#---------------------------
# combined muts / cnas table
#---------------------------

# all mutations + cnas, add pheno columns
variants <-
    bind_rows( muts %>% mutate(class='muts') %>% select(class, sample, chrom, pos, span=gene, effect, ccf, loh, clonal),
               cnas %>% mutate(class='cnas') %>% select(class, sample, chrom, pos=start, span=band, effect) ) %>%
    left_join(sub.sets.pheno, by='sample')


# write variants file
muts %>%
arrange(sample, effect) %>%
write_tsv(opts$variants_table_out)


#----------------------#
#                      #
# per-sub set plotting #
#                      #
#----------------------#

#--------------------
# CN heatmap plotting
#--------------------

for(sub.set.num in 1:length(sub.sets)) {

    sub.set.name <- names(sub.sets)[sub.set.num] %>% FixName
    sub.set <- sub.sets[[sub.set.num]] %>% KeyMod(keys)

    # cut down gene.cn table
    sub.gene.cn <- gene.cn %>% select(gene,chrom,start,end,band,one_of(sub.set))

    # call CN heatmap plotting function
    PlotCNHeatmap(sub.gene.cn, file.name=str_c('summary/cna_heatmap_',sub.set.name,'.pdf'), threshold=FALSE)
}


#----------------#
#                #
# sub group loop #
#                #
#----------------#

for (sub.num in 0:(nrow(sub.groups)-1)) {

    #-----------------
    # setup sub groups
    #-----------------

    # sub group setup
    sub.group <- sub.groups %>% filter(group.id == sub.num)
    sub.group.sets <- sub.group[names(sub.group) %>% list.filter(!. %in% c('group.id', 'comparison', 'extension'))] %>% list.filter(!is.na(.)) %>% unlist
    sub.group.samples <- sub.group.sets %>% list.map(sub.sets[.]) %>% unlist %>% unique %>% KeyMod(keys)

    # create pheno column & filter samples
    sub.variants <-
        variants %>%
        MergePheno(merge.cols=sub.group.sets, col.name='pheno') %>%
        MergePheno(merge.cols=sub.group.samples, col.name='sample.colors') %>%
        select(class,sample,chrom,pos,span,effect,ccf,loh,clonal,pheno, sample.colors) %>%
        filter(!is.na(pheno))

    # melted mutations table for sub.set pair
    sub.muts <-
        sub.variants %>%
        filter(class=='muts') %>%
        rename(gene=span)

    # melted cna table for sub.set pair
    sub.cnas <-
        sub.variants %>%
        filter(class=='cnas') %>%
        rename(band=span)

    #  sub variants tree
    sub.variants.tree <- sub.variants %>% MeltToTree(dist.method='hamming', clust.method='complete', sort.method='distance', span='span', tree.samples=sub.group.samples)

    # sub muts tree
    sub.muts.tree <- sub.muts %>% MeltToTree(dist.method='hamming', clust.method='complete', sort.method='distance', span='gene', tree.samples=sub.group.samples)

    # sub cnas tree
    sub.cnas.tree <- sub.cnas %>% MeltToTree(dist.method='hamming', clust.method='complete', sort.method='distance', span='band', tree.samples=sub.group.samples)


    #-----------------
    # heatmap plotting
    #-----------------

    # mutations (by type)
    sub.muts %>%
    OrgEvents(sample.order=labels(sub.muts.tree$dend), pheno.order=sub.group.sets, sub.sets.pheno=sub.sets.pheno, recurrence=1, allosome='merge', event.type='gene') %>%
    PlotVariants( output.file = str_c('summary/mutation_heatmap_type_', sub.group$extension, '.pdf'),
                  clonal      = FALSE,
                  pathogenic  = FALSE,
                  ccf         = FALSE,
                  loh         = FALSE,
                  event.type  = 'gene',
                  width       = length(unique(.$sample)) + 10,
                  height      = (length(unique(.$gene))/3) + 10,
                  text.size   = 30 )

    # mutations (by type)
    sub.muts %>%
    OrgEvents(sample.order=labels(sub.muts.tree$dend), pheno.order=sub.group.sets, sub.sets.pheno=sub.sets.pheno, recurrence=2, allosome='merge', event.type='gene') %>%
    PlotVariants( output.file = str_c('summary/mutation_heatmap_type_recurrent_', sub.group$extension, '.pdf'),
                  clonal      = FALSE,
                  pathogenic  = FALSE,
                  ccf         = FALSE,
                  loh         = FALSE,
                  event.type  = 'gene',
                  width       = length(unique(.$sample)) + 10,
                  height      = (length(unique(.$gene))/3) + 10,
                  text.size   = 30 )

    # mutations ccf
    sub.muts %>%
    OrgEvents(sample.order=labels(sub.muts.tree$dend), pheno.order=sub.group.sets, sub.sets.pheno=sub.sets.pheno, recurrence=1, allosome='merge', event.type='gene') %>%
    PlotVariants( output.file = str_c('summary/mutation_heatmap_ccf_', sub.group$extension, '.pdf'),
                  clonal      = TRUE,
                  pathogenic  = FALSE,
                  ccf         = TRUE,
                  loh         = TRUE,
                  event.type  = 'gene',
                  width       = length(unique(.$sample)) + 10,
                  height      = (length(unique(.$gene))/3) + 10,
                  text.size   = 30 )

    # mutations ccf [recurrent]
    sub.muts %>%
    OrgEvents(sample.order=labels(sub.muts.tree$dend), pheno.order=sub.group.sets, sub.sets.pheno=sub.sets.pheno, recurrence=2, allosome='merge', event.type='gene') %>%
    PlotVariants( output.file = str_c('summary/mutation_heatmap_ccf_recurrent_', sub.group$extension, '.pdf'),
                  clonal      = TRUE,
                  pathogenic  = FALSE,
                  ccf         = TRUE,
                  loh         = TRUE,
                  event.type  = 'gene',
                  width       = length(unique(.$sample)) + 10,
                  height      = (length(unique(.$gene))/3) + 10,
                  text.size   = 30 )

    # copy number
    sub.cnas %>%
    OrgEvents(sample.order=labels(sub.cnas.tree$dend), pheno.order=sub.group.sets, sub.sets.pheno=sub.sets.pheno, recurrence=1, allosome='merge', event.type='band') %>%
    PlotVariants( output.file = str_c('summary/cna_heatmap_', sub.group$extension, '.pdf'),
                  clonal      = FALSE,
                  pathogenic  = FALSE,
                  ccf         = FALSE,
                  loh         = FALSE,
                  event.type  = 'band',
                  width       = length(unique(.$sample)) + 10,
                  height      = (length(unique(.$band))/3) + 10,
                  text.size   = 30 )

    # all variants
    sub.variants %>%
    OrgEvents(sample.order=labels(sub.variants.tree$dend), pheno.order=sub.group.sets, sub.sets.pheno=sub.sets.pheno, recurrence=1, allosome='merge', event.type='span') %>%
    PlotVariants( output.file = str_c('summary/variant_heatmap_type_', sub.group$extension, '.pdf'),
                  clonal      = FALSE,
                  pathogenic  = FALSE,
                  ccf         = FALSE,
                  loh         = FALSE,
                  event.type  = 'span',
                  width       = length(unique(.$sample)) + 10,
                  height      = (length(unique(.$span))/3) + 10,
                  text.size   = 30 )


    #-----------------------------#
    #                             #
    # per-sub group pair plotting #
    #                             #
    #-----------------------------#

    if(!is.na(sub.group$b)) {

        #---------------------------
        # Fisher's exact CN plotting
        #---------------------------

        samples.a <- sub.sets[sub.group$a] %>% unlist
        samples.b <- sub.sets[sub.group$b] %>% unlist

        # copy number plotting
        gene.cn.a <- gene.cn[c('gene', 'chrom', 'start', 'end', samples.a)]
        gene.cn.b <- gene.cn[c('gene', 'chrom', 'start', 'end', samples.b)]

        # call fishers plots for copy number
        Fisher( plot.type       = 'copy number',
                gene.matrix.a   = gene.cn.a,
                gene.matrix.b   = gene.cn.b,
                plot.title.main = sub.group$comparison,
                plot.title.a    = sub.group$a,
                plot.title.b    = sub.group$b,
                allosome        = allosome,
                targets.file    = targets.file,
                suffix          = sub.group$extension,
                threshold.a     = FALSE,
                threshold.b     = FALSE,
                gene.names      = gene.names )

        #---------------------------------
        # Fisher's exact mutation plotting
        #---------------------------------

        # mutation plotting
        gene.muts.a <- sub.muts %>% select(sample, gene, chrom, pos, effect) %>% mutate(effect=1) %>% filter(sample %in% samples.a) %>% spread(sample, effect, fill=0)
        gene.muts.b <- sub.muts %>% select(sample, gene, chrom, pos, effect) %>% mutate(effect=1) %>% filter(sample %in% samples.b) %>% spread(sample, effect, fill=0)

        # call fishers function for mutations
        Fisher( plot.type       = 'mutation',
                gene.matrix.a   = gene.muts.a,
                gene.matrix.b   = gene.muts.b,
                plot.title.main = sub.group$comparison,
                plot.title.a    = sub.group[1],
                plot.title.b    = sub.group[2],
                allosome        = allosome,
                targets.file    = targets.file,
                suffix          = sub.group$extension,
                gene.names      = gene.names )


        #-------#
        #       #
        # TREES #
        #       #
        #-------#

        # remove labels from tree
        if(tree.labels == FALSE){
            suppressWarnings(sub.muts.tree$dend %<>% set('labels', ''))
            suppressWarnings(sub.cnas.tree$dend %<>% set('labels', ''))
            suppressWarnings(sub.variants.tree$dend %<>% set('labels', ''))
        }


        #--------------
        # distance tree
        #--------------

        pdf(str_c('summary/genomic_distance_tree_ladder_', sub.group$extension, '.pdf'), 140, 38)
            par(mar=c(10, 10, 10, 10))  # bottom, left, top, right
            layout(matrix(c(1, 2), nrow=1), widths=c(10, 1))
            plot(sub.muts.tree$dend, cex.axis=6)

            colored_bars( colors               = sub.muts %>% select(pheno, sample.colors),
                          dend                 = sub.muts.tree$dend,
                          sort_by_labels_order = FALSE,
                          add                  = TRUE,
                          rowLabels            = sub.muts$sample,
                          y_scale              = 10,
                          cex.rowLabels        = 5.8 )

            legend(x=3, y=3, legend=unique(sub.muts$sample), fill=unique(sub.muts$sample.colors), cex=4)
            legend(x=3, y=3, legend=unique(sub.muts$pheno), fill=unique(sub.group$comparison), cex=4)
        dev.off()


        #--------------
        # weighted tree
        #--------------

        pdf(str_c('summary/genomic_distance_tree_weighted_', sub.group$extension, '.pdf'),60,25)
            suppressWarnings(heatmap.2(
                sub.muts.tree$event.matrix,
                trace='none',
                Rowv=FALSE,
                hclustfun=function(x){hclust(x, 'ward.D2')},
                col=c('white','black'),
                reorderfun=function(d,w) rev(reorder(d,w)),
                distfun=function(x) as.dist(Hamming(x)),
                key=FALSE ))
        dev.off()



        pdf(str_c('summary/genomic_distance_tree_bw_', sub.group$extension, '.pdf'),60,25)
            heatmap.2(
                t(sub.muts.tree$event.matrix),
                trace='none',
                dendrogram='column',
                Colv=rev(sub.muts.tree$dend),
                col=c('white','black'),
                key=FALSE
            )
        dev.off()

        #-------------------------------
        # MUTANTS LINEAGE MAP GENERATION
        #-------------------------------

        lineage.matrix = cbind(sub.muts.tree$event.matrix, Parental=FALSE) %>% t

        phylo.data <- phyDat(lineage.matrix, type="USER", levels=c(0, 1))
        phylo.hamming <- dist.hamming(phylo.data)
        phylo.tree <- njs(phylo.hamming)
        phylo.ratchet <- pratchet(phylo.data, start=phylo.tree) %>% acctran(phylo.data)
        lineage.dend <- root(phylo.ratchet, "Parental")

        pdf(str_c('summary/variant_lineage_muts_', sub.group$extension, '.pdf'))
            plot(lineage.dend)
        dev.off()

        #----------------------------
        # CNAS LINEAGE MAP GENERATION
        #----------------------------

        lineage.matrix = cbind(sub.cnas.tree$event.matrix, Parental=FALSE) %>% t

        phylo.data <- phyDat(lineage.matrix, type="USER", levels=c(0, 1))
        phylo.hamming <- dist.hamming(phylo.data)
        phylo.tree <- njs(phylo.hamming)
        phylo.ratchet <- pratchet(phylo.data, start=phylo.tree) %>% acctran(phylo.data)
        lineage.dend <- root(phylo.ratchet, "Parental")

        pdf(str_c('summary/variant_lineage_cnas_', sub.group$extension, '.pdf'))
            plot(lineage.dend)
        dev.off()

        #-------------------------------
        # VARIANT LINEAGE MAP GENERATION
        #-------------------------------

        lineage.matrix = cbind(sub.variants.tree$event.matrix, Parental=FALSE) %>% t

        phylo.data <- phyDat(lineage.matrix, type="USER", levels=c(0, 1))
        phylo.hamming <- dist.hamming(phylo.data)
        phylo.tree <- njs(phylo.hamming)
        phylo.ratchet <- pratchet(phylo.data, start=phylo.tree) %>% acctran(phylo.data)
        lineage.dend <- root(phylo.ratchet, "Parental")

        pdf(str_c('summary/variant_lineage_variants_', sub.group$extension, '.pdf'))
            plot(lineage.dend)
        dev.off()

    }

}
