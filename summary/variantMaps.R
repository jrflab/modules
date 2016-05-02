#!/usr/bin/env Rscript

# script for build mutation heatmaps, trees & other standard plots

#-----------#
#           #
# LIBRARIES #
#           #
#-----------#

pacman::p_load( dplyr,lazyeval,readr,tidyr,magrittr,purrr,stringr,rlist,openxlsx, # base
                crayon,colorspace,RColorBrewer, # coloring
                ggplot2,grid,gridExtra,gplots, # plot layout
                dendextend,dendextendRcpp,dynamicTreeCut,gclus,phangorn, # dentrogram
                gtools,digest ) # permutation & hashing

loadNamespace('plyr')

source('modules/summary/variantFishers.R')

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
    message('Need geneCN file')
    print_help(parser)
    stop()
else if (is.null(opt$mutationSummary)) {
    message('Need mutation summary file')
    print_help(parser)
    stop()
} else if (is.null(opt$mutOutFile)) {
    message('Need mut output file')
    print_help(parser)
    stop()
} else if (is.null(opt$mutRecurrentOutFile)) {
    message('Need mut recurrent output file')
    print_help(parser)
    stop()
} else if (is.null(opt$cnOutFile)) {
    message('Need cn output file')
    print_help(parser)
    message()
} else if (is.null(opt$cnRecurrentOutFile)) {
    message('Need cn recurrent output file')
    print_help(parser)
    stop()
} else if (is.null(opt$cnAmpDeltOutFile)) {
    message('Need cn amp del output file')
    print_help(parser)
    stop()
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

#------------
# run options
#------------

include.silent       = FALSE
dist.method          = 'hamming'   # 'hamming', euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'     | see ?dist for details [hamming method implemented manually]
clust.method         = 'complete' # 'complete', ward', 'single', 'complete', 'average', 'mcquitty', 'median', 'centroid' | see ?hclust for details
exclude.values       = c('', '.', 'Normal', 'Not Performed', 'Performed but Not Available', FALSE)
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
muts.out             = 'summary/mutation_heatmap.tsv'
# muts.file            = 'pub_table.tsv'


#---------------------
# assign config params
#---------------------

# load config yaml file
config <- list.load('subsets_config.yaml')

# sample key values
keys <- config$keys %>% unlist

# define subsets
if(use.keys !=FALSE){
    subsets <- config$subsets %>% map(~ { keys[.x]})
} else {
    subsets <- config$subsets
}

# subset groups for pairwise comparisons
if(!is.null(config$subset_groups)) {
    subset.groups <-
        config$subset_groups %>%
        list.map(data_frame(.) %>% t %>% as_data_frame) %>%
        plyr::rbind.fill(.) %>%
        set_names(letters[1:(config$subset_groups %>% list.mapv(length(.)) %>% max)]) %>%
        tbl_df %>%
        gather(col,subset) %>%
        group_by(col) %>%
        mutate(group.id=row_number()) %>%
        ungroup %>%
        group_by(group.id) %>%
        mutate(subv=subsets[subset]) %>%
        mutate(overlap=anyDuplicated(unlist(subv))) %>%
        ungroup %>%
        select(-subv) %>%
        spread(col,subset)

        if(any(subset.groups$overlap>0)){
            message(red('warning: subset groups contain overlapping samples'))
        }

} else {
    subset.groups <-
        permutations(n=length(config$subsets), r=2, v=names(config$subsets)) %>%
        as_data_frame %>%
        set_names(c('a', 'b')) %>%
        mutate(group.id=row_number()) %>%
        select(group.id, a, b) %>%
        rowwise %>%
        mutate(overlap=length(intersect(unlist(subsets[a]), unlist(subsets[b])))) %>%
        ungroup
}

# stop if subset specifications absent
if(subset.groups %>% select(a, b) %>% unlist %in% names(subsets) %>% all == FALSE) {
    print((subset.groups %>% select(a, b) %>% unlist)[!subset.groups %>% select(a, b) %>% unlist %in% names(subsets)] %>% unname %>% unique)
    stop('missing subsets specified in subset groups')
}

subset.groups %<>%
    filter(overlap==0) %>%
    select(-overlap)

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

#--------------------------
# Hamming distance function
#--------------------------

Hamming <- function(event.matrix) {
    D <- (1 - event.matrix) %*% t(event.matrix)
    D + t(D)
}


#------------------
# character hashing
#------------------

hash <- function(x) {
    hstr <- digest(x, algo='xxhash32')
    as.numeric(paste0('0x', hstr)) %% 80
}


#-------------------------------------
# melted event table into ordered tree
#-------------------------------------

MeltToTree <- function(event.melt, dist.method='hamming', clust.method='complete', sort.method='distance', span='gene') {

    event.melt %<>% rename_(span=span)

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

    # build master event matrix
    event.matrix <- event.melt %>% MeltToMatrix

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

    # function to dummy column if absent
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

    # convert chrom vector to desired format
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

    # re-type colummns
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
        keys.index <- which(events$sample %in% names(keys))
        if(length(keys.index) == 0) {
            message(green('no keys to convert'))
        } else if(force.keys==FALSE) {
            conversions <- data.frame(pre=as.character(events$sample[keys.index]), post=keys[as.character(events$sample[keys.index])], row.names=NULL) %>% unique
            message(green(str_c('key conversions: ',length(conversions$pre),'/',length(unique(events$sample)))))
            print(conversions)
            events %<>% mutate(sample=ifelse(row_number() %in% keys.index, keys[sample], sample))
        } else {
            conversions <- data.frame(pre=as.character(unique(events$sample)), post=keys[unique(as.character(events$sample))], in.keys=unique(as.character(events$sample)) %in% names(keys), row.names=NULL)
            message(yellow('forcing key conversions'))
            print(conversions)
            events %<>% mutate(sample=keys[sample])
        }
    }

    return(events)
}


#----------------------------------
# prepare melted table for plotting
#----------------------------------

OrgEvents <- function(events, sample.order, pheno.order=NULL, subsets.pheno, recurrence=1, allosome='merge', event.type='gene') {

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
    DummyCols <- function(events, col.names) {
        for(col.name in col.names) {
            if(!col.name %in% colnames(events)) {
                message(green(str_c('adding dummy column: ', col.name)))
                events.names <- colnames(events) %>% replace(which(is.na(.)),'empty') %>% list.filter(.!='add.col')
                events$add.col <- NA
                events %<>% set_names(c(events.names,col.name))
            }
        }
        return(events)
    }

    # convert chrom vector to desired format
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

    missing.samples <- sample.order[!sample.order %in% events$sample]

    if(length(missing.samples) > 0){

        message(yellow(str_c('missing specified samples: ',missing.samples,'\n')))

        missing.fill <-
            subsets.pheno %>%
            select(sample,one_of(pheno.order)) %>%
            set_names(c('sample','a','b')) %>%
            mutate(pheno=ifelse(!is.na(a),a,ifelse(!is.na(b),b,NA))) %>%
            filter(!is.na(pheno)) %>%
            select(sample, pheno) %>%
            filter(sample %in% missing.samples) %>%
            full_join(data_frame(missing.samples,unique(unlist(events[,event.type]))) %>%
            set_names(c('sample', event.type)))
    }

    # specify output columns
    col.names <- c('sample','gene','chrom','effect','pheno','band','span','pos','maf','ccf','loh','clonal','pathogenic')

    # facet ordering & color assignment
    if(!is.null(pheno.order)) {
        cols <- subsets.pheno %>% select(one_of(pheno.order)) %>% summarise_each(funs(max(.,na.rm=TRUE))) %>% unlist
        events %<>% mutate(pheno=factor(pheno, levels=cols))
    } else if('pheno' %in% names(events)) {
        cols <- '#666666'
    } else {
        events %<>% mutate(pheno='Variants')
    }

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
            ifelse(ccf>0.00 & ccf<=0.05, '0% < CCF ≤ 5%',
            ifelse(ccf>0.05 & ccf<=0.20, '5% < CCF ≤ 20%',
            ifelse(ccf>0.20 & ccf<=0.40, '20% < CCF ≤ 40%',
            ifelse(ccf>0.40 & ccf<=0.60, '40% < CCF ≤ 60%',
            ifelse(ccf>0.60 & ccf<=0.80, '60% < CCF ≤ 80%',
            ifelse(ccf>0.80,             '80% < CCF ≤ 100%',
            NA)))))))) %>%
        mutate(ccf=ifelse(is.na(effect),NA,ccf)) %>%
        mutate(ccf=factor(ccf, levels=c('CCF = 0%','0% < CCF ≤ 5%','5% < CCF ≤ 20%','20% < CCF ≤ 40%','40% < CCF ≤ 60%','60% < CCF ≤ 80%','80% < CCF ≤ 100%'))) %>%
        mutate(pathogenic=as.factor(ifelse(pathogenic=='Pathogenic', 'darkgoldenrod2', NA))) %>%
        mutate(clonal=ifelse(clonal=='Clonal',clonal,NA)) %>%
        select(-order) %>%
        mutate(sample=factor(sample, levels=sample.order))

    if(event.type=='band') {
        events %<>%
            arrange(!is.na(effect), desc(n.band), desc(chrom), desc(band)) %>%
            mutate(band=factor(band, levels=filter(.,!is.na(effect)) %>% .$band %>% unique))
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
                    `80% < CCF ≤ 100%`='#0b3269' )

    geometry <- c(`Loss of heterozygosity`=3)
    surround.clonal <- c(Clonal='darkgoldenrod2')
    surround.pathogenic <- c(Pathogenic='purple')

    # main plot params
    if(event.type=='gene') {
        hp <- ggplot(events, aes(Sample,Gene))
    } else if(event.type=='band') {
        hp <- ggplot(events, aes(Sample,Band))
    } else {
        hp <- ggplot(events, aes(Sample,Span))
    }

    if(ccf==TRUE){  # CCF coloring
        hp <- hp +
        geom_tile(data=events, aes(fill=CCF), colour='white') +
        scale_fill_manual(breaks=names(palette), values=palette, na.value='gray90', drop=FALSE)
    } else {  # draw tiles and color
        hp <- hp +
        geom_tile(data=events, aes(fill=Effect), colour='white') +
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

    # facet coloring
    names.grobs <- grid.ls(grid.force(hpg),print=FALSE)$name
    strips <- names.grobs[which(grepl('strip.background',names.grobs))]
    #cols <- colorRampPalette(brewer.pal(8,"Dark2"))(length(strips))
    if(length(events$pheno %>% unique) == 1) {
        cols <- '#666666'
    } else {
        cols <- events$pheno %>% unique
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

#--------------------------------
# choose colors for pheno palette
#--------------------------------

# random color selection
set.seed(color.seed)
if(random.pheno.color==TRUE){
    pheno.palette <-
        grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)] %>%
        sample(.,length(subsets)) %>%
        setNames(names(subsets))
} else {
    pheno.palette <-
        colorRampPalette(pheno.palette)(length(subsets)) %>%
        setNames(names(subsets))  # static palette
}

# build pheno bar lookup tables
subsets.pheno <-
    map2(subsets, pheno.palette, ~{ data_frame(sample=.x, .y=.y) %>% set_names(c('sample', .y)) }) %>%
    plyr::join_all(by='sample', type='full', match='all') %>%
    plyr::rename(replace=names(pheno.palette) %>% set_names(pheno.palette)) %>%
    tbl_df %>%
    bind_rows(config$keys %>% unlist %>% list.filter(!. %in% unlist(subsets)) %>% data_frame(sample=.))


#---------------------
# mutations processing
#---------------------

# read mutations file
if(length(grep('xlsx',muts.file))){
    muts <- tryCatch({
        read.xlsx(muts.file, 'MUTATION_SUMMARY') %>% tbl_df},
        error=function(e){read.xlsx(muts.file) %>% tbl_df
    })
} else {
    muts <- read.delim(muts.file, sep='\t', stringsAsFactors=FALSE) %>% tbl_df
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

# join pheno columns to mutation tree
muts %<>% left_join(subsets.pheno, by='sample')


#-------------------
# ABSOLUTE CCF calls
#-------------------

# retreive absolute ccf calls
abs.ccf <-
    list.files('absolute/reviewed/SEG_MAF', pattern='_ABS_MAF.txt', full.names=TRUE) %>%
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


#----------
# LOH calls
#----------

# read cnf & make calls
segs.loh <-
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


# add loh to mutation tables
muts %<>% left_join(muts.loh)


#------------------------
# construct mutation tree
#------------------------

muts.tree <- muts %>% MeltToTree(dist.method='hamming', clust.method='complete', sort.method='distance', span='gene')


#----------------------------------------
# format per-gene copy number information
#----------------------------------------

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
    FormatEvents %>%
    select(sample, band, chrom, start, end, effect)

# write CNA file
cnas %>%
spread(sample, effect) %>%
arrange(chrom, start) %>%
write_tsv(cn.amp.del.out)

# add pheno columns
cnas %<>% left_join(subsets.pheno)


#--------------------------------
# construct copy number event tree
#--------------------------------

cnas.tree <- cnas %>% MeltToTree(dist.method='hamming', clust.method='complete', sort.method='distance', span='band')


#---------------------------
# construct all-variant tree
#---------------------------

# all mutations and cnas
variants <-
    bind_rows( muts %>% select(sample, chrom, pos, span=gene, effect, ccf, loh, clonal),
               cnas %>% select(sample, chrom, pos=start, span=band, effect) )

# all variants tree
variants.tree <- variants %>% MeltToTree(dist.method='hamming', clust.method='complete', sort.method='distance', span='span')


#---------------------------#
#                           #
# build variant & CCF plots #
#                           #
#---------------------------#

# mutations (by type)
muts %>%
mutate(pheno='Mutations by type') %>%
OrgEvents(sample.order=labels(muts.tree$dend), pheno.order=NULL, subsets.pheno=subsets.pheno, recurrence=1, allosome='merge', event.type='gene') %>%
PlotVariants( output.file = 'summary/mutation_heatmap_type.pdf',
              clonal      = FALSE,
              pathogenic  = FALSE,
              ccf         = FALSE,
              loh         = FALSE,
              event.type  = 'gene',
              width       = length(unique(.$sample)) + 10,
              height      = (length(unique(.$gene))/3) + 10,
              text.size   = 30 )

# mutations (by type)
muts %>%
mutate(pheno='Recurrent mutations by type') %>%
OrgEvents(sample.order=labels(muts.tree$dend), pheno.order=NULL, subsets.pheno=subsets.pheno, recurrence=2, allosome='merge', event.type='gene') %>%
PlotVariants( output.file = 'summary/mutation_heatmap_type_recurrent.pdf',
              clonal      = FALSE,
              pathogenic  = FALSE,
              ccf         = FALSE,
              loh         = FALSE,
              event.type  = 'gene',
              width       = length(unique(.$sample)) + 10,
              height      = (length(unique(.$gene))/3) + 10,
              text.size   = 30 )

# mutations ccf
muts %>%
mutate(pheno='Mutations CCF') %>%
OrgEvents(sample.order=labels(muts.tree$dend), pheno.order=NULL, subsets.pheno=subsets.pheno, recurrence=1, allosome='merge', event.type='gene') %>%
PlotVariants( output.file = 'summary/mutation_heatmap_ccf.pdf',
              clonal      = TRUE,
              pathogenic  = FALSE,
              ccf         = TRUE,
              loh         = TRUE,
              event.type  = 'gene',
              width       = length(unique(.$sample)) + 10,
              height      = (length(unique(.$gene))/3) + 10,
              text.size   = 30 )

# mutations ccf [recurrent]
muts %>%
mutate(pheno='Recurrent mutations CCF') %>%
OrgEvents(sample.order=labels(muts.tree$dend), pheno.order=NULL, subsets.pheno=subsets.pheno, recurrence=2, allosome='merge', event.type='gene') %>%
PlotVariants( output.file = 'summary/mutation_heatmap_ccf_recurrent.pdf',
              clonal      = TRUE,
              pathogenic  = FALSE,
              ccf         = TRUE,
              loh         = TRUE,
              event.type  = 'gene',
              width       = length(unique(.$sample)) + 10,
              height      = (length(unique(.$gene))/3) + 10,
              text.size   = 30 )

# copy number
cnas %>%
mutate(pheno='Copy number') %>%
OrgEvents(sample.order=labels(cnas.tree$dend), pheno.order=NULL, subsets.pheno=subsets.pheno, recurrence=1, allosome='merge', event.type='band') %>%
PlotVariants( output.file = 'summary/cna_heatmap.pdf',
              clonal      = FALSE,
              pathogenic  = FALSE,
              ccf         = FALSE,
              loh         = FALSE,
              event.type  = 'band',
              width       = length(unique(.$sample)) + 10,
              height      = (length(unique(.$band))/3) + 10,
              text.size   = 30 )

# all variants
variants %>%
mutate(pheno='Recurrent mutations CCF') %>%
OrgEvents(sample.order=labels(variants.tree$dend), pheno.order=NULL, subsets.pheno=subsets.pheno, recurrence=1, allosome='merge', event.type='span') %>%
PlotVariants( output.file = 'summary/variant_heatmap_type.pdf',
              clonal      = FALSE,
              pathogenic  = FALSE,
              ccf         = FALSE,
              loh         = FALSE,
              event.type  = 'span',
              width       = length(unique(.$sample)) + 10,
              height      = (length(unique(.$span))/3) + 10,
              text.size   = 30 )


#-------------------------#
#                         #
# subset comparison plots #
#                         #
#-------------------------#

system("mkdir summary/subsets &>/dev/null")

for (sub.num in 1:nrow(subset.groups)) {

    # format file name
    sub.name <- subset.groups$group.id[sub.num]
    sub.ext <- gsub(' ', '_', subset.groups$group.id[sub.num])
    sub.pair <- c(subset.groups[sub.num,'a'], subset.groups[sub.num,'b']) %>% unlist

    #-----------------
    # heatmap plotting
    #-----------------

    # melted mutations table for subset pair
    sub.muts <-
        muts %>%
        select(sample,chrom,pos,gene,effect,loh,ccf,clonal,one_of(sub.pair)) %>%
        set_names(c('sample','chrom','pos','gene','effect','loh','ccf','clonal','a','b')) %>%
        mutate(pheno=ifelse(!is.na(a), a,
                     ifelse(!is.na(b), b, NA))) %>%
        filter(!is.na(pheno)) %>%
        select(-a, -b)

    sub.muts.tree <- sub.muts %>% MeltToTree(dist.method='hamming', clust.method='complete', sort.method='distance', span='gene')

    # melted cna table for subset pair
    sub.cnas <-
        cnas %>%
        select(sample,band,chrom,start,end,effect,one_of(sub.pair)) %>%
        set_names(c('sample','band','chrom','start','end','effect','a','b')) %>%
        mutate(pheno=ifelse(!is.na(a), a,
                     ifelse(!is.na(b), b, NA))) %>%
        filter(!is.na(pheno)) %>%
        select(-a, -b)

    sub.cnas.tree <- sub.cnas %>% MeltToTree(dist.method='hamming', clust.method='complete', sort.method='distance', span='band')

    # melted variant table for subset pair
    sub.variants <-
        bind_rows( sub.muts %>% select(sample, chrom, pos, span=gene, effect, ccf, loh, clonal),
                   sub.cnas %>% select(sample, chrom, pos=start, span=band, effect) )

    # all variants tree
    sub.variants.tree <- variants %>% MeltToTree(dist.method='hamming', clust.method='complete', sort.method='distance', span='span')


    # mutations (by type)
    sub.muts %>%
    OrgEvents(sample.order=labels(sub.muts.tree$dend), pheno.order=sub.pair, subsets.pheno=subsets.pheno, recurrence=1, allosome='merge', event.type='gene') %>%
    PlotVariants( output.file = str_c('summary/mutation_heatmap_type',sub.ext,'.pdf'),
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
    OrgEvents(sample.order=labels(sub.muts.tree$dend), pheno.order=sub.pair, subsets.pheno=subsets.pheno, recurrence=2, allosome='merge', event.type='gene') %>%
    PlotVariants( output.file = str_c('summary/mutation_heatmap_type_recurrent',sub.ext,'.pdf'),
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
    OrgEvents(sample.order=labels(sub.muts.tree$dend), pheno.order=sub.pair, subsets.pheno=subsets.pheno, recurrence=1, allosome='merge', event.type='gene') %>%
    PlotVariants( output.file = str_c('summary/mutation_heatmap_ccf',sub.ext,'.pdf'),
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
    OrgEvents(sample.order=labels(sub.muts.tree$dend), pheno.order=sub.pair, subsets.pheno=subsets.pheno, recurrence=2, allosome='merge', event.type='gene') %>%
    PlotVariants( output.file = str_c('summary/mutation_heatmap_ccf_recurrent',sub.ext,'.pdf'),
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
    OrgEvents(sample.order=labels(sub.cnas.tree$dend), pheno.order=sub.pair, subsets.pheno=subsets.pheno, recurrence=1, allosome='merge', event.type='band') %>%
    PlotVariants( output.file = str_c('summary/cna_heatmap',sub.ext,'.pdf'),
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
    OrgEvents(sample.order=labels(sub.variants.tree$dend), pheno.order=sub.pair, subsets.pheno=subsets.pheno, recurrence=1, allosome='merge', event.type='span') %>%
    PlotVariants( output.file = str_c('summary/variant_heatmap_type',sub.ext,'.pdf'),
                  clonal      = FALSE,
                  pathogenic  = FALSE,
                  ccf         = FALSE,
                  loh         = FALSE,
                  event.type  = 'span',
                  width       = length(unique(.$sample)) + 10,
                  height      = (length(unique(.$span))/3) + 10,
                  text.size   = 30 )


    #---------------------------
    # Fisher's exact CN plotting
    #---------------------------

    samples.a <- subsets[unlist(subset.groups[sub.num,'a'])] %>% unlist
    samples.b <- subsets[unlist(subset.groups[sub.num,'b'])] %>% unlist

    # copy number plotting
    gene.cn.a <- gene.cn[c('gene', 'chrom', 'start', 'end', samples.a)]
    gene.cn.b <- gene.cn[c('gene', 'chrom', 'start', 'end', samples.b)]

    Fisher( plot.type        = 'copy number',
            gene.matrix.a   = gene.cn.a,
            gene.matrix.b   = gene.cn.b,
            plot.title.main = sub.name,
            plot.title.a    = sub.pair[1],
            plot.title.b    = sub.pair[2],
            allosome        = allosome,
            targets.file    = targets.file,
            suffix          = '',
            threshold.a     = FALSE,
            threshold.b     = FALSE,
            gene.names      = FALSE )

    #---------------------------------
    # Fisher's exact mutation plotting
    #---------------------------------

    # mutation plotting
    gene.muts.a <- sub.muts %>% select(sample, gene, chrom, pos, effect) %>% mutate(effect=1) %>% filter(sample %in% samples.a) %>% spread(sample, effect, fill=0)
    gene.muts.b <- sub.muts %>% select(sample, gene, chrom, pos, effect) %>% mutate(effect=1) %>% filter(sample %in% samples.b) %>% spread(sample, effect, fill=0)

    Fisher( plot.type       = 'mutation',
            gene.matrix.a   = gene.muts.a,
            gene.matrix.b   = gene.muts.b,
            plot.title.main = sub.name,
            plot.title.a    = sub.pair[1],
            plot.title.b    = sub.pair[2],
            allosome        = allosome,
            targets.file    = targets.file,
            suffix          = '',
            gene.names      = FALSE )

    #----------------
    # TREE GENERATION
    #----------------

    system("mkdir summary/tree &>/dev/null")

    # remove labels from tree
    if(tree.labels == FALSE){
        dend.tree %<>% set('labels', '')
    }

    # distance tree
    pdf('summary/tree/genomic_distance_tree_ladder.pdf',140,38)
        par(mar = c(10,10,10,10))  # bottom, left, top, right
        layout(matrix(c(1,2),nrow=1), widths=c(10,1))
        plot(muts.tree$, cex.axis=6)
        colored_bars( colors = color.table %>% select(-sample, -gene, -exists, -pheno),
                      dend = muts.tree$,
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

}






