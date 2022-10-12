#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("magrittr"))
suppressPackageStartupMessages(library("VariantAnnotation"))
suppressPackageStartupMessages(library("Biostrings"))
suppressPackageStartupMessages(library("GenomicRanges"))
suppressPackageStartupMessages(library("BSgenome.Hsapiens.UCSC.hg19"))
                               
if (!interactive()) {
    options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))
}

optList = list (make_option("--file_in", default = NULL, type = "character", help = "input vcf file to annotate"),
		make_option("--file_out", default = NULL, type= "character", help = "output vcf file to write"),
                make_option("--ensembl_gene", default = NULL, type= "character", help = "path to .RData workspace of Ensembl genes"))
parser = OptionParser(usage = "%prog [options]", option_list = optList)
arguments = parse_args(parser, positional_arguments = T)
opt = arguments$options

ref_genome = BSgenome.Hsapiens.UCSC.hg19
load(as.character(opt$ensembl_gene))
if (!any(grepl("chr", ensgene$Chromosome))) {
	ensgene$Chromosome = paste0("chr", ensgene$Chromosome)
}

'add_mutation_context' <- function (vr = NULL, ref_genome = NULL, k = 3, check_strand = FALSE, num_of_dels = 0)
{
    mid = (k + 1)/2
    gr = granges(vr)
    strand_mut = strand(gr)
    ranges = GenomicRanges::resize(gr, k, fix = "center")
    context = Biostrings::getSeq(ref_genome, ranges)
    ref_base = DNAStringSet(VariantAnnotation::ref(vr))
    alt_base = DNAStringSet(VariantAnnotation::alt(vr))
    ref0 = subseq(context, mid, mid)
    idx_invalid = (ref0 != ref_base)
    wronguns = sum(idx_invalid) - num_of_dels
    if (check_strand) {
        idx_minus = (strand_mut == "-")
        context[idx_minus] = Biostrings::reverseComplement(context[idx_minus])
        ref_base[idx_complement] = Biostrings::reverseComplement(ref_base[idx_complement])
        alt_base[idx_complement] = Biostrings::reverseComplement(alt_base[idx_complement])
        strand_mut[idx_minus] = "+"
    }
    idx_complement = as.character(ref_base) %in% c("A", "G")
    context[idx_complement] = Biostrings::reverseComplement(context[idx_complement])
    ref_base[idx_complement] = Biostrings::reverseComplement(ref_base[idx_complement])
    alt_base[idx_complement] = Biostrings::reverseComplement(alt_base[idx_complement])
    strand_mut[idx_complement] = "-"
    alteration = as.character(xscat(ref_base, alt_base))
    alteration[idx_invalid] = NA
    context = as.character(context)
    context[idx_invalid] = NA
    vr$alteration = alteration
    vr$context = context
    vr@strand = Rle(strand_mut)
    return(invisible(vr))
}

vcf = readr::read_tsv(file=as.character(opt$file_in), col_names = TRUE, comment = "##", col_types = cols(.default = col_character())) %>%
      readr::type_convert() %>%
      dplyr::select(CHROM = `#CHROM`,
                    POS = POS,
                    REF = REF,
                    ALT = ALT) %>%
      dplyr::filter(nchar(REF) == 1 & nchar(ALT) == 1,
                    REF != "-",
                    ALT != "-") %>%
      dplyr::mutate(CHROM = as.character(gsub("chr", "", CHROM, fixed=TRUE))) %>%
      dplyr::filter(CHROM != "M",
                    CHROM != "MT") %>%
      dplyr::mutate(CHROM = ifelse(CHROM=="X", "23", CHROM)) %>%
      dplyr::mutate(CHROM = ifelse(CHROM=="Y", "24", CHROM)) %>%
      dplyr::mutate(CHROM = as.numeric(CHROM)) %>%
      dplyr::arrange(CHROM, POS, REF, ALT) %>%
      dplyr::mutate(CHROM = ifelse(CHROM==23, "X", CHROM)) %>%
      dplyr::mutate(CHROM = ifelse(CHROM=="24", "Y", CHROM)) %>%
      dplyr::mutate(CHROM = paste0("chr", CHROM)) %>%
      dplyr::mutate(strand_mut = "+",
                    strand_gene = NA,
                    gene_symbol = NA)
ensgene_split = split(ensgene, ensgene %>% .[["Chromosome"]])
vcf_split = split(vcf, vcf %>% .[["CHROM"]])
chr_listy = gtools::mixedsort(intersect(names(ensgene_split),names(vcf_split)))
total = length(chr_listy)
for(chr in chr_listy) {
    ind = unlist(sapply(vcf_split[[chr]]$POS, function(pos) {
                tmp = which(ensgene_split[[chr]]$`Gene Start` <= pos & ensgene_split[[chr]]$`Gene End` >= pos)
                if (length(tmp)==1) {
                    res = tmp
                } else {
                    res = NA
                }
                return(res)
            }))
    if (all(is.na(ind))) {
        vcf_split[[chr]]$gene_symbol = NA 
    } else {
        vcf_split[[chr]]$gene_symbol = ensgene_split[[chr]][ind,"Associated Gene Symbol"]
    }
}
vcf = do.call(rbind, vcf_split) %>%
      dplyr::mutate(strand_gene = c("-", NA, "+")[ensgene[match(gene_symbol, ensgene$`Associated Gene Symbol`), "Strand"] + 2]) %>%
      dplyr::mutate(strand_gene = ifelse(is.na(strand_gene) | strand_gene == "*", "+", strand_gene))
vr = VRanges(seqnames = vcf %>% .[["CHROM"]],
             ranges = IRanges(start = vcf %>% .[["POS"]], end = vcf %>% .[["POS"]]),
             ref = vcf %>% .[["REF"]],
             alt = vcf %>% .[["ALT"]])
vr@strand = Rle(strand(vcf$strand_mut))
vr3 = add_mutation_context(vr = vr, ref = ref_genome, k = 3, num_of_dels = 0)
vr5 = add_mutation_context(vr = vr, ref = ref_genome, k = 5, num_of_dels = 0)
vcf = vcf %>%
      dplyr::mutate(strand_mut = as.character(vr3@strand),
                    strand_ts = NA,
                    strand_ts = ifelse(strand_gene == "+" & strand_mut == "-", "ts", strand_ts),
                    strand_ts = ifelse(strand_gene == "+" & strand_mut == "+", "nt", strand_ts),
                    strand_ts = ifelse(strand_gene == "-" & strand_mut == "+", "ts", strand_ts),
                    strand_ts = ifelse(strand_gene == "-" & strand_mut == "-", "nt", strand_ts))
index = do.call(rbind, vcf_split) %>%
        dplyr::mutate(strand_gene = c("-", NA, "+")[ensgene[match(gene_symbol, ensgene$`Associated Gene Symbol`), "Strand"] + 2]) %>%
        dplyr::mutate(index = is.na(strand_gene) | strand_gene == "*", "+") %>%
        .[["index"]]
vcf = vcf %>%
      dplyr::mutate(strand_gene = ifelse(index, NA, strand_gene),
                    strand_ts = ifelse(index, NA, strand_ts),
                    substype = as.character(vr3$alteration),
                    context3 = as.character(vr3$context),
                    sbs_cat3 = paste0(substype, "_", context3),
                    context5 = as.character(vr5$context),
                    sbs_cat5 = paste0(substype, "_", context5)) %>%
      dplyr::select(chrom = CHROM,
                    pos = POS,
                    ref = REF,
                    alt = ALT,
                    gene_symbol,
                    strand_mut,
                    strand_gene,
                    strand_ts,
                    substype,
                    context3,
                    sbs_cat3,
                    context5,
                    sbs_cat5) %>%
    dplyr::arrange(pos, ref, alt) %>%
    dplyr::mutate(chrom = gsub("chr", "", x=chrom))
index = gtools::mixedorder(vcf$chrom)
vcf = vcf[index,,drop=FALSE]
write_tsv(vcf, path=as.character(opt$file_out), append = FALSE, col_names = TRUE)
