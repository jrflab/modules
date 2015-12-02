# vim: set ft=make:
ifndef GRCM38_INC
REF = GRCm38

REF_FASTA := $(HOME)/share/reference/Mus_musculus_GRCm38/Mus_musculus.GRCm38.71.dna.chromosome.genome.fa
REF_DICT := $(HOME)/share/reference/Mus_musculus_GRCm38/Mus_musculus.GRCm38.71.dna.chromosome.genome.dict
SNP_EFF_GENOME = GRCm38.79 # == mm10
DBSNP := $(HOME)/share/reference/mgp.v5.merged.snps_all.dbSNP142.vcf.gz
MGP_SNP_DBSNP := $(DBSNP)
MGP_INDEL_DBSNP := $(HOME)/share/reference/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz

EXOME_BED = $(HOME)/share/reference/mus_musculus_known_genes_exons_GRCm38.bed
EXOME_BED_NOHEADER = $(HOME)/share/reference/mus_musculus_known_genes_exons_GRCm38_noheader.bed

FREEC_REF := $(HOME)/share/reference/Mus_musculus_GRCm38/freec
CHR_LEN = $(HOME)/share/reference/Mus_musculus_GRCm38/mm10.len

BOWTIE_REF = $(HOME)/share/reference/Mus_musculus_GRCm38/Mus_musculus.GRCm38.71.dna.chromosome.genome

ENSEMBL_TXDB = $(HOME)/share/reference/mus_musculus_ensembl_biomart.2014-04-28.sqlite

#TXDB = $(HOME)/ensmusg70.08032013.sqlite

RIBOSOMAL_INTERVALS = $(HOME)/share/reference/mm10_rrna_intervals.txt
GENE_REF_FLAT = $(HOME)/share/reference/mm10_genes.refFlat.txt

EXOME ?= false
ifeq ($(EXOME),true)
TARGETS_FILE = $(EXOME_BED_NOHEADER)
QUALIMAP_TARGETS_FILE = $(TARGETS_FILE)
endif

INCLUDE_CHR_Y ?= true
ifneq ($(and $(TARGETS_FILE),$(findstring false,$(EXOME))),)
CHROMOSOMES := $(shell grep -v '@' $(TARGETS_FILE) | cut -f1 | sort | uniq)
else
CHROMOSOMES := $(shell seq 1 19) X $(if $(findstring true,$(INCLUDE_CHR_Y)),Y) MT
endif


endif
GRCM38_INC = true
