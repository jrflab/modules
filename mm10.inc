# vim: set ft=make:
ifndef MM10_INC
REF = mm10

REF_FASTA := $(HOME)/share/reference/Mus_musculus_GRCm38/Mus_musculus.GRCm38.71.dna.chromosome.genome.fa
REF_DICT := $(HOME)/share/reference/Mus_musculus_GRCm38/Mus_musculus.GRCm38.71.dna.chromosome.genome.dict
SNP_EFF_GENOME = GRCm38.70 # == mm10
DBSNP := $(HOME)/share/reference/mus_musculus_dbSNP_GRCm38/dbsnp_138.vcf

EXOME_BED = $(HOME)/share/reference/mus_musculus_known_genes_exons_GRCm38.bed

FREEC_REF := $(HOME)/share/reference/Mus_musculus_GRCm38/freec
CHR_LEN = $(HOME)/share/reference/Mus_musculus_GRCm38/mm10.len

BOWTIE_REF = $(HOME)/share/reference/Mus_musculus_GRCm38/Mus_musculus.GRCm38.71.dna.chromosome.genome

ENSEMBL_TXDB = $(HOME)/share/reference/mus_musculus_ensembl_biomart.2014-04-28.sqlite

#TXDB = $(HOME)/ensmusg70.08032013.sqlite

RIBOSOMAL_INTERVALS = $(HOME)/share/reference/mm10_rrna_intervals.txt
GENE_REF_FLAT = $(HOME)/share/reference/mm10_genes.refFlat.txt

CHR_CHROMOSOMES ?= false
ifdef TARGETS_FILE
   CHROMOSOMES := $(shell grep -v '@' $(TARGETS_FILE) | cut -f1 | sort | uniq)
else
   ifeq ($(CHR_CHROMOSOMES),true)
      CHROMOSOMES := chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrM
      ifeq ($(INCLUDE_CHR_Y),true)
	    CHROMOSOMES += chrY
      endif
   else
      CHROMOSOMES := 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X MT
      ifeq ($(INCLUDE_CHR_Y),true)
	    CHROMOSOMES += Y
      endif
   endif
endif
INCLUDE_CHR_Y ?= true

endif
MM10_INC = true
