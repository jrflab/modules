# vim: set ft=make:
REF = mm10

REF_FASTA := $(HOME)/share/reference/Mus_musculus_GRCm38/Mus_musculus.GRCm38.71.dna.chromosome.genome.fa
REF_DICT := $(HOME)/share/reference/Mus_musculus_GRCm38/Mus_musculus.GRCm38.71.dna.chromosome.genome.dict
SNP_EFF_GENOME = GRCm38.70 # == mm10
DBSNP := $(HOME)/share/reference/mus_musculus_dbSNP_GRCm38/dbsnp_138.vcf

EXOME_BED = $(HOME)/share/reference/mus_musculus_known_genes_exons_GRCm38.bed

FREEC_REF := $(HOME)/share/reference/Mus_musculus_GRCm38/freec
CHR_LEN = $(HOME)/share/reference/Mus_musculus_GRCm38/mm10.len

#TXDB = $(HOME)/ensmusg70.08032013.sqlite
