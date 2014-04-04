# vim: set ft=make:
REF = hg19
GSNAP_REF_DIR = $(HOME)/share/references/gmapdb
#UCSC_REF_FASTA := $(HOME)/share/references/genomes/ucsc.hg19.fasta
#UCSC_REF_DICT = $(HOME)/share/references/genomes/ucsc.hg19.dict
REF_FASTA ?= $(HOME)/share/reference/GATK_bundle/2.3/human_g1k_v37.fasta
REF_NAME := 1000GenomesPilot-NCBI37
REF_DATE := 2010517
#REF_CHR_ORD = $(HOME)/share/references/genomes/ucsc.hg19.chr_order
#REF_FA = $(HOME)/share/references/genomes/hg19.fa
REF_DICT = $(REF_FASTA:.fasta=.dict)
SNP_EFF_GENOME = GRCh37.70 # == hg19
DBSNP := $(HOME)/share/reference/GATK_bundle/2.3/dbsnp_137.b37.vcf
DBSNP1PC := $(HOME)/share/reference/GATK_bundle/2.3/dbsnp_137.b37.G5.GMAF1pc.vcf
KNOWN_INDELS := $(HOME)/share/reference/GATK_bundle/2.3/Mills_and_1000G_gold_standard.indels.b37.vcf
OMNI := $(HOME)/share/reference/GATK_bundle/2.3/1000G_omni2.5.b37.vcf
HAPMAP := $(HOME)/share/reference/GATK_bundle/2.3/hapmap_3.3.b37.vcf
COSMIC := $(HOME)/share/reference/CosmicCodingMuts_v64_02042013_noLimit.vcf
#TXDB = $(HOME)/ensg69.biomart.13012013.sqlite

#GC_WIG = $(HOME)/share/references/gc_hg19.wig
#MAP_WIG = $(HOME)/share/references/map_hgi19.wig

#BOWTIE_UCSC_REF = $(HOME)/share/references/genomes/ucsc.hg19
BOWTIE_REF = $(HOME)/share/reference/GATK_bundle/2.3/human_g1k_v37

#CUFFLINKS_REF = GRCh37
#CUFFLINKS_REF_GTF = $(HOME)/share/references/iGenomes/Homo_sapiens/Ensembl/${CUFFLINKS_REF}/Annotation/Genes/genes.gtf

#UCSC_REFSEQ_GTF = $(HOME)/share/references/genomes/ucsc.refseq.gtf
#REFSEQ_GTF := $(UCSC_REFSEQ_GTF)

EXOME_BED = $(HOME)/share/reference/SureSelect_50MB_S02972011_Regions_nochr_noAnn.bed
EXOME_BED_NOHEADER = $(HOME)/share/reference/SureSelect_50MB_S02972011_Regions_nochr_noAnn.bed
FREEC_REF = $(HOME)/share/reference/GATK_bundle/2.3/FREEC_ref
CHR_LEN = $(HOME)/share/reference/GATK_bundle/2.3/hg19.len

SNP_TXT = $(HOME)/share/reference/hg19_snp137.SingleDiNucl.1based.nochr.txt

GEM_MAP_FILE = $(HOME)/share/reference/gem/out76_hg19.gem

CHIMERASCAN_REF = $(HOME)/share/reference/chimerascan_index

RIBOSOMAL_INTERVALS = $(HOME)/share/reference/rrna_intervals.txt
GENE_REF_FLAT = $(HOME)/share/reference/hg19_genes.refFlat.txt

CENTROMERE_TABLE = $(HOME)/share/reference/centromere.txt
CENTROMERE_TABLE2 = $(HOME)/share/reference/centromere2.txt

MUT_ASS_RDATA = $(HOME)/share/reference/MA.hg19/MA.Rdata

ENSEMBL_TXDB = $(HOME)/share/reference/ensembl_biomart.sqlite

MAPSPLICE_REF_DIR = $(HOME)/share/reference/GATK_bundle/2.3/mapsplice_ref/chr
MAPSPLICE_REF_BASENAME = $(HOME)/share/reference/GATK_bundle/2.3/mapsplice_ref/human_g1k_v37
