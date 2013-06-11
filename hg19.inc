# vim: set ft=make:
REF = hg19
GSNAP_REF_DIR = $(HOME)/share/references/gmapdb
#UCSC_REF_FASTA := $(HOME)/share/references/genomes/ucsc.hg19.fasta
#UCSC_REF_DICT = $(HOME)/share/references/genomes/ucsc.hg19.dict
REF_FASTA := $(HOME)/share/reference/GATK_bundle/2.3/human_g1k_v37.fasta
#REF_CHR_ORD = $(HOME)/share/references/genomes/ucsc.hg19.chr_order
#REF_FA = $(HOME)/share/references/genomes/hg19.fa
REF_DICT := $(HOME)/share/reference/GATK_bundle/2.3/human_g1k_v37.dict
SNP_EFF_GENOME = GRCh37.70 # == hg19
DBSNP := $(HOME)/share/reference/GATK_bundle/2.3/dbsnp_137.b37.vcf
KNOWN_INDELS := $(HOME)/share/reference/GATK_bundle/2.3/Mills_and_1000G_gold_standard.indels.b37.vcf
OMNI := $(HOME)/share/reference/GATK_bundle/2.3/1000G_omni2.5.b37.vcf
HAPMAP := $(HOME)/share/reference/GATK_bundle/2.3/hapmap_3.3.b37.vcf
COSMIC := $(HOME)/share/reference/CosmicCodingMuts_v64_02042013_noLimit.vcf
#TXDB = $(HOME)/ensg69.biomart.13012013.sqlite

#GC_WIG = $(HOME)/share/references/gc_hg19.wig
#MAP_WIG = $(HOME)/share/references/map_hgi19.wig

#BOWTIE_UCSC_REF = $(HOME)/share/references/genomes/ucsc.hg19
#BOWTIE_REF = $(BOWTIE_UCSC_REF)

#CUFFLINKS_REF = GRCh37
#CUFFLINKS_REF_GTF = $(HOME)/share/references/iGenomes/Homo_sapiens/Ensembl/${CUFFLINKS_REF}/Annotation/Genes/genes.gtf

#UCSC_REFSEQ_GTF = $(HOME)/share/references/genomes/ucsc.refseq.gtf
#REFSEQ_GTF := $(UCSC_REFSEQ_GTF)
