# vim: set ft=make:

ifndef HG19_INC
REF = hg19
REF_DIR = $(HOME)/share/reference

TOPHAT_REF = GRCh37
GENES_GTF = $(REF_DIR)/Homo_sapiens/Ensembl/$(TOPHAT_REF)/Annotation/Genes/genes.gtf

#UCSC_REFSEQ_GTF = $(REF_DIR)s/genomes/ucsc.refseq.gtf
#REFSEQ_GTF := $(UCSC_REFSEQ_GTF)

FREEC_REF = $(REF_DIR)/GATK_bundle/2.3/FREEC_ref
CHR_LEN = $(REF_DIR)/GATK_bundle/2.3/hg19.len

SNP_TXT = $(REF_DIR)/hg19_snp137.SingleDiNucl.1based.nochr.txt

GEM_MAP_FILE = $(REF_DIR)/gem/out76_hg19.gem

CHIMERASCAN_REF = $(REF_DIR)/chimerascan_index

RIBOSOMAL_INTERVALS = $(REF_DIR)/rrna_intervals.txt
GENE_REF_FLAT = $(REF_DIR)/hg19_genes.refFlat.txt

CENTROMERE_TABLE = $(REF_DIR)/centromere.txt
CENTROMERE_TABLE2 = $(REF_DIR)/centromere2.txt

MUT_ASS_RDATA = $(REF_DIR)/MA.hg19/MA.Rdata

ENSEMBL_TXDB = $(REF_DIR)/hsapiens_ensembl_biomart.sqlite

MAPSPLICE_REF_DIR = $(REF_DIR)/GATK_bundle/2.3/mapsplice_ref/chr
MAPSPLICE_REF_BASENAME = $(REF_DIR)/GATK_bundle/2.3/mapsplice_ref/hg19

GISTIC_REF ?= $(HOME)/share/usr/gistic_2_0_21/refgenefiles/hg19.mat

ifndef UCSC_REF # regular gatk bundle
   REF_FASTA ?= $(REF_DIR)/GATK_bundle/2.3/human_g1k_v37.fasta
   REF_2BIT ?= $(REF_DIR)/GATK_bundle/2.3/human_g1k_v37.2bit
   DBSNP ?= $(REF_DIR)/GATK_bundle/2.3/dbsnp_137.b37.vcf
   DBSNP1PC ?= $(REF_DIR)/GATK_bundle/2.3/dbsnp_137.b37.G5.GMAF1pc.include_cosmic.vcf
   KNOWN_INDELS ?= $(REF_DIR)/GATK_bundle/2.3/Mills_and_1000G_gold_standard.indels.b37.vcf
   OMNI ?= $(REF_DIR)/GATK_bundle/2.3/1000G_omni2.5.b37.vcf
   HAPMAP ?= $(REF_DIR)/GATK_bundle/2.3/hapmap_3.3.b37.vcf
   COSMIC ?= $(REF_DIR)/CosmicCodingMuts_v69.vcf
   BOWTIE_REF = $(REF_DIR)/GATK_bundle/2.3/human_g1k_v37
   EXOME_BED ?= $(REF_DIR)/SureSelect_50MB_S02972011_Regions_nochr_noAnn.bed
   EXOME_BED_NOHEADER ?= $(REF_DIR)/SureSelect_50MB_S02972011_Regions_nochr_noAnn_noheader.bed
   MAP_BIGWIG ?= $(REF_DIR)/wgEncodeCrgMapabilityAlign100mer.nochr.bigWig
else # ucsc reference
   REF_FASTA ?= $(REF_DIR)/ucsc_gatk_bundle_2.8/ucsc.hg19.fasta
   DBSNP ?= $(REF_DIR)/ucsc_gatk_bundle_2.8/dbsnp_138.hg19.vcf
   DBSNP1PC ?= $(REF_DIR)/ucsc_gatk_bundle_2.8/dbsnp_138.hg19.vcf
   KNOWN_INDELS ?= $(REF_DIR)/ucsc_gatk_bundle_2.8/Mills_and_1000G_gold_standard.indels.hg19.vcf
   OMNI ?= $(REF_DIR)/ucsc_gatk_bundle_2.8/1000G_omni2.5.hg19.vcf
   HAPMAP ?= $(REF_DIR)/ucsc_gatk_bundle_2.8/hapmap_3.3.hg19.vcf
   EXOME_BED ?= $(REF_DIR)/SureSelect_50MB_S02972011_Regions.bed
   EXOME_BED_NOHEADER ?= $(REF_DIR)/SureSelect_50MB_S02972011_Regions.bed
   MAP_BIGWIG ?= $(REF_DIR)/wgEncodeCrgMapabilityAlign100mer.bigWig
endif

EXOME ?= false
ifeq ($(EXOME),true)
TARGETS_FILE = $(EXOME_BED_NOHEADER)
QUALIMAP_TARGETS_FILE = $(REF_DIR)/SureSelect_50MB_S02972011_Regions_nochr.bed
EMU_TARGETS_FILE = $(REF_DIR)/SureSelect_50MB_S02972011_Regions_nochr_noAnn_noheader_numchr_spacedelim.bed
endif

INCLUDE_CHR_Y ?= false
ifdef TARGETS_FILE
CHROMOSOMES := $(shell grep -v '@' $(TARGETS_FILE) | cut -f1 | sort | uniq)
else
CHRS = $(shell seq 1 22) $(if $(findstring true,$(INCLUDE_CHR_Y)),Y) $(if $(UCSC_REF),M,MT)
CHROMOSOMES ?= $(strip $(if $(UCSC_REF),$(addprefix chr,$(CHRS)),$(CHRS)))
endif

KNOWN_CNVS = $(REF_DIR)/GRCh37_hg19_variants_2013-07-23.txt


REF_NAME := 1000GenomesPilot-NCBI37
REF_DATE := 2010517
REF_DICT = $(REF_FASTA:.fasta=.dict)
SNP_EFF_GENOME = hg19

EMU_REF_DIR ?= $(REF_DIR)/GATK_bundle/2.3/emu_ref
EMU_REF_MUTATIONS ?= $(REF_DIR)/emu_ref/exome_ref/Breast_clean_somatic_mutations_for_signature_analysis.emu.mutations.txt
EMU_REF_CNV ?= $(REF_DIR)/emu_ref/exome_ref/Breast_clean_somatic_mutations_for_signature_analysis.emu.copynum.txt

HMMCOPY_GC_WIG = $(REF_DIR)/GATK_bundle/2.3/gc_hg19.wig
HMMCOPY_MAP_WIG = $(REF_DIR)/map_hg19.wig

REF_HOUSEKEEPING_GENE_BED = $(HOME)/share/reference/hg19.HouseKeepingGenes.bed

INTEGRATE_ANN = $(REF_DIR)/annot.ucsc.txt
INTEGRATE_BWTS = $(REF_DIR)/GATK_bundle/2.3/bwts
ONCOTATOR_DB ?= $(REF_DIR)/oncotator_v1_ds_Jan262014

endif
HG19_INC = true

