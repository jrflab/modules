ifndef B37_INC
REF = b37
REF_DIR = $(HOME)/share/reference

STAR_REF_DIR = $(REF_DIR)/GATK_bundle/2.3/STAR_genome
STAR_CTAT_DIR= $(REF_DIR)/GRCh37_gencode_v19_CTAT_lib_July272016

TOPHAT_REF = GRCh37
GENES_GTF = $(REF_DIR)/Homo_sapiens/Ensembl/$(TOPHAT_REF)/Annotation/Genes/genes.gtf

ANNOVAR_REF = hg19

CLINVAR = $(REF_DIR)/clinvar_20160531.vcf.gz
EXAC_NONTCGA = $(REF_DIR)/ExAC.r0.3.nonTCGA.sites.vep.2.vcf.gz

#UCSC_REFSEQ_GTF = $(REF_DIR)s/genomes/ucsc.refseq.gtf
#REFSEQ_GTF := $(UCSC_REFSEQ_GTF)

FREEC_REF = $(REF_DIR)/GATK_bundle/2.3/FREEC_ref
CHR_LEN = $(REF_DIR)/GATK_bundle/2.3/b37.len

SNP_TXT = $(REF_DIR)/b37_snp137.SingleDiNucl.1based.nochr.txt

GEM_MAP_FILE = $(REF_DIR)/gem/out76_hg19.gem

CHIMERASCAN_REF = $(REF_DIR)/chimerascan_index

RIBOSOMAL_INTERVALS = $(REF_DIR)/rrna_intervals.txt
GENE_REF_FLAT = $(REF_DIR)/b37_genes.refFlat.txt

CENTROMERE_TABLE = $(REF_DIR)/centromere.txt
CENTROMERE_TABLE2 = $(REF_DIR)/centromere2.txt

MUT_ASS_RDATA = $(REF_DIR)/MA.hg19/MA.Rdata

ENSEMBL_TXDB = $(REF_DIR)/hsapiens_ensembl_biomart.sqlite

MAPSPLICE_REF_DIR = $(REF_DIR)/GATK_bundle/2.3/mapsplice_ref/chr
MAPSPLICE_REF_BASENAME = $(REF_DIR)/GATK_bundle/2.3/mapsplice_ref/hg19

GISTIC_REF ?= $(HOME)/share/usr/gistic_2_0_21/refgenefiles/hg19.mat

# align to combined reference if PDX, filter mouse chromosomes afterwards
PDX ?= false
ifeq ($(PDX),true)
BWAMEM_REF_FASTA ?= $(REF_DIR)/pdx/combined_Mus_musculus_GRCm38_human_g1k_v37.fasta
endif

REF_FASTA ?= $(REF_DIR)/GATK_bundle/2.3/human_g1k_v37.fasta
NOCHR_REF_FASTA ?= $(REF_DIR)/GATK_bundle/2.3/human_g1k_v37.fasta
CHR_REF_FASTA ?= $(REF_DIR)/ucsc_gatk_bundle_2.8/ucsc.hg19.fasta
REF_2BIT ?= $(REF_DIR)/GATK_bundle/2.3/human_g1k_v37.2bit
DBSNP ?= $(REF_DIR)/dbsnp_138.b37.gmaf.vcf.gz
KNOWN_INDELS ?= $(REF_DIR)/GATK_bundle/2.3/Mills_and_1000G_gold_standard.indels.b37.vcf.gz
OMNI ?= $(REF_DIR)/GATK_bundle/2.3/1000G_omni2.5.b37.vcf.gz
HAPMAP ?= $(REF_DIR)/GATK_bundle/2.3/hapmap_3.3.b37.vcf.gz
COSMIC ?= $(REF_DIR)/cosmic_v82/CosmicCodingMuts.vcf.gz
COSMIC_NONCODING ?= $(REF_DIR)/cosmic_v82/CosmicNonCodingVariants.vcf.gz
HAPMAP_POP_FILE ?= $(REF_DIR)/b37_gatk_bundle_2017-04-26/hapmap_3.3_b37_pop_stratified_af.vcf.gz
BOWTIE_REF = $(REF_DIR)/GATK_bundle/2.3/human_g1k_v37
EXOME_BED ?= $(REF_DIR)/target_panels/sure_select_exome_v4_b37.bed
MAP_BIGWIG ?= $(REF_DIR)/wgEncodeCrgMapabilityAlign100mer.nochr.bigWig

EXOME ?= false
ifeq ($(EXOME),true)
TARGETS_FILE = $(EXOME_BED)
QUALIMAP_TARGETS_FILE = $(REF_DIR)/SureSelect_50MB_S02972011_Regions_nochr.bed
EMU_TARGETS_FILE = $(REF_DIR)/SureSelect_50MB_S02972011_Regions_nochr_noAnn_noheader_numchr_spacedelim.bed
endif

INCLUDE_CHR_Y ?= false
ifneq ($(and $(TARGETS_FILE),$(findstring false,$(EXOME))),)
CHROMOSOMES := $(shell grep -v '@' $(TARGETS_FILE) | cut -f1 | sort | uniq)
else
CHROMOSOMES ?= $(shell seq 1 22) X $(if $(findstring true,$(INCLUDE_CHR_Y)),Y) MT
endif

KNOWN_CNVS = $(REF_DIR)/GRCh37_hg19_variants_2013-07-23.txt


REF_NAME := 1000GenomesPilot-NCBI37
REF_DATE := 20100517
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

ANNOVAR_DB = $(HOME)/share/usr/annovar/humandb

HISAT_REF = $(HOME)/share/reference/hisat/grch37_snp_tran/genome_snp_tran

HOTSPOT_VCF_DIR = modules/reference/hotspots
HOTSPOT_VCF.ext = $(HOTSPOT_VCF_DIR)/hotspot-v3.vcf.gz
HOTSPOT_VCF.int = $(HOTSPOT_VCF_DIR)/hotspot-v1.vcf.gz

HAPLOTYPE_INSUF_BED ?= modules/reference/gene_lists/haplo_insuff_genes.bed
CANCER_GENE_CENSUS_BED ?= modules/reference/gene_lists/cancer_gene_census.b37.2017-05-25.bed
KANDOTH_BED ?= modules/reference/gene_lists/Kandoth_127genes.bed
LAWRENCE_BED ?= modules/reference/gene_lists/Lawrence_cancer5000-S.bed

CN_ENDOMETRIAL_SUBTYPES = CN_high CN_low Endometrioid_MSI_H Endometrioid_MSS Endometrioid MSI POLE Serous
CN_BREAST_SUBTYPES = ER_negative ER_positive HER2_postitive Pam50_Basal Pam50_Her2 Pam50_LumA Pam50_LumB Pam50_Normal Triple_negative
CN_ENDOMETRIAL_BED = $(foreach set,$(CN_ENDOMETRIAL_SUBTYPES), $(HOME)/share/reference/annotation_gene_lists/cn_reg/endometrial/copy_number_regulated_genes_subtype_$(set)_spearmanrsquare0.4_fdrbh_adjp_lt0.05.HUGO.bed)
CN_BREAST_BED = $(foreach set,$(CN_BREAST_SUBTYPES), $(HOME)/share/reference/annotation_gene_lists/cn_reg/breast/metabric_subtype_$(set)_copy_number_regulated_genes_std0.5_spearmanrsquare0.4_fdrbh_adjp_lt0.05.HUGO.bed)

PARSSNP_RESOURCES = $(REF_DIR)/ParsSNP_resources.Rdata

REF_MSI = $(REF_DIR)/GATK_bundle/2.3/b37_msi.list

ONCOTATOR_REF = hg19
ONCOTATOR_DB = $(REF_DIR)/oncotator_v1_ds_April052016/
ONCOTATOR_TX_OVERRIDES = $(REF_DIR)/tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt

CMO_HOTSPOT_FILE = $(REF_DIR)/hotspot_list.txt

VEP_REF = GRCh37
VEP_DATA = $(REF_DIR)/vep/v86
VEP_OVERRIDES = $(REF_DIR)/isoform_overrides_at_mskcc
VEP_SPECIES = homo_sapiens


FUENTES_BED = $(REF_DIR)/fuentes_blacklist.b37.bed
DGD_BED = $(REF_DIR)/dgd.b37.bed

ONCOKB = $(REF_DIR)/allActionableVariants.txt

GNOMAD_DB_DIR = $(REF_DIR)/gnomad_data/vcf/genomes/
GNOMAD_PREFIX = gnomad.genomes.r2.0.1.sites
ALL_GNOMAD_CHROMOSOMES = 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X
GNOMAD_CHROMOSOMES = $(shell comm -12 <(echo "$(CHROMOSOMES)" | tr " " "\n" | sort) <(echo "$(ALL_GNOMAD_CHROMOSOMES)" | tr " " "\n" | sort))

GENE_DB = $(REF_DIR)/genes_homo_sapiens.GRCh37.75.db

ERICSCRIPT_DB = $(REF_DIR)/ericscript_db_hg38_84

endif
B37_INC = true

