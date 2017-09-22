# vim: set ft=make:

ifndef HG38_INC

REF?= hg38
REF_FASTA ?= $(REF_DIR)/hg38_gatk_bundle/Homo_sapiens_assembly38.fasta

DBSNP ?= $(REF_DIR)/hg38_gatk_bundle/dbsnp_146.hg38.vcf.gz
HAPMAP ?= $(REF_DIR)/hg38_gatk_bundle/hapmap_3.3_grch38_pop_stratified_af.vcf.gz
KNOWN_INDELS ?= $(REF_DIR)/hg38_gatk_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
OMNI ?= $(REF_DIR)/hg38_gatk_bundle/1000G_omni2.5.hg38.vcf.gz

TOPHAT_REF = $(REF)

ANNOVAR_REF = hg38

ERICSCRIPT_DB = $(REF_DIR)/ericscript_db_hg38_84
ERICSCRIPT_SPECIES = homo_sapiens

endif
HG38_INC = true
