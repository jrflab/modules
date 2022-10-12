ifndef HG19_INC
REF = hg19

REF_DIR = $(HOME)/share/reference/
REF_FASTA ?= $(REF_DIR)/ucsc_gatk_bundle_2.8/ucsc.hg19.fasta
DBSNP ?= $(REF_DIR)/ucsc_gatk_bundle_2.8/dbsnp_138.hg19.vcf
KNOWN_INDELS ?= $(REF_DIR)/ucsc_gatk_bundle_2.8/Mills_and_1000G_gold_standard.indels.hg19.vcf
OMNI ?= $(REF_DIR)/ucsc_gatk_bundle_2.8/1000G_omni2.5.hg19.vcf
HAPMAP ?= $(REF_DIR)/ucsc_gatk_bundle_2.8/hapmap_3.3.hg19.vcf

endif
HG19_INC = true

