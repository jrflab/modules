ifndef B37_INC
REF = b37

REF_DIR = $(HOME)/share/lib/ref_files
REF_FASTA ?= $(REF_DIR)/b37/human_g1k_v37.fasta
DMP_FASTA ?= $(REF_DIR)/b37/human_g1k_v37_dmp.fasta
DECOY_FASTA ?= $(REF_DIR)/b37/human_g1k_v37_decoy.fasta
DBSNP ?= $(REF_DIR)/b37/dbsnp_138.b37.vcf
KNOWN_INDELS ?= $(REF_DIR)/b37/Mills_and_1000G_gold_standard.indels.b37.vcf
OMNI ?= $(REF_DIR)/b37/1000G_omni2.5.b37.vcf
HAPMAP ?= $(REF_DIR)/b37/hapmap_3.3.b37.vcf
HAPMAP_POP_FILE ?= $(REF_DIR)/b37/hapmap_3.3_b37_pop_stratified_af.vcf.gz
CHROMOSOMES := $(shell grep -v '@' $(TARGETS_FILE) | cut -f1 | sort | uniq)
REF_MSI = $(REF_DIR)/b37/human_g1k_v37_msi.list
ENSEMBL = $(HOME)/share/lib/resource_files/Ensemble_b37.RData
STAR_REF = $(REF_DIR)/STAR_b37
KALLISTO_INDEX = $(REF_DIR)/b37/ensembl_v75-0.43.0_kallisto_index
BISMARK_GENOME = $(REF_DIR)/b37/bismark_bt2/
BWAMETH_GENOME = $(REF_DIR)/b37/bwa_meth/human_g1k_v37.fasta
SALMON_INDEX = $(REF_DIR)/b37/salmon_idx
REF_MSI = $(REF_DIR)/b37/human_g1k_v37.msi

endif
B37_INC = true
