# run star fusion on fastqs
include modules/Makefile.inc

LOGDIR = log/star_fusion.$(NOW)

$(if $(STAR_CTAT_DIR),,$(error no STAR CTAT dir))
STAR_FUSION = $(HOME)/share/usr/STAR-Fusion/STAR-Fusion
STAR_FUSION_OPTS = --genome_lib_dir $(STAR_CTAT_DIR)

PHONY += star_fusion
star_fusion : $(foreach sample,$(SAMPLES),star_fusion/$(sample).star_fusion_timestamp)
	
star_fusion/%.star_fusion_timestamp : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	$(call LSCRIPT_MEM,8G,12G,"$(STAR_FUSION) \
		--genome_lib_dir $(STAR_CTAT_DIR) \
		-J $< --verbose_level 2 \
		--left_fq $^ --right_fq $(<<) && touch $@")

.PHONY: $(PHONY)
.SECONDARY: 
.DELETE_ON_ERROR:

include modules/fastq_tools/mergeSplitFastq.mk
