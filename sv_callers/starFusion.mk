# run star fusion on fastqs
include modules/Makefile.inc

LOGDIR = log/star_fusion.$(NOW)

$(if $(STAR_CTAT_DIR),,$(error no STAR CTAT dir))
STAR_FUSION = STAR-Fusion
STAR_FUSION_ENV = $(HOME)/share/usr/anaconda-envs/star-fusion-1.0.0
STAR_FUSION_OPTS = --genome_lib_dir $(STAR_CTAT_DIR)

STAR_FUSION_TO_USV = python modules/sv_callers/starfusion2usv.py


PHONY += star_fusion
star_fusion : $(foreach sample,$(SAMPLES),usv/$(sample).star_fusion.tsv)
	
star_fusion/%.star_fusion_timestamp : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	$(call RUN,-v $(STAR_FUSION_ENV) -n 8 -s 2G -m 5G,"$(STAR_FUSION) \
		--CPU 8 \
		--output_dir $(@D)/$* \
		--genome_lib_dir $(STAR_CTAT_DIR) \
		--verbose_level 2 \
		--left_fq $< --right_fq $(<<) && touch $@")

usv/%.star_fusion.tsv : star_fusion/%.star_fusion_timestamp
	$(call RUN,,"$(STAR_FUSION_TO_USV) < $(<D)/$*/star-fusion.fusion_candidates.final > $@")

.PHONY: $(PHONY)
.SECONDARY: 
.DELETE_ON_ERROR:

include modules/fastq_tools/mergeSplitFastq.mk
