# star alignment + star fusion
#
include modules/Makefile.inc

STAR_CHIMERIC = true

STAR_FUSION = $(HOME)/share/usr/STAR-Fusion/STAR-Fusion

$(if $(STAR_CTAT_DIR),,$(error no STAR CTAT dir))

PHONY += star_fusion
star_fusion : $(foreach sample,$(SAMPLES),star_fusion/$(sample).star_fusion_timestamp)

star_fusion/%.star_fusion_timestamp : star/%.Chimeric.out.junction
	$(call RUN,-s 8G -m 12G,"$(STAR_FUSION) --genome_lib_dir $(STAR_CTAT_DIR) -J $< --output_dir $(@D)/$*.star_fusion && touch $@")

include modules/aligners/starAligner.mk
