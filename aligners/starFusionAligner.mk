include modules/Makefile.inc

STAR_CHIMERIC = true

STAR_FUSION = STAR-Fusion

STAR_FUSION_ENV = $(HOME)/share/usr/anaconda-envs/star-fusion-1.0.0

STAR_FUSION_TO_USV = python modules/sv_callers/starfusion2usv.py

$(if $(STAR_CTAT_DIR),,$(error no STAR CTAT dir))

PHONY += star_fusion
star_fusion : $(foreach sample,$(SAMPLES),star_fusion/$(sample).star_fusion_timestamp)

star_fusion/%.star_fusion_timestamp : star/%.Chimeric.out.junction
	$(call RUN,-v $(STAR_FUSION_ENV) -s 8G -m 12G,"$(STAR_FUSION) --genome_lib_dir $(STAR_CTAT_DIR) -J $< --output_dir $(@D)/$* && touch $@")

usv/%.star_fusion.tsv : star_fusion/%.star_fusion_timestamp
	$(call RUN,,"$(STAR_FUSION_TO_USV) < $(<D)/$*/star-fusion.fusion_candidates.final > $@")

include modules/aligners/starAligner.mk
