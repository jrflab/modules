include modules/Makefile.inc

LOGDIR ?= log/star_fish.$(NOW)

MIN_SIZE = 1
MAX_SIZE = 10000000000000000

star_fish :  $(foreach pair,$(SAMPLE_PAIRS),star_fish/$(pair)/$(pair).merged.bed) \
	     $(foreach pair,$(SAMPLE_PAIRS),star_fish/$(pair)/$(pair).merged.bedpe) \
	     $(foreach pair,$(SAMPLE_PAIRS),star_fish/$(pair)/$(pair).merged.txt)
		
define starfish-sv
star_fish/$1_$2/$1_$2.merged.bed : vcf/$1_$2.merged_sv.vcf
	$$(call RUN,-c -n 1 -s 4G -m 8G -v $(SURVIVOR_ENV),"set -o pipefail && \
							    SURVIVOR vcftobed \
							    $$(<) \
							    $(MIN_SIZE) \
							    $(MAX_SIZE) \
							    $$(@)")
							    
star_fish/$1_$2/$1_$2.merged.bedpe : star_fish/$1_$2/$1_$2.merged.bed
	$$(call RUN,-c -n 1 -s 4G -m 8G -v $(STARFISH_ENV),"set -o pipefail && \
							    $(RSCRIPT) $(SCRIPTS_DIR)/star_fish.R \
							    --option 1 \
							    --sample_name $1_$2 \
							    --input_file $$(<) \
							    --output_file $$(@)")
							    
star_fish/$1_$2/$1_$2.merged.txt : facets/cncf/$1_$2.txt
	$$(call RUN,-c -n 1 -s 4G -m 8G -v $(STARFISH_ENV),"set -o pipefail && \
							    $(RSCRIPT) $(SCRIPTS_DIR)/star_fish.R \
							    --option 2 \
							    --sample_name $1_$2 \
							    --input_file $$(<) \
							    --output_file $$(@)")
							    
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call starfish-sv,$(tumor.$(pair)),$(normal.$(pair)))))

..DUMMY := $(shell mkdir -p version; \
	     $(STARFISH_ENV)/bin/R --version &> version/star_fish.txt;)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: star_fish
