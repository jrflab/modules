include modules/Makefile.inc

LOGDIR ?= log/star_fish.$(NOW)

MIN_SIZE = 1
MAX_SIZE = 10000000000000000000

star_fish :  $(foreach pair,$(SAMPLE_PAIRS),star_fish/$(pair)/$(pair).merged_sv.bed) \
	     $(foreach pair,$(SAMPLE_PAIRS),star_fish/$(pair)/$(pair).merged_sv.bedpe) \
	     $(foreach pair,$(SAMPLE_PAIRS),star_fish/$(pair)/$(pair).merged_cn.txt) \
	     $(foreach pair,$(SAMPLE_PAIRS),star_fish/$(pair)/$(pair).taskcomplete) \
	     star_fish/summary/taskcomplete
		
define starfish-sv
star_fish/$1_$2/$1_$2.merged_sv.bed : vcf/$1_$2.merged_sv.vcf
	$$(call RUN,-c -n 1 -s 4G -m 8G -v $(SURVIVOR_ENV),"set -o pipefail && \
							    SURVIVOR vcftobed \
							    $$(<) \
							    $(MIN_SIZE) \
							    $(MAX_SIZE) \
							    $$(@)")
							    
star_fish/$1_$2/$1_$2.merged_sv.bedpe : star_fish/$1_$2/$1_$2.merged_sv.bed
	$$(call RUN,-c -n 1 -s 4G -m 8G -v $(STARFISH_ENV),"set -o pipefail && \
							    $(RSCRIPT) $(SCRIPTS_DIR)/star_fish.R \
							    --option 1 \
							    --sample_name $1_$2 \
							    --input_file $$(<) \
							    --output_file $$(@)")
							    
star_fish/$1_$2/$1_$2.merged_cn.txt : facets/cncf/$1_$2.txt
	$$(call RUN,-c -n 1 -s 4G -m 8G -v $(STARFISH_ENV),"set -o pipefail && \
							    $(RSCRIPT) $(SCRIPTS_DIR)/star_fish.R \
							    --option 2 \
							    --sample_name $1_$2 \
							    --input_file $$(<) \
							    --output_file $$(@)")
							    
star_fish/$1_$2/$1_$2.taskcomplete : star_fish/$1_$2/$1_$2.merged_sv.bedpe star_fish/$1_$2/$1_$2.merged_cn.txt
	$$(call RUN,-c -n 1 -s 4G -m 8G -v $(STARFISH_ENV),"set -o pipefail && \
							    $(RSCRIPT) $(SCRIPTS_DIR)/star_fish.R \
							    --option 3 \
							    --sample_name $1_$2")
							    
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call starfish-sv,$(tumor.$(pair)),$(normal.$(pair)))))
		
star_fish/summary/taskcomplete : $(foreach pair,$(SAMPLE_PAIRS),star_fish/$(pair)/$(pair).merged_sv.bedpe) $(foreach pair,$(SAMPLE_PAIRS),star_fish/$(pair)/$(pair).merged_cn.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G -v $(STARFISH_ENV),"set -o pipefail && \
							     $(RSCRIPT) $(SCRIPTS_DIR)/star_fish.R \
							     --option 4 \
							     --sample_name '$(SAMPLE_PAIRS)'")

..DUMMY := $(shell mkdir -p version; \
	     $(STARFISH_ENV)/bin/R --version &> version/star_fish.txt;)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: star_fish
