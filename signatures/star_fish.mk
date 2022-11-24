include modules/Makefile.inc

LOGDIR ?= log/star_fish.$(NOW)

MIN_SIZE = 1
MAX_SIZE = 10000000000000000

star_fish :  $(foreach pair,$(SAMPLE_PAIRS),star_fish/$(pair)/$(pair).merged.bed)
		
define starfish-sv
star_fish/$1_$2/$1_$2.merged.bed : vcf/$1_$2.merged_sv.vcf
	$$(call RUN,-c -n 1 -s 4G -m 8G -v $(SURVIVOR_ENV),"set -o pipefail && \
							    SURVIVOR vcftobed \
							    $$(<) \
							    $(MIN_SIZE) \
							    $(MAX_SIZE) \
							    $$(@)")
							    
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call starfish-sv,$(tumor.$(pair)),$(normal.$(pair)))))

..DUMMY := $(shell mkdir -p version; \
	     $(STARFISH_ENV)/bin/R --version &> version/star_fish.txt;)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: star_fish
