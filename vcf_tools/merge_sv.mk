include modules/Makefile.inc

LOGDIR ?= log/merge_sv.$(NOW)

SV_CALLERS = svaba manta gridss
MAX_DIST = 500
NUM_CALLERS = 2
TYPE = 1
STRAND = 1
MIN_SIZE = 30

merge_sv :  $(foreach pair,$(SAMPLE_PAIRS),merge_sv/$(pair)/samples.txt) \
	    $(foreach pair,$(SAMPLE_PAIRS),merge_sv/$(pair)/$(pair).merged_sv.vcf) \
	    $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).merged_sv.vcf)
	   
define merge-sv
merge_sv/$1_$2/samples.txt : $(foreach caller,$(SV_CALLERS),vcf/$1_$2.$(caller)_sv.vcf)
	mkdir -p merge_sv/$1_$2 && \
	$(foreach caller,$(SV_CALLERS),echo vcf/$1_$2.$(caller)_sv.vcf >> $$(@);)

merge_sv/$1_$2/$1_$2.merged_sv.vcf : merge_sv/$1_$2/samples.txt
	$$(call RUN,-c -n 1 -s 4G -m 8G -v $(SURVIVOR_ENV),"set -o pipefail && \
							    SURVIVOR merge $$(<) \
							    $(MAX_DIST) $(NUM_CALLERS) $(TYPE) $(STRAND) 0 $(MIN_SIZE) $$(@)")

vcf/$1_$2.merged_sv.vcf : merge_sv/$1_$2/$1_$2.merged_sv.vcf
	$$(INIT) cat $$(<) > $$(@)
	
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call merge-sv,$(tumor.$(pair)),$(normal.$(pair)))))

..DUMMY := $(shell mkdir -p version; \
	     $(SURVIVOR_ENV)/bin/SURVIVOR --version &> version/merge_sv.txt;)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: merge_sv
