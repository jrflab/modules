include modules/Makefile.inc

LOGDIR ?= log/merge_sv.$(NOW)

SURVIVOR_ENV ?= $(HOME)/share/usr/env/survivor-1.0.7
SV_CALLERS = svaba manta

merge_sv :  $(foreach pair,$(SAMPLE_PAIRS),merge_sv/$(pair)/sample_list_sv.txt) \
	    $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).merged_sv.vcf)
#	    $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).merged_indels.vcf) \
#	    $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).merged_candidate_sv.vcf)
	   
define merge-sv
merge_sv/$1_$2/sample_list_sv.txt : $(foreach caller,$(SV_CALLERS),vcf/$1_$2.$(caller)_sv.vcf)
	echo vcf/$1_$2.svaba_sv.vcf > $$(@)
	echo vcf/$1_$2.manta_sv.vcf >> $$(@)
	
vcf/$1_$2.merged_sv.vcf : merge_sv/$1_$2/sample_list_sv.txt
	$$(call RUN,-c -n 1 -s 4G -m 8G -v $(SURVIVOR_ENV),"set -o pipefail && \
							    SURVIVOR merge $$(<) \
							    500 2 1 1 0 30 $$(@)")
	
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call merge-sv,$(tumor.$(pair)),$(normal.$(pair)))))
	
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: merge_sv
