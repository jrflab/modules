include modules/Makefile.inc

LOGDIR ?= log/merge_sv.$(NOW)

SURVIVOR_ENV ?= $(HOME)/share/usr/env/survivor-1.0.7
SV_CALLERS = svaba manta

merge_sv :  $(foreach pair,$(SAMPLE_PAIRS),merge_sv/$(pair)/sample_list.txt)
#	    $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).merged_sv.vcf)
#	    $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).merged_indels.vcf) \
#	    $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).merged_candidate_sv.vcf)
	   
define merge-sv
merge_sv/$1_$2/sample_list.txt : vcf/$1_$2.$3_sv.vcf
	echo vcf/$1_$2.$3_sv.vcf >> $(@); done
	
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach caller,$(SV_CALLERS), \
		$(eval $(call merge-sv,$(tumor.$(pair)),$(normal.$(pair)),$(caller)))))
	
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: merge_sv
