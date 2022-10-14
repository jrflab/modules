include modules/Makefile.inc

LOGDIR ?= log/merge_sv.$(NOW)

SURVIVOR_ENV ?= $(HOME)/share/usr/env/survivor-1.0.7
SV_CALLERS = svaba manta

merge_sv :  $(foreach pair,$(SAMPLE_PAIRS),merge_sv/$(pair)/sample_list.txt)
#	    $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).merged_sv.vcf)
#	    $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).merged_indels.vcf) \
#	    $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).merged_candidate_sv.vcf)
	   
define merge-sv
merge_sv/$1_$2/sample_list.txt : $(foreach caller,$(SV_CALLERS),vcf/$1_$2.$(caller)_sv.vcf)
	for i in $(SV_CALLERS); do \
		echo vcf/$1_$2.$i_sv.vcf >> $@;
	done
	
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call merge-sv,$(tumor.$(pair)),$(normal.$(pair)))))
	
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: merge_sv
