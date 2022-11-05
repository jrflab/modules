include modules/Makefile.inc

LOGDIR ?= log/anotate_sv.$(NOW)

SV_CALLERS = svaba manta gridss merged

annotate_sv :  $(foreach pair,$(SAMPLE_PAIRS), \
			$(foreach caller,$(SV_CALLERS),annotate_sv/$(pair)/$(pair).$(caller)_sv.txt))
			
define annotate-sv
annotate_sv/$1/$1.$2_sv.txt : vcf/$1.$2_sv.vcf
	$$(call RUN,-c -n 1 -s 4G -m 8G -v $(SURVIVOR_ENV),"set -o pipefail")

endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach caller,$(SV_CALLER), \
		$(eval $(call annotate-sv,$(pair),$(caller)))))
		
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: annotate_sv
