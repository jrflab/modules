include modules/Makefile.inc

LOGDIR ?= log/medicc.$(NOW)
PHONY += medicc medicc/mad

medicc : $(foreach set,$(SAMPLE_SETS),medicc/mad/$(set).RData)

define combine-samples
medicc/mad/%.RData : $(wildcard $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).Rdata))
	$$(call RUN,-c -s 8G -m 12G,"$(RSCRIPT) modules/test/phylogeny/combinesamples.R --sample_set $$* --normal_samples '$(NORMAL_SAMPLES)'")

endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call combine-samples,$(set))))

#define
# segment all Log2 ratio and BAF for a give patient
# 
#endef
#$(foreach set,$(SAMPLE_SETS),\
#		$(eval $(call combine-samples-pdx,$(set))))
#
#define
# predict total and parental copy number aberrations
# 
#endef
#$(foreach set,$(SAMPLE_SETS),\
#		$(eval $(call combine-samples-pdx,$(set))))
#		
#define
# initial run of MEDICC
# 
#endef
#$(foreach set,$(SAMPLE_SETS),\
#		$(eval $(call combine-samples-pdx,$(set))))
#
#define
# bootstrapped runs of MEDICC
# 
#endef
#$(foreach set,$(SAMPLE_SETS),\
#		$(eval $(call combine-samples-pdx,$(set))))
		
		
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

