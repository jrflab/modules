include modules/Makefile.inc

LOGDIR ?= log/medicc.$(NOW)
PHONY += medicc medicc/mad medicc/mpcf medicc/ascat

medicc : $(foreach set,$(SAMPLE_SETS),medicc/mad/$(set).RData) $(foreach set,$(SAMPLE_SETS),medicc/mpcf/$(set).RData)

define combine-samples
medicc/mad/%.RData : $(wildcard $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).Rdata))
	$$(call RUN,-c -s 8G -m 12G,"$(RSCRIPT) modules/test/phylogeny/combinesamples.R --sample_set $$* --normal_samples '$(NORMAL_SAMPLES)'")

endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call combine-samples,$(set))))

define ascat-mpcf
medicc/mpcf/%.RData : medicc/mad/%.RData
	$$(call RUN,-c -s 8G -m 12G -v ~/share/usr/anaconda-envs/ascat,"$(RSCRIPT) modules/test/phylogeny/segmentsamples.R --sample_set $$* --normal_samples '$(NORMAL_SAMPLES)'")

endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call ascat-mpcf,$(set))))

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

