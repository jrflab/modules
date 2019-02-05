include modules/Makefile.inc

LOGDIR ?= log/medicc.$(NOW)
PHONY += medicc

medicc : $(foreach set,$(SAMPLE_SETS),medicc/$(set)/.tsv)

define 
# gather all Log2 ratio and BAF for a given patient

endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call combine-samples-pdx,$(set))))

define
# segment all Log2 ratio and BAF for a give patient
 
endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call combine-samples-pdx,$(set))))

define
# predict total and parental copy number aberrations
 
endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call combine-samples-pdx,$(set))))
		
define
# initial run of MEDICC
 
endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call combine-samples-pdx,$(set))))

define
# bootstrapped runs of MEDICC
 
endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call combine-samples-pdx,$(set))))
		
		
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

