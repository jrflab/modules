include modules/Makefile.inc

LOGDIR ?= log/medicc.$(NOW)

ALLELE_SPECIFIC_COPY ?= false

ifeq ($(ALLELE_SPECIFIC_COPY),true)

PHONY += medicc medicc/allele_specific medicc/allele_specific/mad medicc/allele_specific/ascat medicc/allele_specific/aspcf medicc/allele_specific/medicc

medicc : $(foreach set,$(SAMPLE_SETS),medicc/allele_specific/medicc/$(set)/tree_final.new) $(foreach set,$(SAMPLE_SETS),medicc/allele_specific/medicc/$(set)/tree_final.pdf)

define allele-specific-medicc
medicc/allele_specific/mad/%.RData : $(wildcard $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).Rdata))
	$$(call RUN,-c -s 8G -m 12G -v $(ASCAT_ENV),"mkdir -p medicc/allele_specific && \
												 mkdir -p medicc/allele_specific/mad && \
												 $(RSCRIPT) modules/test/phylogeny/combinesamples.R --sample_set $$* --normal_samples '$(NORMAL_SAMPLES)' --type allele_specific")

medicc/allele_specific/aspcf/%.RData : medicc/allele_specific/mad/%.RData
	$$(call RUN,-c -s 8G -m 12G -v $(ASCAT_ENV),"mkdir -p medicc/allele_specific/ascat && \
												 mkdir -p medicc/allele_specific/aspcf && \
												 $(RSCRIPT) modules/test/phylogeny/segmentsamples.R --sample_set $$* --normal_samples '$(NORMAL_SAMPLES)' --gamma '$${mpcf_gamma}' --nlog2 '$${mpcf_nlog2}' --nbaf '$${mpcf_nbaf}'  --type allele_specific")

medicc/allele_specific/medicc/%/desc.txt : medicc/allele_specific/aspcf/%.RData
	$$(call RUN,-c -s 8G -m 12G -v $(ASCAT_ENV),"mkdir -p medicc/allele_specific/medicc && \
												 mkdir -p medicc/allele_specific/medicc/$$* && \
												 $(RSCRIPT) modules/test/phylogeny/initmedicc.R --sample_set $$* --type allele_specific")

medicc/allele_specific/medicc/%/tree_final.new : medicc/allele_specific/medicc/%/desc.txt
	$$(call RUN,-c -s 8G -m 12G -v $(MEDICC_ENV),"source $(MEDICC_VAR) && \
												  $(MEDICC_BIN)/medicc.py medicc/allele_specific/medicc/$$*/desc.txt medicc/allele_specific/medicc/$$* -v")
												  
medicc/allele_specific/medicc/%/tree_final.pdf : medicc/allele_specific/medicc/%/tree_final.new
	$$(call RUN,-c -n 12 -s 1G -m 2G -v $(PHYLO_ENV),"$(RSCRIPT) modules/test/phylogeny/plotmedicc.R --sample_set $$(*) --type allele_specific")

endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call allele-specific-medicc,$(set))))
		
else

PHONY += medicc medicc/total_copy medicc/total_copy/mad medicc/total_copy/mpcf medicc/total_copy/medicc

medicc : $(foreach set,$(SAMPLE_SETS),medicc/total_copy/medicc/$(set)/tree_final.new) $(foreach set,$(SAMPLE_SETS),medicc/total_copy/medicc/$(set)/tree_final.pdf)

define total-copy-medicc
medicc/total_copy/mad/%.RData : $(wildcard $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).Rdata))
	$$(call RUN,-c -s 8G -m 12G -v $(ASCAT_ENV),"mkdir -p medicc/total_copy && \
												 mkdir -p medicc/total_copy/mad && \
												 $(RSCRIPT) modules/test/phylogeny/combinesamples.R --sample_set $$* --normal_samples '$(NORMAL_SAMPLES)' --type total_copy")
												 
medicc/total_copy/mpcf/%.RData : medicc/total_copy/mad/%.RData
	$$(call RUN,-c -s 8G -m 12G -v $(ASCAT_ENV),"mkdir -p medicc/total_copy/mpcf && \
												 $(RSCRIPT) modules/test/phylogeny/segmentsamples.R --sample_set $$* --normal_samples '$(NORMAL_SAMPLES)' --gamma '$${mpcf_gamma}' --nlog2 '$${mpcf_nlog2}' --type total_copy")
												 
medicc/total_copy/medicc/%/desc.txt : medicc/total_copy/mpcf/%.RData
	$$(call RUN,-c -s 8G -m 12G -v $(ASCAT_ENV),"mkdir -p medicc/total_copy/medicc && \
												 mkdir -p medicc/total_copy/medicc/$$* && \
												 $(RSCRIPT) modules/test/phylogeny/initmedicc.R --sample_set $$* --type total_copy")
												 
medicc/total_copy/medicc/%/tree_final.new : medicc/total_copy/medicc/%/desc.txt
	$$(call RUN,-c -s 8G -m 12G -v $(MEDICC_ENV),"source $(MEDICC_VAR) && \
												  $(MEDICC_BIN)/medicc.py medicc/total_copy/medicc/$$*/desc.txt medicc/total_copy/medicc/$$* -t -v && \
												  cp medicc/total_copy/medicc/$$*/tree_fitch_nc.xml medicc/total_copy/medicc/$$*/tree_final.xml && \
												  cp medicc/total_copy/medicc/$$*/tree_fitch_nc.graph medicc/total_copy/medicc/$$*/tree_final.graph && \
												  cp medicc/total_copy/medicc/$$*/tree_fitch_nc.new medicc/total_copy/medicc/$$*/tree_final.new")
												  
medicc/total_copy/medicc/%/tree_final.pdf : medicc/total_copy/medicc/%/tree_final.new
	$$(call RUN,-c -n 12 -s 1G -m 2G -v $(PHYLO_ENV),"$(RSCRIPT) modules/test/phylogeny/plotmedicc.R --sample_set $$(*) --type total_copy")

endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call total-copy-medicc,$(set))))

endif

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
