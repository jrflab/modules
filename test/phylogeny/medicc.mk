include modules/Makefile.inc

LOGDIR ?= log/medicc.$(NOW)
PHONY += medicc medicc/mad medicc/mpcf medicc/medicc medicc/boot medicc/plots

medicc : $(foreach set,$(SAMPLE_SETS),medicc/mad/$(set).RData) $(foreach set,$(SAMPLE_SETS),medicc/mpcf/$(set).RData) $(foreach set,$(SAMPLE_SETS),medicc/medicc/$(set)/desc.txt) $(foreach set,$(SAMPLE_SETS),medicc/medicc/$(set)/tree_final.new) $(foreach set,$(SAMPLE_SETS),medicc/boot/$(set)) $(foreach set,$(SAMPLE_SETS),medicc/boot/$(set)/init.timestamp) $(foreach set,$(SAMPLE_SETS),medicc/plots/$(set).pdf)

define combine-samples
medicc/mad/%.RData : $(wildcard $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).Rdata))
	$$(call RUN,-c -s 8G -m 12G -v $(ASCAT_ENV),"$(RSCRIPT) modules/test/phylogeny/combinesamples.R --sample_set $$* --normal_samples '$(NORMAL_SAMPLES)'")

endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call combine-samples,$(set))))

define ascat-mpcf
medicc/mpcf/%.RData : medicc/mad/%.RData
	$$(call RUN,-c -s 8G -m 12G -v $(ASCAT_ENV),"if [ ! -d medicc/ascat ]; then mkdir medicc/ascat; fi && \
												 $(RSCRIPT) modules/test/phylogeny/segmentsamples.R --sample_set $$* --normal_samples '$(NORMAL_SAMPLES)' --gamma '$${mpcf_gamma}' --nlog2 '$${mpcf_nlog2}' --nbaf '$${mpcf_nbaf}'")
endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call ascat-mpcf,$(set))))

define init-medicc
medicc/medicc/%/desc.txt : medicc/mpcf/%.RData
	$$(call RUN,-c -s 8G -m 12G -v $(ASCAT_ENV),"$(RSCRIPT) modules/test/phylogeny/initmedicc.R --sample_set $$*")

endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call init-medicc,$(set))))

define run-medicc
medicc/medicc/%/tree_final.new : medicc/medicc/%/desc.txt
	$$(call RUN,-c -s 8G -m 12G -v $(MEDICC_ENV),"source $(MEDICC_VAR) && \
												  $(MEDICC_BIN)/medicc.py medicc/medicc/$$*/desc.txt medicc/medicc/$$* -v")

endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call run-medicc,$(set))))

define boot-medicc
medicc/boot/% medicc/boot/%/init.timestamp : medicc/mpcf/%.RData
	$$(call RUN,-c -s 8G -m 12G -v $(ASCAT_ENV),"$(RSCRIPT) modules/test/phylogeny/bootstrapmedicc.R --sample_set $$* && \
												 touch medicc/boot/$$*/init.timestamp")

medicc/boot/%/bootstrap.timestamp : medicc/boot/%/init.timestamp
	$$(call RUN,-s 2G -m 4G -n 12 -v $(MEDICC_ENV) -w 3600,"source $(MEDICC_VAR) && \
											  	   	   		seq -f '%03g' 1 100 | parallel -j 12 'if [ ! -f medicc/boot/$$*/{}/tree_final.new ]; then $(MEDICC_BIN)/medicc.py medicc/boot/$$*/{}/desc.txt medicc/boot/$$*/{}/ -v; fi' && \
											  	   	   		seq -f '%03g' 1 100 | parallel -j 12 'if [ -f medicc/boot/$$*/{}/tree_final.new ]; then rm -rf medicc/boot/$$*/{}/desc.txt; fi' && \
											  	   	   		seq -f '%03g' 1 100 | parallel -j 12 'if [ -f medicc/boot/$$*/{}/tree_final.new ]; then rm -rf medicc/boot/$$*/{}/*.fasta; fi' && \
											  	   	   		seq -f '%03g' 1 100 | parallel -j 12 'if [ -f medicc/boot/$$*/{}/tree_final.new ]; then rm -rf medicc/boot/$$*/{}/chrom*; fi' && \
											  	   	   		touch medicc/boot/$$*/bootstrap.timestamp")
endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call boot-medicc,$(set))))
		
define plot-medicc
medicc/plots/%.pdf : medicc/medicc/%/tree_final.new medicc/boot/%/init.timestamp
	$$(call RUN,-c -n 12 -s 1G -m 2G -v $(PHYLO_ENV),"$(RSCRIPT) modules/test/phylogeny/plotmedicc.R --sample_set $$*)

endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call plot-medicc,$(set))))


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
