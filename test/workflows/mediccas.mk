include modules/Makefile.inc

LOGDIR ?= log/medicc.$(NOW)
PHONY += medicc medicc/allele_specific medicc/allele_specific/mad medicc/allele_specific/ascat medicc/allele_specific/aspcf

medicc : $(foreach set,$(SAMPLE_SETS),medicc/allele_specific/aspcf/$(set).RData)

define combine-samples
medicc/allele_specific/mad/%.RData : $(wildcard $(foreach pair,$(SAMPLE_PAIRS),facets/cncf/$(pair).Rdata))
	$$(call RUN,-c -s 8G -m 12G -v $(ASCAT_ENV),"mkdir -p medicc/allele_specific && \
												 mkdir -p medicc/allele_specific/mad && \
												 $(RSCRIPT) modules/test/phylogeny/combinesamples.R --sample_set $$* --normal_samples '$(NORMAL_SAMPLES)'")

endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call combine-samples,$(set))))

define ascat-aspcf
medicc/allele_specific/aspcf/%.RData : medicc/allele_specific/mad/%.RData
	$$(call RUN,-c -s 8G -m 12G -v $(ASCAT_ENV),"mkdir -p medicc/allele_specific/ascat && \
												 mkdir -p medicc/allele_specific/aspcf && \
												 $(RSCRIPT) modules/test/phylogeny/segmentsamples.R --sample_set $$* --normal_samples '$(NORMAL_SAMPLES)' --gamma '$${mpcf_gamma}' --nlog2 '$${mpcf_nlog2}' --nbaf '$${mpcf_nbaf}'")
endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call ascat-aspcf,$(set))))

#define init-medicc
#medicc/medicc/allele_specific/%/desc.txt medicc/medicc/total_copy/%/desc.txt : medicc/mpcf/%.RData
#	$$(call RUN,-c -s 8G -m 12G -v $(ASCAT_ENV),"if [ ! -d medicc/medicc ]; then mkdir medicc/medicc; fi && \
#												 if [ ! -d medicc/medicc/allele_specific ]; then mkdir medicc/medicc/allele_specific; fi && \
#												 if [ ! -d medicc/medicc/total_copy ]; then mkdir medicc/medicc/total_copy; fi && \
#												 if [ ! -d medicc/medicc/allele_specific/$$* ]; then mkdir medicc/medicc/allele_specific/$$*; fi && \
#												 if [ ! -d medicc/medicc/total_copy/$$* ]; then mkdir medicc/medicc/total_copy/$$*; fi && \
#												 $(RSCRIPT) modules/test/phylogeny/initmedicc.R --sample_set $$*")
#
#endef
#$(foreach set,$(SAMPLE_SETS),\
#		$(eval $(call init-medicc,$(set))))
#
#define run-medicc
#medicc/medicc/allele_specific/%/tree_final.new : medicc/medicc/allele_specific/%/desc.txt
#	$$(call RUN,-c -s 8G -m 12G -v $(MEDICC_ENV),"source $(MEDICC_VAR) && \
#												  $(MEDICC_BIN)/medicc.py medicc/medicc/allele_specific/$$*/desc.txt medicc/medicc/allele_specific/$$* -v")
#
#endef
#$(foreach set,$(SAMPLE_SETS),\
#		$(eval $(call run-medicc,$(set))))
#
#define boot-medicc
#medicc/boot/allele_specific/%/init.timestamp : medicc/aspcf/%.RData
#	$$(call RUN, -s 8G -m 12G -v $(ASCAT_ENV),"$(RSCRIPT) modules/test/phylogeny/bootstrapmedicc.R --sample_set $$*")
#
#medicc/boot/allele_specific/%/bootstrap.timestamp : medicc/boot/allele_specific/%/init.timestamp
#	$$(call RUN,-s 2G -m 4G -n 12 -v $(MEDICC_ENV) -w 3600,"source $(MEDICC_VAR) && \
#											  	   	   		seq -f '%03g' 1 100 | parallel -j 12 'if [ ! -f medicc/boot/allele_specific/$$*/{}/tree_final.new ]; then $(MEDICC_BIN)/medicc.py medicc/boot/allele_specific/$$*/{}/desc.txt medicc/boot/allele_specific/$$*/{}/ -v; fi' && \
#											  	   	   		seq -f '%03g' 1 100 | parallel -j 12 'if [ -f medicc/boot/allele_specific/$$*/{}/tree_final.new ]; then rm -rf medicc/boot/allele_specific/$$*/{}/desc.txt; fi' && \
#											  	   	   		seq -f '%03g' 1 100 | parallel -j 12 'if [ -f medicc/boot/allele_specific/$$*/{}/tree_final.new ]; then rm -rf medicc/boot/allele_specific/$$*/{}/*.fasta; fi' && \
#											  	   	   		seq -f '%03g' 1 100 | parallel -j 12 'if [ -f medicc/boot/allele_specific/$$*/{}/tree_final.new ]; then rm -rf medicc/boot/allele_specific/$$*/{}/chrom*; fi' && \
#											  	   	   		touch medicc/boot/allele_specific/$$*/bootstrap.timestamp")
#endef
#$(foreach set,$(SAMPLE_SETS),\
#		$(eval $(call boot-medicc,$(set))))
#		
#define plot-medicc
#medicc/plots/%_allele_specific.pdf medicc/plots/%_total_copy.pdf : medicc/medicc/allele_specific/%/tree_final.new medicc/medicc/total_copy/%/tree_final.new
#	$$(call RUN,-c -n 12 -s 1G -m 2G -v $(PHYLO_ENV),"$(RSCRIPT) modules/test/phylogeny/plotmedicc.R --sample_set $$(*)")
#
#endef
#$(foreach set,$(SAMPLE_SETS),\
#		$(eval $(call plot-medicc,$(set))))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
