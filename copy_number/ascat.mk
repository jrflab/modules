include modules/Makefile.inc

LOGDIR ?= log/ascat.$(NOW)

ascat : $(foreach pair,$(SAMPLE_PAIRS),ascat/log2/$(pair).pdf) \
	$(foreach pair,$(SAMPLE_PAIRS),ascat/bafall/$(pair).pdf) \
	$(foreach pair,$(SAMPLE_PAIRS),ascat/bafhet/$(pair).pdf) \
	$(foreach pair,$(SAMPLE_PAIRS),ascat/mad/$(pair).RData) \
	$(foreach pair,$(SAMPLE_PAIRS),ascat/aspcf/$(pair).pdf) \
	$(foreach pair,$(SAMPLE_PAIRS),ascat/ascat/$(pair).RData) \
	$(foreach pair,$(SAMPLE_PAIRS),ascat/total/$(pair).pdf)

define ascat-plot-log2
ascat/log2/$1_$2.pdf : facets/cncf/$1_$2.RData
	$$(call RUN,-c -v $(ASCAT_ENV) -s 1G -m 2G,"$(RSCRIPT) $(RSCRIPT_ASCAT) --option 1 --file_in $$(<) --file_out $$(@)")
	
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call ascat-plot-log2,$(tumor.$(pair)),$(normal.$(pair)))))
		
define ascat-plot-bafall
ascat/bafall/$1_$2.pdf : facets/cncf/$1_$2.RData
	$$(call RUN,-c -v $(ASCAT_ENV) -s 1G -m 2G,"$(RSCRIPT) $(RSCRIPT_ASCAT) --option 2 --file_in $$(<) --file_out $$(@)")
	
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call ascat-plot-bafall,$(tumor.$(pair)),$(normal.$(pair)))))

define ascat-plot-bafhet
ascat/bafhet/$1_$2.pdf : facets/cncf/$1_$2.RData
	$$(call RUN,-c -v $(ASCAT_ENV) -s 1G -m 2G,"$(RSCRIPT) $(RSCRIPT_ASCAT) --option 3 --file_in $$(<) --file_out $$(@)")
	
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call ascat-plot-bafhet,$(tumor.$(pair)),$(normal.$(pair)))))
		
define ascat-aspcf
ascat/mad/$1_$2.RData : facets/cncf/$1_$2.RData
	$$(call RUN,-c -v $(ASCAT_ENV) -s 3G -m 6G,"$(RSCRIPT) $(RSCRIPT_ASCAT) --option 4 --file_in $$(<) --file_out $$(@) --gamma '$${aspcf_gamma.$1}'")
	
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call ascat-aspcf,$(tumor.$(pair)),$(normal.$(pair)))))

define ascat-plot-aspcf
ascat/aspcf/$1_$2.pdf : ascat/mad/$1_$2.RData
	$$(call RUN,-c -v $(ASCAT_ENV) -s 3G -m 6G,"$(RSCRIPT) $(RSCRIPT_ASCAT) --option 5 --file_in $$(<) --file_out $$(@) --nlog2 '$${aspcf_nlog2.$1}' --nbaf '$${aspcf_nbaf.$1}'")
	
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call ascat-plot-aspcf,$(tumor.$(pair)),$(normal.$(pair)))))
		
define ascat-run-ascat
ascat/ascat/$1_$2.RData : ascat/mad/$1_$2.RData
	$$(call RUN,-c -v $(ASCAT_ENV) -s 3G -m 6G,"$(RSCRIPT) $(RSCRIPT_ASCAT) --option 6 --file_in $$(<) --file_out $$(@) --rho '$${ascat_rho.$1}' --psi '$${ascat_psi.$1}' --nlog2 '$${aspcf_nlog2.$1}' --nbaf '$${aspcf_nbaf.$1}'")
	
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call ascat-run-ascat,$(tumor.$(pair)),$(normal.$(pair)))))

define ascat-run-total
ascat/total/$1_$2.pdf : facets/cncf/$1_$2.RData
	$$(call RUN,-c -v $(ASCAT_ENV) -s 6G -m 12G,"$(RSCRIPT) $(RSCRIPT_ASCAT) --option 7 --file_in $$(<) --file_out $$(@)")

endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call ascat-run-total,$(tumor.$(pair)),$(normal.$(pair)))))


define ascat-plot-chr
ascat/bychr/$1_$2/timestamp : facets/cncf/$1_$2.RData
	$$(call RUN, -v $(ASCAT_ENV) -s 6G -m 12G,"mkdir -p ascat/bychr/ && \
						   mkdir -p ascat/bychr/$1_$2 && \
						   $(RSCRIPT) $(RSCRIPT_ASCAT) --option 8 --file_in $$(<) --file_out $$(@)")
		
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call ascat-plot-chr,$(tumor.$(pair)),$(normal.$(pair)))))
		
.DUMMY := $(shell mkdir -p version; \
	    $(ASCAT_ENV)/bin/R --version &> version/ascat.txt)
.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY: ascat
