include modules/Makefile.inc

LOGDIR ?= log/ascat.$(NOW)

ascat : $(foreach pair,$(SAMPLE_PAIRS),ascat/log2/$(pair).pdf) \
	$(foreach pair,$(SAMPLE_PAIRS),ascat/bafall/$(pair).pdf) \
	$(foreach pair,$(SAMPLE_PAIRS),ascat/bafhet/$(pair).pdf)
#	$(foreach pair,$(SAMPLE_PAIRS),ascat/mad/$(pair).RData) \
#	$(foreach pair,$(SAMPLE_PAIRS),ascat/log2nbaf/$(pair).pdf) \
#	$(foreach pair,$(SAMPLE_PAIRS),ascat/ascat/$(pair).pdf) \
#	$(foreach pair,$(SAMPLE_PAIRS),ascat/total/$(pair).pdf) \
#	$(foreach pair,$(SAMPLE_PAIRS),ascat/bychr/$(pair)/timestamp)

define ascat-plot-log2
ascat/log2/$1_$2.pdf : facets/cncf/$1_$2.Rdata
	$$(call RUN,-c -v $(ASCAT_ENV) -s 1G -m 2G,"set -o pipefail && \
						    $(RSCRIPT) modules/copy_number/ascat.R \
						    --type log2 \
						    --file_in $$(<) \
						    --file_out ascat/log2/$1_$2.pdf")
	
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call ascat-plot-log2,$(tumor.$(pair)),$(normal.$(pair)))))
		
define ascat-plot-bafall
ascat/bafall/$1_$2.pdf : facets/cncf/$1_$2.Rdata
	$$(call RUN,-c -v $(ASCAT_ENV) -s 1G -m 2G,"set -o pipefail && \
						    $(RSCRIPT) modules/copy_number/ascat.R \
						    --type bafall \
						    --file_in $$(<) \
						    --file_out ascat/bafall/$1_$2.pdf")
	
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call ascat-plot-bafall,$(tumor.$(pair)),$(normal.$(pair)))))

define ascat-plot-bafhet
ascat/bafhet/$1_$2.pdf : facets/cncf/$1_$2.Rdata
	$$(call RUN,-c -v $(ASCAT_ENV) -s 1G -m 2G,"set -o pipefail && \
						    $(RSCRIPT) modules/copy_number/ascat.R \
						    --type bafhet \
						    --file_in $$(<) \
						    --file_out ascat/bafhet/$1_$2.pdf")
	
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call ascat-plot-bafhet,$(tumor.$(pair)),$(normal.$(pair)))))
		
define ascat-aspcf
ascat/mad/$1_$2.RData : facets/cncf/$1_$2.Rdata
	$$(call RUN,-c -v $(ASCAT_ENV) -s 3G -m 6G,"$(RSCRIPT) modules/copy_number/ascat.R --type aspcf --file_in $$< --file_out ascat/mad/$1_$2.RData --gamma '$${aspcf_gamma.$1}'")
	
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call ascat-aspcf,$(tumor.$(pair)),$(normal.$(pair)))))

define ascat-plot-aspcf
ascat/log2nbaf/$1_$2.pdf : ascat/mad/$1_$2.RData
	$$(call RUN,-c -v $(ASCAT_ENV) -s 3G -m 6G,"$(RSCRIPT) modules/copy_number/ascat.R --type plot-aspcf --file_in $$< --file_out ascat/log2nbaf/$1_$2.pdf --nlog2 '$${aspcf_nlog2.$1}' --nbaf '$${aspcf_nbaf.$1}'")
	
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call ascat-plot-aspcf,$(tumor.$(pair)),$(normal.$(pair)))))
		
define ascat-run-ascat
ascat/ascat/$1_$2.pdf : ascat/mad/$1_$2.RData
	$$(call RUN,-c -v $(ASCAT_ENV) -s 3G -m 6G,"$(RSCRIPT) modules/copy_number/ascat.R --type run-ascat --file_in $$< --file_out ascat/ascat/$1_$2.pdf --rho '$${ascat_rho.$1}' --psi '$${ascat_psi.$1}' --nlog2 '$${aspcf_nlog2.$1}' --nbaf '$${aspcf_nbaf.$1}'")
	
ascat/total/$1_$2.pdf : facets/cncf/$1_$2.Rdata ascat/ascat/$1_$2.pdf
	$$(call RUN,-c -v $(ASCAT_ENV) -s 6G -m 12G,"$(RSCRIPT) modules/copy_number/ascat.R --type total-copy --file_in $$< --file_out ascat/total/$1_$2.pdf")	

endef

$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call ascat-run-ascat,$(tumor.$(pair)),$(normal.$(pair)))))

define ascat-plot-chr
ascat/bychr/$1_$2/timestamp : facets/cncf/$1_$2.Rdata ascat/ascat/$1_$2.pdf
	$$(call RUN, -v $(ASCAT_ENV) -s 6G -m 12G,"mkdir -p ascat/bychr/ && \
											   mkdir -p ascat/bychr/$1_$2 && \
											   $(RSCRIPT) modules/copy_number/ascat.R --type plot-chr --file_in $$< --file_out ascat/bychr/$1_$2")
		
endef

$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call ascat-plot-chr,$(tumor.$(pair)),$(normal.$(pair)))))

..DUMMY := $(shell mkdir -p version; \
	     R --version > version/ascat.txt;)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: ascat
