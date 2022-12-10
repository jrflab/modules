include modules/Makefile.inc

LOGDIR ?= log/medicc2.$(NOW)

medicc : $(foreach sample,$(TUMOR_SAMPLES),medicc2/$(sample)/$(sample).txt) \
	 $(foreach set,$(SAMPLE_SETS),medicc2/$(set)/$(set).txt) \
	 $(foreach set,$(SAMPLE_SETS),medicc2/$(set)/$(set).tsv)

define collect-copy-number
medicc2/$1/$1.txt : facets/cncf/$1_$2.Rdata
	$$(call RUN,-c -n 1 -s 1G -m 2G -v $(MEDICC_ENV),"set -o pipefail && \
							  $(RSCRIPT) modules/copy_number/medicc2.R \
							  --option 1 \
							  --tumor_sample_name $1 \
							  --normal_sample_name $2 \
							  --file_in $$(<) \
							  --file_out $$(@)")
	
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call collect-copy-number,$(tumor.$(pair)),$(normal.$(pair)))))


define aggregate-copy-number
medicc2/$1/$1.txt : $(foreach sample,$(TUMOR_SAMPLES),medicc2/$(sample)/$(sample).txt)
	$$(call RUN,-c -n 1 -s 2G -m 4G -v $(MEDICC_ENV),"set -o pipefail && \
							  $(RSCRIPT) modules/copy_number/medicc2.R \
							  --option 2 \
							  --tumor_sample_name '$(tumors.$1)' \
							  --normal_sample_name '$(normal.$1)' \
							  --file_out $$(@)")

medicc2/$1/$1.tsv : medicc2/$1/$1.txt
	$$(call RUN,-c -n 1 -s 2G -m 4G -v $(MEDICC_ENV),"set -o pipefail && \
							  $(RSCRIPT) modules/copy_number/medicc2.R \
							  --option 3 \
							  --tumor_sample_name '$(tumors.$1)' \
							  --normal_sample_name '$(normal.$1)' \
							  --file_in $$(<) \
							  --file_out $$(@)")

endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call aggregate-copy-number,$(set))))


..DUMMY := $(shell mkdir -p version; \
	     $(MEDICC_ENV)/bin/R --version > version/medicc2.txt; \
	     $(MEDICC_ENV)/bin/medicc2 --help >> version/medicc2.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: medicc
