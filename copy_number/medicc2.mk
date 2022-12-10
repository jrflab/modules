include modules/Makefile.inc

LOGDIR ?= log/medicc2.$(NOW)

medicc : $(foreach sample,$(TUMOR_SAMPLES),medicc2/$(sample)/$(sample).txt)

define aggregate-copynumber
medicc2/$1/$1.txt : facets/cncf/$1_$2.Rdata
	$$(call RUN,-c -s 1G -m 2G -v $(MEDICC_ENV),"set -o pipefail && \
						    $(RSCRIPT) modules/copy_number/medicc.R \
						    --option 1 \
						    --tumor_sample $1 \
						    --normal_sample $2 \
						    --file_in $$(<) \
						    --file_out $$(@)")
	
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call aggregate-copynumber,$(tumor),$(normal))))


..DUMMY := $(shell mkdir -p version; \
	     $(MEDICC_ENV)/bin/R --version > version/medicc2.txt; \
	     $(MEDICC_ENV)/bin/medicc2 --help >> version/medicc2.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: medicc
