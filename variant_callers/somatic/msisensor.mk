include modules/Makefile.inc

LOGDIR ?= log/msisensor.$(NOW)

msisensor: $(foreach pair,$(SAMPLE_PAIRS),msisensor/$(pair).msi) \
	   msisensor/msi.tsv

MSISENSOR_OPTS ?= -d $(REF_MSI) $(if $(TARGETS_FILE),-e $(TARGETS_FILE))

define msisensor-tumor-normal
msisensor/$1_$2.msi : bam/$1.bam bam/$2.bam bam/$1.bam.bai bam/$2.bam.bai
	$$(call RUN,-c -n 8 -s 1G -m 2G -v $(MSISENSOR_ENV),"set -o pipefail && \
							     msisensor msi $$(MSISENSOR_OPTS) \
							     -n $$(<<) \
							     -t $$(<) \
							     -b 8 \
							     -o $$(@)")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call msisensor-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

msisensor/msi.tsv : $(foreach pair,$(SAMPLE_PAIRS),msisensor/$(pair).msi)
	$(INIT) (head -1 $< | sed 's/^/sample\t/'; for x in $^; do sed "1d; s/^/$$(basename $$x)\t/" $$x; done | sed 's/_.*msi//' ) > $@

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: msisensor
