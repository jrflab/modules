include modules/Makefile.inc

LOGDIR = log/split_rg.$(NOW)

NUM_BARCODES = 5
BARCODE_NUM = $(shell seq 1 $(NUM_BARCODES))

splitrg :  $(foreach sample,$(SAMPLES), \
		  	$(foreach n,$(BARCODE_NUM),split_rg/$(sample)/$(BARCODE)$(n).bam)) \
	   $(foreach sample,$(SAMPLES), \
		  	$(foreach n,$(BARCODE_NUM),split_rg/$(sample)/$(BARCODE)$(n).bam.bai)) \
	   $(foreach sample,$(SAMPLES), \
		  	$(foreach n,$(BARCODE_NUM),split_rg/$(sample)/$(BARCODE)$(n).bai))

define split-rg
split_rg/$1/$2.bam : bam/$1.bam
	$$(call RUN,-n 1 -s 2G -m 4G,"set -o pipefail && \
				      $$(SAMTOOLS) view -b -r $1 $$(<) > $$(@)")

split_rg/$1/$2.bam.bai : split_rg/$1/$2.bam
	$$(call RUN,-n 1 -s 2G -m 4G,"set -o pipefail && \
				      $$(SAMTOOLS) index $$(<)")

split_rg/$1/$2.bai : split_rg/$1/$2.bam.bai
	$$(call RUN,-n 1 -s 2G -m 4G,"set -o pipefail && \
				      $$(CP) $$(<) $$(@)")

endef
$(foreach sample,$(SAMPLES), \
	$(foreach n,$(BARCODE_NUM), \
		$(eval $(call split-rg,$(sample),$(n)))))

..DUMMY := $(shell mkdir -p version; \
	     $(SAMTOOLS) --version > version/split_rg.txt;)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: splitrg
