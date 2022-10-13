include modules/Makefile.inc

LOGDIR = log/split_rg.$(NOW)

splitrg :  $(foreach sample,$(SAMPLES), \
		  	$(foreach barcode,$(BARCODES),split_rg/$(sample)/$(barcode).bam))

define split-rg
split_rg/$1/$2.bam : bam/$1.bam
	$$(call RUN,-n 1 -s 2G -m 4G,"set -o pipefail && \
				      $$(SAMTOOLS) view -b -r $2 $$(<) > $$(@)")

endef
$(foreach sample,$(SAMPLES), \
	$(foreach barcode,$(BARCODES), \
		$(eval $(call split-rg,$(sample),$(barcode)))))

..DUMMY := $(shell mkdir -p version; \
	     $(SAMTOOLS) --version > version/split_rg.txt;)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: splitrg
