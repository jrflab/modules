include modules/Makefile.inc

LOGDIR = log/split_rg.$(NOW)

splitrg :  $(foreach sample,$(SAMPLES), \
		  	$(foreach barcode,$(BARCODES),split_rg/$(sample)/$(barcode).bam))

#$(foreach sample,$(SAMPLES), \
#		  	$(foreach n,$(BARCODE_NUM),split_rg/$(sample)/$(BARCODES).$(n).bam.bai)) \
#	   $(foreach sample,$(SAMPLES), \
#		  	$(foreach n,$(BARCODE_NUM),split_rg/$(sample)/$(BARCODES).$(n).bai))

define split-rg
split_rg/$1/$2.bam : bam/$1.bam
	$$(call RUN,-n 1 -s 2G -m 4G,"set -o pipefail && \
				      $$(SAMTOOLS) view -b -r $2 $$(<) > $$(@)")

#split_rg/$1/$(BARCODES).$2.bam.bai : split_rg/$1/$(BARCODES).$2.bam
#	$$(call RUN,-n 1 -s 2G -m 4G,"set -o pipefail && \
#				      $$(SAMTOOLS) index $$(<)")

#split_rg/$1/$(BARCODES).$2.bai : split_rg/$1/$(BARCODES).$2.bam.bai
#	$$(call RUN,-n 1 -s 2G -m 4G,"set -o pipefail && \
#				      $$(CP) $$(<) $$(@)")

endef
$(foreach sample,$(SAMPLES), \
	$(foreach barcode,$(BARCODES), \
		$(eval $(call split-rg,$(sample),$(barcode)))))

..DUMMY := $(shell mkdir -p version; \
	     $(SAMTOOLS) --version > version/split_rg.txt;)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: splitrg
