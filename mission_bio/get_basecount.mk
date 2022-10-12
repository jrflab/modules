include modules/Makefile.inc

LOGDIR ?= log/get_basecount.$(NOW)

MAPQ := 15
BAQ := 20

getbasecount : $(foreach sample,$(SAMPLES),gbc/$(sample).txt)

define get-basecount
gbc/$1.txt : bam/$1.bam vcf/$1.vcf
	$$(call RUN,-n 6 -s 3G -m 6G,"set -o pipefail && \
				      $(GBC) --fasta $(REF_FASTA) \
				      --bam $$(<) \
				      --vcf $$(<<) \
				      --output $$(@) \
				      --maq $(MAPQ) \
				      --baq $(BAQ) \
				      --filter_duplicate 1 \
				      --filter_improper_pair 1 \
				      --filter_qc_failed 1 \
				      --thread 6")
						    
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call get-basecount,$(sample))))

..DUMMY := $(shell mkdir -p version; \
	     ${GBC} &> version/get_basecount.txt;)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: getbasecount
