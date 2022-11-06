include modules/Makefile.inc

LOGDIR ?= log/get_basecount.$(NOW)

MAPQ := 0
BAQ := 0
COV := 0

getbasecount : $(foreach sample,$(SAMPLES),gbc/$(sample).txt.gz)

define get-basecount
gbc/$1.txt.gz : bam/$1.bam vcf/dataSilentNoPoleNotTertPromot.vcf
	$$(call RUN,-n 6 -s 3G -m 6G,"set -o pipefail && \
				      $(GBC) --fasta $(REF_FASTA) \
				      --bam $$(<) \
				      --vcf $$(<<) \
				      --output $$(@) \
				      --thread 6 \
				      --sort_output \
				      --compress_output \
				      --maq $(MAPQ) \
				      --baq $(BAQ) \
				      --cov $(COV) \
				      --filter_duplicate 0 \
				      --filter_improper_pair 0 \
				      --filter_qc_failed 1 \
				      --filter_indel 0 \
				      --filter_non_primary 1")
						    
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call get-basecount,$(sample))))

..DUMMY := $(shell mkdir -p version; \
	     ${GBC} &> version/get_basecount.txt;)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: getbasecount
