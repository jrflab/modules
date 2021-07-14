include modules/Makefile.inc

LOGDIR ?= log/get_basecount.$(NOW)
PHONY += getbasecount

GBC_ENV = $(HOME)/share/data/common/eec_sc_split/etc/conda
GBC_EXE = $(HOME)/share/data/common/eec_sc_split/etc/GetBaseCounts/GetBaseCounts
MAPQ := 10
BAQ := 15

getbasecount : $(foreach sample,$(SAMPLES),gbc/EEC128/$(sample).tsv)

define get-basecount
gbc/EEC128/$1.txt : bam/EEC128/$1.bam
	$$(call RUN,-n 6 -s 1G -m 2G -v $(GBC_ENV),"set -o pipefail && \
				      		    mkdir -p gbc/EEC128 && \
						    $(GBC_EXE) --fasta ~/share/reference/ucsc_gatk_bundle_2.8/ucsc.hg19.fasta \
						    --bam $$(<) \
						    --vcf etc/vcf/EEC128.vcf \
						    --output $$(@) \
						    --maq $(MAPQ) \
						    --baq $(BAQ) \
						    --filter_duplicate 0 \
						    --filter_improper_pair 0 \
						    --filter_qc_failed 1 \
						    --thread 6")
						    
gbc/EEC128/$1.tsv : gbc/EEC128/$1.txt
	$$(call RUN,-n 1 -s 2G -m 4G,"set -o pipefail && \
				      $(RSCRIPT) modules/variant_callers/getBaseCount.R --file_name $$(<)")


endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call get-basecount,$(sample))))

..DUMMY := $(shell mkdir -p version; \
	     /lila/home/brownd7/share/data/common/eec_sc_split/etc/GetBaseCounts/GetBaseCounts &> version/get_basecount.txt;)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: getbasecount
