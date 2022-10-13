include modules/Makefile.inc

LOGDIR ?= log/get_basecount.$(NOW)

MAPQ := 0
BAQ := 0
COV := 0

getbasecount : $(foreach sample,$(SAMPLES), \
		  	$(foreach barcode,$(BARCODES),get_bc/$(sample)/$(barcode).txt.gz)) \

define get-basecount
get_bc/$1/$2.txt.gz : split_rg/$1/$2.bam vcf/MSKCC_Weigelt_Mission_Bio_11132018.vcf
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
				      
#summary/$1_sum_alt.txt.gz : gbc/$1.txt.gz
#	$$(call RUN,-n 1 -s 4G -m 8G,"set -o pipefail && \
#				      $(RSCRIPT) $(SCRIPTS_DIR)/summary/get_basecount.R --option 1 --sample_name $1 && \
#				      gzip summary/$1_sum_alt.txt")
#
#summary/$1_ins_del.txt.gz : gbc/$1.txt.gz
#	$$(call RUN,-n 1 -s 4G -m 8G,"set -o pipefail && \
#				      $(RSCRIPT) $(SCRIPTS_DIR)/summary/get_basecount.R --option 2 --sample_name $1 && \
#				      gzip summary/$1_ins_del.txt")
#
#summary/$1_all_alt.txt.gz : gbc/$1.txt.gz
#	$$(call RUN,-n 1 -s 4G -m 8G,"set -o pipefail && \
#				      $(RSCRIPT) $(SCRIPTS_DIR)/summary/get_basecount.R --option 3 --sample_name $1 && \
#				      gzip summary/$1_all_alt.txt")

endef
$(foreach sample,$(SAMPLES), \
	$(foreach barcode,$(BARCODES), \
		$(eval $(call get-basecount,$(sample),$(barcode)))))

..DUMMY := $(shell mkdir -p version; \
	     ${GBC} &> version/get_basecount.txt;)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: getbasecount
