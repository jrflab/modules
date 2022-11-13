include modules/Makefile.inc

LOGDIR ?= log/sv_signature.$(NOW)

MIN_SIZE = 1
MAX_SIZE = 10000000000000000
VIOLA_ENV = $(HOME)/share/usr/env/viola-sv-1.0.2

signature_sv :  $(foreach pair,$(SAMPLE_PAIRS),sv_signature/$(pair)/$(pair).merged.bed) \
		$(foreach pair,$(SAMPLE_PAIRS),sv_signature/$(pair)/$(pair).merged.bedpe) \
		$(foreach pair,$(SAMPLE_PAIRS),sv_signature/$(pair)/$(pair).merged.txt)
		
define signature-sv
sv_signature/$1_$2/$1_$2.merged.bed : vcf/$1_$2.merged_sv.vcf
	$$(call RUN,-c -n 1 -s 4G -m 8G -v $(SURVIVOR_ENV),"set -o pipefail && \
							    SURVIVOR vcftobed \
							    $$(<) \
							    $(MIN_SIZE) \
							    $(MAX_SIZE) \
							    $$(@)")
							    
sv_signature/$1_$2/$1_$2.merged.bedpe : sv_signature/$1_$2/$1_$2.merged.bed
	$$(call RUN,-c -n 1 -s 4G -m 8G,"set -o pipefail && \
					 echo \"chrom1	start1	end1	chrom2	start2	end2	sv_id	pe_support	strand1	strand2	svclass\" > \
					 $$(@) && \
					 cat $$(<) >> $$(@)")
					 
sv_signature/$1_$2/$1_$2.merged.txt : sv_signature/$1_$2/$1_$2.merged.bedpe
	$$(call RUN,-c -n 1 -s 4G -m 8G -v $(VIOLA_ENV),"set -o pipefail && \
							 python $(SCRIPTS_DIR)/sv_signature.py \
							 --bedpe_infile $$(<) \
							 --fragile_bed /data/reis-filho/lib/resource_files/viola/annotation/fragile_site.hg19.bed \
							 --timing_bedgraph /data/reis-filho/lib/resource_files/viola/annotation/replication_timing.bedgraph \
							 --sv_definitions /data/reis-filho/lib/resource_files/viola/definitions/sv_class_default.txt \
							 --text_outfile $$(@)")

endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call signature-sv,$(tumor.$(pair)),$(normal.$(pair)))))

..DUMMY := $(shell mkdir -p version; \
	     $(SURVIVOR_ENV)/bin/SURVIVOR --version &> version/sv_signature.txt; \
	     )
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: signature_sv
