# This module runs ContEst on snp vcf files from gatk
# Author: inodb, limr

##### MAKE INCLUDES #####
include modules/Makefile.inc
include modules/variant_callers/gatk.inc

LOGDIR ?= log/contest.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: contest

contest : contest/all_contest.txt

# ContEst doing on-the-fly genotyping
define contest-tumor-normal
contest/$1_$2.contest.txt : bam/$1.bam bam/$2.bam
	$$(call LSCRIPT_MEM,12G,12G,"$$(call GATK_MEM2,4G) -T ContEst -I:eval $$(<) -I:genotype $$(<<) \
		-pf $$(HAPMAP_POP_FILE) -o $$(@) -R $$(REF_FASTA)")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call contest-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

contest/all_contest.txt : $(foreach pair,$(SAMPLE_PAIRS),contest/$(pair).contest.txt)
	( \
		head -1 $< | sed "s/^/sample\t/"; \
		for s in $(^); do \
			grep -P "META\t" $$s | sed "s/^/`basename $$s _contamination.txt`/"; \
		done | sort -rnk5,5; \
	) > $@
