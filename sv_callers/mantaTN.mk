# run manta on tumour-normal matched pairs

include modules/Makefile.inc
include modules/sv_callers/manta.inc

LOGDIR ?= log/manta.$(NOW)
PHONY += manta manta_vcfs

manta : manta_vcfs 

manta_vcfs : $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).manta_sv.eff.vcf vcf/$(pair).manta_indels.eff.vcf vcf/$(pair).manta_candidate_sv.eff.vcf)

define manta-tumor-normal
manta/$1_$2/runWorkflow.py : bam/$1.bam bam/$2.bam bam/$1.bam.bai bam/$2.bam.bai
	$$(INIT) $$(CONFIG_MANTA) $$(CONFIG_MANTA_OPTS) --tumorBam $$< --normalBam $$(<<) --runDir $$(@D) 

manta/$1_$2.manta_timestamp : manta/$1_$2/runWorkflow.py
	$$(call RUN,-n 8 -s 2G -m 2G,"python $$< -m local -j 8 && touch $$@")

manta/$1_$2/results/variants/somaticSV.vcf.gz : manta/$1_$2.manta_timestamp

manta/$1_$2/results/variants/candidateSmallIndels.vcf.gz : manta/$1_$2.manta_timestamp

manta/$1_$2/results/variants/candidateSV.vcf.gz : manta/$1_$2.manta_timestamp

vcf/$1_$2.manta_indels.vcf : manta/$1_$2/results/variants/candidateSmallIndels.vcf.gz
	$$(INIT) zcat $$< > $$@

vcf/$1_$2.manta_sv.vcf : manta/$1_$2/results/variants/somaticSV.vcf.gz
	$$(INIT) zcat $$< > $$@

vcf/$1_$2.manta_candidate_sv.vcf : manta/$1_$2/results/variants/candidateSV.vcf.gz
	$$(INIT) zcat $$< > $$@
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call manta-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

.PHONY: $(PHONY)

include modules/vcf_tools/vcftools.mk
