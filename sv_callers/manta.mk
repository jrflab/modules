# run manta on tumour-normal matched pairs

include modules/Makefile.inc
include modules/sv_callers/manta.inc

LOGDIR ?= log/manta.$(NOW)
PHONY += manta manta_vcfs

manta : manta_vcfs 

manta_vcfs: $(foreach sample,$(SAMPLES),vcf/$(sample).manta_sv.eff.vcf vcf/$(sample).manta_indels.eff.vcf vcf/$(sample).manta_candidate_sv.eff.vcf)

manta/%/runWorkflow.py : bam/%.bam bam/%.bam.bai
	$(INIT) $(CONFIG_MANTA) $(CONFIG_MANTA_OPTS) --tumorBam $< --runDir $(@D) 

manta/%/results/variants/tumorSV.vcf.gz manta/%/results/variants/candidateSmallIndels.vcf.gz manta/%/results/variants/candidateSV.vcf.gz : manta/%/runWorkflow.py
	$(call RUN,-n 8 -s 2G -m 2G,"python $< -m local -j 8")

vcf/%.manta_sv.vcf : manta/%/results/variants/tumorSV.vcf.gz
	$(INIT) zcat $< > $@

vcf/%.manta_indels.vcf : manta/%/results/variants/candidateSmallIndels.vcf.gz
	$(INIT) zcat $< > $@

vcf/%.manta_candidate_sv.vcf : manta/%/results/variants/candidateSV.vcf.gz
	$(INIT) zcat $< > $@

.PHONY: $(PHONY)

include modules/vcf_tools/vcftools.mk
