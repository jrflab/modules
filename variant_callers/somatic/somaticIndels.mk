# Run multiple indel callers and then merge them
# taking calls found in >= 2 callers
include modules/Makefile.inc
LOGDIR ?= log/somatic_indels.$(NOW)

INDEL_TYPES ?= varscan_indels strelka_indels scalpel_indels lancet_indels platypus_indels pindel #mutect_indels

.PHONY : somatic_indels
somatic_indels : $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).somatic_indels.vcf)

UPS_INDEL_DIR = $(HOME)/share/usr/ups-indel
UPS_INDEL = if [ ! -f ext/SnpSift.jar ]; then mkdir -p ext; ln -fs $(UPS_INDEL_DIR)/ext/SnpSift.jar ext/SnpSift.jar; fi ; LD_LIBRARY_PATH=$(UPS_INDEL_DIR):$(HOME)/share/usr/lib $(UPS_INDEL_DIR)/ups_indel
UPS_SPLIT_CHR ?= true
MERGE_UVCF_VCF = python modules/vcf_tools/merge_uvcf_vcf.py
MERGE_INDEL_VCF = python modules/vcf_tools/merge_indel_vcf.py

vcf/%.somatic_indels.vcf : $(foreach type,$(INDEL_TYPES),vcf/%.$(type).uvcf.vcf)
	$(call RUN,-s 9G -m 12G -c,"$(MERGE_INDEL_VCF) $^ | $(VCF_SORT) $(REF_DICT) - > $@.tmp && $(call VERIFY_VCF,$@.tmp,$@)")

ifeq ($(UPS_SPLIT_CHR),true)
chr_vcf/%.chr_timestamp : vcf/%.vcf
	mkdir -p $(@D); for i in $(CHROMOSOMES); do (grep '^#' $<; grep -P "^$$i\t" $< || true) > $(@D)/$*.$$i.vcf ; done && touch $@

define uvcf-chr
chr_vcf/%.$1.vcf : chr_vcf/%.chr_timestamp
	if ! grep -q '^#CHROM' $$@; then rm -f $$< $$@; fi

chr_uvcf/%.$1.uvcf : chr_vcf/%.$1.vcf
	$$(call CHECK_UVCF,$$(call RUN,-s 4G -m 12G -c,"$$(UPS_INDEL) $$(REF_FASTA) $$< $$(@D)/$$*.$1 -hd=true"))
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call uvcf-chr,$(chr))))

uvcf/%.uvcf : $(foreach chr,$(CHROMOSOMES),chr_uvcf/%.$(chr).uvcf)
	$(INIT) h=`find $^ -size +1 | head -n1`; \
		if [[ ! -z $$h ]]; then \
		(grep -s '^#' $$h; grep -P '^CHROM\t' $$h; \
		for x in $^; do grep -v '^#' $$x || true; done) > $@; else touch $@; fi
else
uvcf/%.uvcf : vcf/%.vcf
	$(call RUN,-s 4G -m 12G -c,"$(UPS_INDEL) $(REF_FASTA) $< $(@D)/$*")
endif

vcf/%.uvcf.vcf : vcf/%.vcf uvcf/%.uvcf
	$(call CHECK_VCF,$(call RUN,-s 4G -m 6G,"$(MERGE_UVCF_VCF) $(<<) $(<) > $@"))

include modules/variant_callers/somatic/mutect2.mk
include modules/variant_callers/somatic/strelka.mk
include modules/variant_callers/somatic/scalpel.mk
include modules/variant_callers/somatic/lancet.mk
include modules/variant_callers/somatic/varscanTN.mk
include modules/variant_callers/somatic/platypus.mk
include modules/variant_callers/somatic/pindelTN.mk
