LOGDIR = log/sufam_ss.$(NOW)

include modules/Makefile.inc

SUFAM_ENV = $(HOME)/share/usr/anaconda-envs/sufam-dev
SUFAM_OPTS = --format vcf --mpileup-parameters='-A -q 15 -Q 15 -d 15000'
SOMATIC_VCF2TSV = python modules/vcf_tools/somatic_vcf2tsv.py

ANNOTATE_SUFAM_GT_VCF = python modules/vcf_tools/annotate_sufam_gt_vcf.py

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY:

sufam_sampleset: $(foreach set,$(SAMPLE_SETS),vcf_ann/$(set).sufam.vcf) tsv/sufam_variants.tsv

define sufam-set
sufam/$1.set.vcf : $$(foreach pair,$$(pairs.$1),vcf_ann/$$(pair).somatic_variants.vcf.gz vcf_ann/$$(pair).somatic_variants.vcf.gz.tbi)
	$$(call RUN,-s 6G -m 8G,"bcftools merge -O v --force-samples $$(filter %.vcf.gz,$$^) > $$@")

vcf/$1.sufam.vcf : sufam/$1.set.vcf $$(foreach sample,$$(set.$1),bam/$$(sample).bam)
	$$(call RUN,-v $$(SUFAM_ENV) -s 2G -m 3G,"sufam --sample_name $$(set.$1) $$(SUFAM_OPTS) $$(REF_FASTA) $$^ > $$@")

vcf_ann/$1.sufam.vcf : vcf/$1.sufam.vcf $$(foreach pair,$$(pairs.$1),vcf_ann/$$(pair).somatic_variants.vcf.gz)
	$$(call RUN,-s 2G -m 3G,"$$(ANNOTATE_SUFAM_GT_VCF) $$^ > $$@")

tsv/$1.sufam.tsv : vcf_ann/$1.sufam.vcf.gz
	$$(call RUN,-s 4G -m 6G,"$$(SOMATIC_VCF2TSV) --normal $$(normal.$1) $$< > $$@")

endef
$(foreach set,$(SAMPLE_SETS),$(eval $(call sufam-set,$(set))))

tsv/sufam_variants.tsv : $(foreach set,$(SAMPLE_SETS),tsv/$(set).sufam.tsv)
	$(call RUN,-s 4G -m 6G,"(sed -n 1p $<; for x in $^; do sed 1d \$$x; done) > $@")


include modules/vcf_tools/vcftools.mk
