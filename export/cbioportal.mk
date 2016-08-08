include modules/Makefile.inc
include modules/variant_callers/gatk.inc
include modules/variant_callers/variantCaller.inc
include modules/variant_callers/somatic/somaticVariantCaller.inc
SHELL = /bin/bash

LOGDIR = log/cbioportal.$(NOW)

.PHONY: mafs

mafs : $(foreach pair,$(SAMPLE_PAIRS),$(foreach caller,mutect strelka_varscan_indels,export/cbioportal/$(pair).$(caller).maf))

define CBIOPORTAL_VCF_RULE
unset PERL5LIB PERL_MB_OPT PERLBREW_ROOT PERL_LOCAL_LIB_ROOT PERL_MM_OPT && \
source activate /ifs/e63data/reis-filho/usr/anaconda-envs/vcf2maf && \
mkdir -p $(@D) && \
/opt/common/CentOS_6-dev/perl/perl-5.22.0/bin/perl \
	/ifs/e63data/reis-filho/usr/vcf2maf/vcf2maf.pl \
	--vep-path /opt/common/CentOS_6/vep/v84 \
	--vep-data /opt/common/CentOS_6/vep/v84 \
	--ref-fasta $(REF_FASTA) \
	--input-vcf $< \
	--output-maf $@ \
	--tumor-id $(tumor.$*) \
	--normal-id $(normal.$*)
endef

export/cbioportal/%.mutect.maf: vcf/%.$(call SOMATIC_VCF_SUFFIXES,mutect).vcf
	$(CBIOPORTAL_VCF_RULE)

export/cbioportal/%.strelka_varscan_indels.maf: vcf/%.strelka_varscan_indels.vcf
	$(CBIOPORTAL_VCF_RULE)

include modules/vcf_tools/vcftools.mk
