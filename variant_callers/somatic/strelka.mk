# Run strelka on tumour-normal matched pairs

include modules/Makefile.inc
include modules/variant_callers/gatk.inc
##### DEFAULTS ######


LOGDIR ?= log/strelka.$(NOW)
PHONY += strelka strelka_vcfs strelka_mafs

CONFIGURE_STRELKA = $(PERL) $(HOME)/share/usr/bin/configureStrelkaWorkflow.pl
STRELKA_CONFIG = $(HOME)/share/usr/etc/strelka_config.ini
STRELKA_SOURCE_ANN_VCF = python modules/vcf_tools/annotate_source_vcf.py --source strelka

strelka : strelka_vcfs #strelka_mafs
	
STRELKA_VARIANT_TYPES := strelka_snps strelka_indels
strelka_vcfs : $(foreach type,$(STRELKA_VARIANT_TYPES),$(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).$(type).vcf))
strelka_mafs : $(foreach type,$(STRELKA_VARIANT_TYPES),$(foreach pair,$(SAMPLE_PAIRS),maf/$(pair).$(type).maf))

define strelka-tumor-normal
strelka/$1_$2/Makefile : bam/$1.bam bam/$2.bam
	$$(call RUN,-N strelka_$1_$2,"rm -rf $$(@D) && $$(CONFIGURE_STRELKA) --tumor=$$< --normal=$$(<<) --ref=$$(REF_FASTA) --config=$$(STRELKA_CONFIG) --output-dir=$$(@D)")

strelka/$1_$2/task.complete : strelka/$1_$2/Makefile
	$$(call RUN,-N $1_$2.strelka -n 10 -s 1G -m 1.5G,"make -j 10 -C $$(<D)")

strelka/vcf/$1_$2.%.vcf.tmp : strelka/vcf/$1_$2.%.vcf
	$$(call RUN,-s 1G -m 2G,"$$(RSCRIPT) modules/scripts/swapvcf.R --file $$< --tumor $1 --normal $2")

vcf/$1_$2.%.vcf : strelka/vcf/$1_$2.%.vcf.tmp
	$$(INIT) perl -ne 'if (/^#CHROM/) { s/NORMAL/$2/; s/TUMOR/$1/; } print;' $$< > $$@ && $$(RM) $$<
	
strelka/vcf/$1_$2.strelka_snps.vcf : strelka/$1_$2/task.complete
	$$(INIT) $$(STRELKA_SOURCE_ANN_VCF) < strelka/$1_$2/results/all.somatic.snvs.vcf > $$@

strelka/vcf/$1_$2.strelka_indels.vcf : strelka/$1_$2/task.complete
	$$(INIT) $$(STRELKA_SOURCE_ANN_VCF) < strelka/$1_$2/results/all.somatic.indels.vcf > $$@

endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call strelka-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

include modules/vcf_tools/vcftools.mk

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

