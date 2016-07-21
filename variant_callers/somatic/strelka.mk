# Run strelka on tumour-normal matched pairs

include modules/Makefile.inc
include modules/variant_callers/gatk.inc
include modules/variant_callers/somatic/strelka.inc
include modules/variant_callers/somatic/somaticVariantCaller.inc
##### DEFAULTS ######


LOGDIR ?= log/strelka.$(NOW)
PHONY += strelka strelka_vcfs strelka_mafs

strelka : strelka_vcfs #strelka_mafs
	
VARIANT_TYPES := strelka_snps strelka_indels
strelka_vcfs : $(foreach type,$(VARIANT_TYPES),$(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).$(type).vcf))
strelka_mafs : $(foreach type,$(VARIANT_TYPES),$(foreach pair,$(SAMPLE_PAIRS),maf/$(pair).$(type).maf))

define strelka-tumor-normal
strelka/$1_$2/Makefile : bam/$1.bam bam/$2.bam
	$$(call LSCRIPT_NAMED,strelka_$1_$2,"rm -rf $$(@D) && $$(CONFIGURE_STRELKA) --tumor=$$< --normal=$$(<<) --ref=$$(REF_FASTA) --config=$$(STRELKA_CONFIG) --output-dir=$$(@D)")

strelka/$1_$2/task.complete : strelka/$1_$2/Makefile
	$$(call LSCRIPT_NAMED_PARALLEL_MEM,$1_$2.strelka,10,1G,1.5G,"make -j 10 -C $$(<D)")

vcf/$1_$2.%.vcf : strelka/vcf/$1_$2.%.vcf
	$$(INIT) perl -ne 'if (/^#CHROM/) { s/NORMAL/$2/; s/TUMOR/$1/; } print;' $$< > $$@ && $$(RM) $$<

strelka/vcf/$1_$2.strelka_snps.vcf : strelka/$1_$2/task.complete
	$$(INIT) cp -f strelka/$1_$2/results/all.somatic.snvs.vcf $$@

strelka/vcf/$1_$2.strelka_indels.vcf : strelka/$1_$2/task.complete
	$$(INIT) cp -f strelka/$1_$2/results/all.somatic.indels.vcf $$@

endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call strelka-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

include modules/vcf_tools/vcftools.mk

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

