# Run strelka on tumour-normal matched pairs

include ~/share/modules/Makefile.inc
include ~/share/modules/variant_callers/gatk.inc
include ~/share/modules/variant_callers/strelka.inc
##### DEFAULTS ######


LOGDIR = log/strelka.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: all vcfs tables alltables


all : vcfs tables alltables
vcfs : $(foreach pair,$(SAMPLE_PAIRS),$(foreach type,$(VARIANT_TYPES),vcf/$(pair).$(type).$(STRELKA_FILTER_SUFFIX.$(type)).vcf))
tables : $(foreach suff,$(STRELKA_TABLE_SUFFIXES),$(foreach pair,$(SAMPLE_PAIRS),tables/$(pair).$(suff).txt) alltables/allTN.$(suff).txt) 

define strelka-tumor-normal
strelka/$1_$2/Makefile : bam/$1.bam bam/$2.bam
	$$(INIT) rm -rf $$(@D) && $$(CONFIGURE_STRELKA) --tumor=$$< --normal=$$(<<) --ref=$$(REF_FASTA) --config=$$(STRELKA_CONFIG) --output-dir=$$(@D) &> $$(LOG)

#$$(INIT) qmake -inherit -q jrf.q -- -j 20 -C $$< > $$(LOG) && touch $$@
strelka/$1_$2/task.complete : strelka/$1_$2/Makefile
	$$(call LSCRIPT_NAMED_PARALLEL_MEM,$1_$2.strelka,12,1G,1.5G,"make -j 12 -C $$(<D)")

vcf/$1_$2.strelka_snps.vcf : strelka/$1_$2/task.complete
	$$(INIT) cp -f strelka/$1_$2/results/all.somatic.snvs.vcf $$@

vcf/$1_$2.strelka_indels.vcf : strelka/$1_$2/task.complete
	$$(INIT) cp -f strelka/$1_$2/results/all.somatic.indels.vcf $$@
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call strelka-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

include ~/share/modules/vcf_tools/vcftools.mk
