# Run strelka on tumour-normal matched pairs

include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc
##### DEFAULTS ######


LOGDIR = log/strelka.$(NOW)

CONFIGURE_STRELKA = $(HOME)/share/usr/bin/configureStrelkaWorkflow.pl
STRELKA_CONFIG = $(HOME)/share/usr/etc/strelka_config.ini

VCF_GEN_IDS = DP FDP SDP SUBDP AU CU GU TU

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: all vcfs tables alltables

VARIANT_TYPES = strelka_snps strelka_indels
EFF_TYPES = silent missense nonsilent_cds nonsilent
FILTER_SUFFIX := target_ft.pass.dbsnp
FILTER_SUFFIX.strelka_snps := $(FILTER_SUFFIX).nsfp.eff.chasm.fathmm.transfic
FILTER_SUFFIX.strelka_indels := $(FILTER_SUFFIX).eff

all : vcfs tables alltables
vcfs : $(foreach pair,$(SAMPLE_PAIRS),$(foreach type,$(VARIANT_TYPES),vcf/$(pair).$(type).$(FILTER_SUFFIX.$(type)).vcf))
tables : $(foreach pair,$(SAMPLE_PAIRS),$(foreach type,$(VARIANT_TYPES),$(foreach eff,$(EFF_TYPES),tables/$(pair).$(type).$(FILTER_SUFFIX.$(type)).tab.$(eff).txt)))
alltables : $(foreach type,$(VARIANT_TYPES),$(foreach eff,$(EFF_TYPES),alltables/allTN.$(type).$(FILTER_SUFFIX.$(type)).tab.$(eff).txt))

define strelka-tumor-normal
strelka/$1_$2 : bam/$1.bam bam/$2.bam
	$$(INIT) $$(CONFIGURE_STRELKA) --tumor=$$< --normal=$$(<<) --ref=$$(REF_FASTA) --config=$$(STRELKA_CONFIG) --output-dir=strelka/$1_$2 &> $$(LOG)

#$$(INIT) qmake -inherit -q jrf.q -- -j 20 -C $$< > $$(LOG) && touch $$@
strelka/$1_$2/task.complete : strelka/$1_$2
	$$(call LSCRIPT_NAMED_PARALLEL_MEM,$1_$2.strelka,12,1G,1.5G,"make -j 12 -C $$<")

vcf/$1_$2.strelka_snps.vcf : strelka/$1_$2/task.complete
	$$(INIT) cp -f strelka/$1_$2/results/all.somatic.snvs.vcf $$@

vcf/$1_$2.strelka_indels.vcf : strelka/$1_$2/task.complete
	$$(INIT) cp -f strelka/$1_$2/results/all.somatic.indels.vcf $$@
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call strelka-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

include ~/share/modules/vcftools.mk
