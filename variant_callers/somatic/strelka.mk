# Run strelka on tumour-normal matched pairs

include ~/share/modules/Makefile.inc
include ~/share/modules/variant_callers/gatk.inc
##### DEFAULTS ######


LOGDIR = log/strelka.$(NOW)

CONFIGURE_STRELKA = $(PERL) $(HOME)/share/usr/bin/configureStrelkaWorkflow.pl
STRELKA_CONFIG = $(HOME)/share/usr/etc/strelka_config.ini

VCF_GEN_IDS = DP FDP SDP SUBDP AU CU GU TU TAR TIR TOR
VCF_FIELDS += QSS TQSS NT QSS_NT TQSS_NT SGT SOMATIC

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: all vcfs tables alltables

VARIANT_TYPES = strelka_snps strelka_indels
EFF_TYPES = silent missense nonsilent_cds nonsilent
FILTER_SUFFIX := pass.dbsnp
ifdef TARGETS_FILE
FILTER_SUFFIX := target_ft.$(FILTER_SUFFIX)
endif
FILTER_SUFFIX.strelka_snps := $(FILTER_SUFFIX).eff.nsfp.chasm.transfic#.fathmm
FILTER_SUFFIX.strelka_indels := $(FILTER_SUFFIX).eff

TABLE_SUFFIXES := $(foreach type,$(VARIANT_TYPES), $(type).$(FILTER_SUFFIX.$(type)).tab \
	$(foreach eff,$(EFF_TYPES),$(type).$(FILTER_SUFFIX.$(type)).tab.$(eff)))
TABLE_SUFFIXES := $(TABLE_SUFFIXES) $(addsuffix .novel,$(TABLE_SUFFIXES))

all : vcfs tables alltables
vcfs : $(foreach pair,$(SAMPLE_PAIRS),$(foreach type,$(VARIANT_TYPES),vcf/$(pair).$(type).$(FILTER_SUFFIX.$(type)).vcf))
tables : $(foreach suff,$(TABLE_SUFFIXES),$(foreach pair,$(SAMPLE_PAIRS),tables/$(pair).$(suff).txt) alltables/allTN.$(suff).txt) 

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
