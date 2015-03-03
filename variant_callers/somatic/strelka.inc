CONFIGURE_STRELKA = $(PERL) $(HOME)/share/usr/bin/configureStrelkaWorkflow.pl
STRELKA_CONFIG = $(HOME)/share/usr/etc/strelka_config.ini

VCF_GEN_IDS = DP FDP SDP SUBDP AU CU GU TU TAR TIR TOR
VCF_FIELDS += QSS TQSS NT QSS_NT TQSS_NT SGT SOMATIC

VARIANT_TYPES = strelka_snps strelka_indels
EFF_TYPES = silent missense nonsilent_cds nonsilent
STRELKA_FILTER_SUFFIX := pass.dbsnp.cosmic.nsfp
ifdef TARGETS_FILE
STRELKA_FILTER_SUFFIX := target_ft.$(STRELKA_FILTER_SUFFIX)
endif
STRELKA_FILTER_SUFFIX.strelka_snps := $(STRELKA_FILTER_SUFFIX).eff.chasm.fathmm
STRELKA_FILTER_SUFFIX.strelka_indels := $(STRELKA_FILTER_SUFFIX).eff

STRELKA_TABLE_SUFFIXES := $(foreach type,$(VARIANT_TYPES), $(type).$(STRELKA_FILTER_SUFFIX.$(type)).tab \
	$(foreach eff,$(EFF_TYPES),$(type).$(STRELKA_FILTER_SUFFIX.$(type)).tab.$(eff)))
STRELKA_TABLE_SUFFIXES := $(STRELKA_TABLE_SUFFIXES) $(addsuffix .novel,$(STRELKA_TABLE_SUFFIXES))
