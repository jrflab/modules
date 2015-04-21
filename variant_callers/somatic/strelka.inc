ifndef STRELKA_INC
CONFIGURE_STRELKA = $(PERL) $(HOME)/share/usr/bin/configureStrelkaWorkflow.pl
STRELKA_CONFIG = $(HOME)/share/usr/etc/strelka_config.ini

VCF_GEN_IDS = DP FDP SDP SUBDP AU CU GU TU TAR TIR TOR
VCF_FIELDS += QSS TQSS NT QSS_NT TQSS_NT SGT SOMATIC
endif
STRELKA_INC = true
