ifndef MANTA_INC
CONFIG_MANTA = $(HOME)/share/usr/manta-0.29.6.centos5_x86_64/bin/configManta.py
CONFIG_MANTA_OPTS ?= --referenceFasta $(REF_FASTA) --config modules/sv_callers/manta_config.py.ini $(if $(TARGETS_FILE),--exome)
MANTA_HS_CONFIG = modules/sv_callers/manta_hs_config.py.ini
MANTA_CONFIG = modules/sv_callers/manta_config.py.ini
MANTA_HIGH_SENS ?= false
CONFIG_MANTA_OPTS = --referenceFasta $(REF_FASTA) \
					--config $(if $(findstring true,$(MANTA_HIGH_SENS)),\
                            $(MANTA_HS_CONFIG),$(MANTA_CONFIG)) \
                        $(if $(TARGETS_FILE),--exome) \
                        $(if $(MANTA_REGION),--region $(MANTA_REGION))
endif
MANTA_INC = true

