ifndef SCALPEL_INC
SCALPEL_MIN_COV ?= 5

SCALPEL_DIR = $(HOME)/share/usr/scalpel-0.3.2
SCALPEL = export LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):$(SCALPEL_DIR)/bamtools-2.3.0/lib/; $(PERL) $(SCALPEL_DIR)/scalpel --mincov $(SCALPEL_MIN_COV)
SCALPEL_OPTS = --ref $(REF_FASTA) --validate --format annovar
ifeq ($(EXOME),true)
BED_DIR = $(HOME)/share/reference/splitExonBed/
BED_FILES = $(shell ls $(BED_DIR))
endif
ifdef TARGETS_FILE
SCALPEL_OPTS += --bed $(TARGETS_FILE)
endif


SCALPEL2VCF = $(PERL) scripts/scalpelToVcf.pl
endif
SCALPEL_INC = true
