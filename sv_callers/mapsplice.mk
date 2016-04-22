# Run mapsplice
##### DEFAULTS ######

LOGDIR = log/mapsplice.$(NOW)

##### MAKE INCLUDES #####
include modules/Makefile.inc

MAPSPLICE = $(PYTHON) $(HOME)/share/usr/MapSplice-v2.1.7/mapsplice.py
MAPSPLICE_OPTS = -c $(MAPSPLICE_REF_DIR) -x $(MAPSPLICE_REF_BASENAME) --bam

ifeq ($(BAM_PHRED64),true)
MAPSPLICE_OPTS += --qual-scal phred64
else
MAPSPLICE_OPTS += --qual-scal phred33
endif

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: all

all : $(foreach sample,$(SAMPLES),mapsplice/$(sample).mapsplice_timestamp)

mapsplice/%.mapsplice_timestamp : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	$(call LSCRIPT_PARALLEL_MEM,6,2G,3G,"TMP1=`mktemp --tmpdir=$(TMPDIR)`.1.fastq; \
		TMP2=`mktemp --tmpdir=$(TMPDIR)`.2.fastq; \
		gzip -c $< > \$$TMP1; gzip -c $(<<) > \$$TMP2; \
		mkdir -p mapsplice/$*; \
		$(MAPSPLICE) $(MAPSPLICE_OPTS) -p 6 -o mapsplice/$* -1 \$$TMP1 -2 \$$TMP2 && touch $@; \
		rm \$$TMP1 \$$TMP2")
