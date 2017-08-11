# Run mapsplice
##### DEFAULTS ######

LOGDIR = log/mapsplice.$(NOW)

##### MAKE INCLUDES #####
include modules/Makefile.inc

MAPSPLICE_ENV = $(HOME)/anaconda-envs/mapsplice-2.2.1
MAPSPLICE = mapsplice.py
MAPSPLICE_OPTS = -c $(MAPSPLICE_REF_DIR) -x $(MAPSPLICE_REF_BASENAME) --bam --gene-gtf $(GENES_GTF) --fusion

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
	$(call RUN,-n 6 -s 2G -m 3G,"TMP1=`mktemp --tmpdir=$(TMPDIR)`.1.fastq; \
		TMP2=`mktemp --tmpdir=$(TMPDIR)`.2.fastq; \
		gzip -c $< | sed 's:^@\(.*\)/1\$$:@\1:' > \$$TMP1; gzip -c $(<<) | 's:^@\(.*\)/2\$$:@\1:' > \$$TMP2; \
		mkdir -p mapsplice/$*; \
		$(MAPSPLICE) $(MAPSPLICE_OPTS) -p 6 -o mapsplice/$* -1 \$$TMP1 -2 \$$TMP2 && touch $@; \
		rm \$$TMP1 \$$TMP2")
