# Run mapsplice
##### DEFAULTS ######

LOGDIR = log/mapsplice.$(NOW)

##### MAKE INCLUDES #####
include modules/Makefile.inc

MAPSPLICE_TO_USV = python modules/sv_callers/mapsplice2usv.py

MAPSPLICE_ENV = $(HOME)/share/usr/anaconda-envs/mapsplice-2.2.1
MAPSPLICE = mapsplice.py
MAPSPLICE_OPTS = -c $(MAPSPLICE_REF_DIR) -x $(MAPSPLICE_REF_BASENAME) --bam --gene-gtf $(GENES_GTF) --fusion

ifeq ($(BAM_PHRED64),true)
MAPSPLICE_OPTS += --qual-scal phred64
else
MAPSPLICE_OPTS += --qual-scal phred33
endif

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: mapsplice

mapsplice : $(foreach sample,$(SAMPLES),usv/$(sample.mapsplice.tsv)

mapsplice/%_mapsplice.timestamp : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	$(call RUN,-n 6 -s 2G -m 3G,"TMP1=`mktemp --tmpdir=$(TMPDIR)`.1.fastq; \
		TMP2=`mktemp --tmpdir=$(TMPDIR)`.2.fastq; \
		gzip -c $< > \$$TMP1; \
		gzip -c $(<<) > \$$TMP2; \
		mkdir -p mapsplice/$*; \
		$(MAPSPLICE) $(MAPSPLICE_OPTS) -p 6 -o mapsplice/$* -1 \$$TMP1 -2 \$$TMP2 && touch $@; \
		rm \$$TMP1 \$$TMP2")

usv/%.mapsplice.tsv : mapsplice/%_mapsplice.timestamp
	$(RUN,,"$(MAPSPLICE_TO_USV) < mapsplice/$*/fusions_not_well_annotated.txt mapsplice/%*/fusions_well_annotated.txt > $@")
