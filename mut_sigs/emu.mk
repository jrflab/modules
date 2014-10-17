# Run emu on tumour/normal data
# Detect mutation signatures using mutect calls
##### DEFAULTS ######

include ~/share/modules/Makefile.inc
include ~/share/modules/variant_callers/somatic/mutect.mk

LOGDIR = log/emu.$(NOW)

EMU_PREPARE = $(HOME)/usr/bin/EMu-prepare
EMU = $(HOME)/usr/bin/EMu
EMU_REF_DIR = $(MAPSPLICE_REF_DIR)

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: all

all : emu/mutations.txt emu/cnv.txt emu/emu_results_bic.txt

emu/mutations.txt : alltables/allTN.mutect.$(FILTER_SUFFIX).tab.txt
	$(INIT) awk 'BEGIN {OFS = "\t"} NR > 1 { sub("X", "23", $$3); sub("Y", "24", $$3); sub("MT", "25", $$3); print $$1 "_" $$2, $$3, $$4, $$6 ">" $$7 }' $< > $@

emu/cnv.txt : $(foreach pair,$(SAMPLE_PAIRS),freec/$(pair)/$(tumor.$(pair)).bam_CNVs)
	$(INIT) for x in $$^; do \
		sample=`echo $$x | sed 's:freec/::; s:/.*::'`; \
		awk -v sample=$$sample 'BEGIN {OFS = "\t"} NR > 1 { sub("chr", "", $$2); sub("X", "23" , $$2); sub("Y", "24", $$2); sub("MT", "25", $$2); print sample, $$2, $$3, $$4, $$11; }' $$x; \
	done

emu/mutations.txt.mut.matrix : emu/mutations.txt emu/cnv.txt
	$(EMU_PREPARE) --chr $(EMU_REF_DIR) --cnv $(<<) --mut $< --pre $(@D) -bin 10000 --regions $(TARGETS_FILE)

emu/emu_results_bic.txt : emu/mutations.txt.mut.matrix
	$(EMU) --mut $< --opp human-exome --pre emu/emu_results
