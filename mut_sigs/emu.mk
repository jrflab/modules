# Run emu on tumour/normal data
# Detect mutation signatures using mutect calls and control FreeC
##### DEFAULTS ######

include ~/share/modules/Makefile.inc

LOGDIR = log/emu.$(NOW)

EMU_PREPARE = $(HOME)/usr/bin/EMu-prepare
EMU = $(HOME)/usr/bin/EMu

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: all

all : emu/mutations.txt emu/cnv.txt emu/emu_results_bic.txt

include ~/share/modules/variant_callers/somatic/mutect.mk

emu/mutations.txt : alltables/allTN.mutect.$(FILTER_SUFFIX).tab.txt
	$(INIT) awk 'BEGIN {OFS = "\t"} NR > 1 { sub("X", "23", $$3); sub("Y", "24", $$3); sub("MT", "25", $$3); print $$1 "_" $$2, $$3, $$4, $$6 ">" $$7 }' $< > $@

emu/cnv.txt : $(foreach pair,$(SAMPLE_PAIRS),freec/$(pair)/$(tumor.$(pair)).bam_CNVs)
	$(INIT) rm -f $@; for x in $^; do \
		sample=`echo $$x | sed 's:freec/::; s:/.*::'`; \
		awk -v sample=$$sample 'BEGIN {OFS = "\t"} NR > 1 { sub("chr", "", $$1); sub("X", "23" , $$1); sub("Y", "24", $$1); sub("MT", "25", $$1); print sample, $$1, $$2, $$3, $$4; }' $$x >> $@; \
	done

emu/mutations.txt.mut.matrix : emu/mutations.txt emu/cnv.txt
	$(call LSCRIPT_MEM,4G,8G,"$(EMU_PREPARE) --chr $(EMU_REF_DIR) --cnv $(<<) --mut $< --pre $(@D) --regions $(TARGETS_FILE)")

emu/emu_results_bic.txt : emu/mutations.txt.mut.matrix
	$(call LSCRIPT_MEM,4G,8G,"$(EMU) --mut $< --opp human-exome --pre emu/emu_results")

