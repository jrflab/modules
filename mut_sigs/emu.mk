# Run emu on tumour/normal data
# Detect mutation signatures using mutect calls and control FreeC
##### DEFAULTS ######

include ~/share/modules/Makefile.inc

LOGDIR = log/emu.$(NOW)

EMU_PREPARE = $(HOME)/usr/bin/EMu-prepare
EMU = $(HOME)/usr/bin/EMu

PLOT_EMU = $(HOME)/share/scripts/plotEmuSignatures.R

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: all

ALL := emu/mutations.txt emu/cnv.txt emu/emu_results_bic.txt emu/report/index.html
ifdef NUM_SPECTRA
ALL += emu/emu_$(NUM_SPECTRA).timestamp
endif

all : $(ALL)

include ~/share/modules/variant_callers/somatic/mutect.mk

emu/mutations.txt : alltables/allTN.mutect.$(MUTECT_FILTER_SUFFIX).tab.txt
	$(INIT) awk 'NR > 1 { sub("X", "23", $$3); sub("Y", "24", $$3); sub("MT", "25", $$3); print $$1 "_" $$2, $$3, $$4, $$6 ">" $$7 }' $< | cat - $(EMU_REF_MUTATIONS) > $@

emu/cnv.txt : $(foreach pair,$(SAMPLE_PAIRS),freec/$(pair)/$(tumor.$(pair)).bam_CNVs)
	$(INIT) rm -f $@; for x in $^; do \
		sample=`echo $$x | sed 's:freec/::; s:/.*::'`; \
		awk -v sample=$$sample 'NR > 1 { sub("chr", "", $$1); sub("X", "23" , $$1); sub("Y", "24", $$1); sub("MT", "25", $$1); print sample, $$1, $$2, $$3, $$4; }' $$x >> $@; \
	done && cat $(EMU_REF_CNV) >> $@

emu/mutations.txt.mut.matrix : emu/mutations.txt emu/cnv.txt
	$(call LSCRIPT_MEM,4G,8G,"$(EMU_PREPARE) --chr $(EMU_REF_DIR) --cnv $(<<) --mut $< --pre $(@D) --regions $(EMU_TARGETS_FILE)")

emu/emu_results_bic.txt : emu/mutations.txt.mut.matrix
	$(call LSCRIPT_MEM,4G,8G,"$(EMU) --mut $< --opp human-exome --pre emu/emu_results")

ifdef NUM_SPECTRA
emu/emu_$(NUM_SPECTRA).timestamp : emu/mutations.txt.mut.matrix
	$(call LSCRIPT_MEM,4G,8G,"$(EMU) --force $(NUM_SPECTRA) --mut $< --opp human-exome --pre emu/emu_results && touch $@")
endif

emu/sample_pairs.txt : 
	$(INIT) echo "$(SAMPLE_PAIRS)" | sed 's/ /\n/' > $@

emu/report/index.html : emu/emu_results_bic.txt emu/sample_pairs.txt emu/mutations.txt emu/mutations.txt.mut.matrix
	$(call LSCRIPT_MEM,4G,8G,"$(PLOT_EMU) --inDir $(<D) --outDir $(@D) --sampleSubset $(<<) --mutations $(<<<) --samples $(<<<).samples")
