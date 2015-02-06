# Run emu on tumour/normal data using absolute calls
# Detect mutation signatures using absolute results and control free-c results
##### DEFAULTS ######

include ~/share/modules/Makefile.inc

LOGDIR = log/emu_absolute.$(NOW)

EMU_PREPARE = $(HOME)/usr/bin/EMu-prepare
EMU = $(HOME)/usr/bin/EMu

PLOT_EMU = $(RSCRIPT) $(HOME)/share/scripts/plotEmuSignatures.R

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: all

ALL := emu_absolute/mutations.txt emu_absolute/emu_results_bic.txt emu_absolute/report/index.html emu_absolute/cnv.txt
ifdef NUM_SPECTRA
ALL += emu_absolute/emu_$(NUM_SPECTRA).timestamp
endif

all : $(ALL)

emu_absolute/mutations.txt : $(foreach pair,$(SAMPLE_PAIRS),absolute/tables/$(pair).absolute.txt)
	$(INIT) rm -f $@ && \
		for x in $^; do \
		sed 1d $$x | awk '$$139 == "snv" && $$141 == "TRUE" {print $$1_subclonal,$$6,$$7,$$140}' >> $@ && \
		sed 1d $$x | awk '$$139 == "snv" && $$142 == "TRUE" {print $$1_clonal,$$6,$$7,$$140}' >> $@; \
		done

emu_absolute/cnv.txt : $(foreach pair,$(SAMPLE_PAIRS),freec/$(pair)/$(tumor.$(pair)).bam_CNVs)
	$(INIT) rm -f $@; for x in $^; do \
		sample=`echo $$x | sed 's:freec/::; s:/.*::'`; \
		awk -v sample=$${sample}_clonal 'NR > 1 { sub("chr", "", $$1); sub("X", "23" , $$1); sub("Y", "24", $$1); sub("MT", "25", $$1); print sample, $$1, $$2, $$3, $$4; }' $$x >> $@; \
		awk -v sample=$${sample}_subclonal 'NR > 1 { sub("chr", "", $$1); sub("X", "23" , $$1); sub("Y", "24", $$1); sub("MT", "25", $$1); print sample, $$1, $$2, $$3, $$4; }' $$x >> $@; \
	done && cat $(EMU_REF_CNV) >> $@

emu_absolute/mutations.txt.mut.matrix : emu_absolute/mutations.txt emu_absolute/cnv.txt
	$(call LSCRIPT_MEM,4G,8G,"$(EMU_PREPARE) --chr $(EMU_REF_DIR) --cnv $(<<) --mut $< --pre $(@D) --regions $(EMU_TARGETS_FILE)")

emu_absolute/emu_results_bic.txt : emu_absolute/mutations.txt.mut.matrix
	$(call LSCRIPT_MEM,4G,8G,"$(EMU) --mut $< --opp human-exome --pre emu_absolute/emu_results")

RESULT_TIMESTAMPS = 
ifdef NUM_SPECTRA
emu_absolute/emu_$(NUM_SPECTRA).timestamp : emu_absolute/mutations.txt.mut.matrix
	$(call LSCRIPT_MEM,4G,8G,"$(EMU) --force $(NUM_SPECTRA) --mut $< --opp human-exome --pre emu_absolute/emu_results && touch $@")

RESULT_TIMESTAMPS += emu_absolute/emu_$(NUM_SPECTRA).timestamp
endif

emu_absolute/sample_pairs.txt : 
	$(INIT) echo "$(SAMPLE_PAIRS)" | sed 's/ /_subclonal\n/g' > $@ \
		&& echo "$(SAMPLE_PAIRS)" | sed 's/ /_clonal\n/g' >> $@

emu_absolute/report/index.html : emu_absolute/emu_results_bic.txt emu_absolute/sample_pairs.txt emu_absolute/mutations.txt $(RESULT_TIMESTAMPS)
	$(call LSCRIPT_MEM,4G,16G,"$(PLOT_EMU) --inPrefix $(<D)/emu_results --outDir $(@D) --sampleSubset $(<<) --mutations $(<<<) --samples $(<<<).samples")
