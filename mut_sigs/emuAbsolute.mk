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

SUBCLONAL := emu_absolute/subclonal_mutations.txt emu_absolute/subclonal_emu_results_bic.txt emu_absolute/report/subclonal.html
CLONAL += emu_absolute/clonal_mutations.txt emu_absolute/clonal_emu_results_bic.txt emu_absolute/report/clonal.html
ifdef NUM_SPECTRA
SUBCLONAL += emu_absolute/subclonal_emu_$(NUM_SPECTRA).timestamp
CLONAL += emu_absolute/clonal_emu_$(NUM_SPECTRA).timestamp
endif

all : $(SUBCLONAL) $(CLONAL) emu/cnv.txt 

emu_absolute/subclonal_mutations.txt : $(foreach pair,$(SAMPLE_PAIRS),absolute/tables/$(pair).absolute.txt)
	$(INIT) rm -f $@ && for x in $^; do sed 1d $$x | awk '$$139 == "snv" && $$141 == "TRUE" {print $$1,$$6,$$7,$$140}' >> $@; done

emu_absolute/clonal_mutations.txt : $(foreach pair,$(SAMPLE_PAIRS),absolute/tables/$(pair).absolute.txt)
	$(INIT) rm -f $@ && for x in $^; do sed 1d $$x | awk '$$139 == "snv" && $$142 == "TRUE" {print $$1,$$6,$$7,$$140}' >> $@; done

emu_absolute/cnv.txt : $(foreach pair,$(SAMPLE_PAIRS),freec/$(pair)/$(tumor.$(pair)).bam_CNVs)
	$(INIT) rm -f $@; for x in $^; do \
		sample=`echo $$x | sed 's:freec/::; s:/.*::'`; \
		awk -v sample=$$sample 'NR > 1 { sub("chr", "", $$1); sub("X", "23" , $$1); sub("Y", "24", $$1); sub("MT", "25", $$1); print sample, $$1, $$2, $$3, $$4; }' $$x >> $@; \
	done && cat $(EMU_REF_CNV) >> $@

emu_absolute/%_mutations.txt.mut.matrix : emu_absolute/%_mutations.txt emu_absolute/cnv.txt
	$(call LSCRIPT_MEM,4G,8G,"$(EMU_PREPARE) --chr $(EMU_REF_DIR) --cnv $(<<) --mut $< --pre $(@D) --regions $(EMU_TARGETS_FILE)")

emu_absolute/%_emu_results_bic.txt : emu_absolute/%_mutations.txt.mut.matrix
	$(call LSCRIPT_MEM,4G,8G,"$(EMU) --mut $< --opp human-exome --pre emu_absolute/$*_emu_results")

RESULT_TIMESTAMPS = 
ifdef NUM_SPECTRA
emu_absolute/%_emu_$(NUM_SPECTRA).timestamp : emu_absolute/%_mutations.txt.mut.matrix
	$(call LSCRIPT_MEM,4G,8G,"$(EMU) --force $(NUM_SPECTRA) --mut $< --opp human-exome --pre emu_absolute/$*_emu_results && touch $@")

RESULT_TIMESTAMPS += emu_absolute/%_emu_$(NUM_SPECTRA).timestamp
endif

emu_absolute/sample_pairs.txt : 
	$(INIT) echo "$(SAMPLE_PAIRS)" | sed 's/ /\n/g' > $@

emu_absolute/report/%.html : emu_absolute/%_emu_results_bic.txt emu_absolute/sample_pairs.txt emu_absolute/%_mutations.txt $(RESULT_TIMESTAMPS)
	$(call LSCRIPT_MEM,4G,16G,"$(PLOT_EMU) --inDir $(<D) --outDir $(@D) --sampleSubset $(<<) --mutations $(<<<) --samples $(<<<).samples")
