# Run wtsi NMF mutation sig on tumour/normal data
# Detect mutation signatures using mutect calls
##### DEFAULTS ######

include ~/share/modules/Makefile.inc

LOGDIR = log/nmf_mutsig.$(NOW)

EMU_PREPARE = $(HOME)/usr/bin/EMu-prepare
MATLAB = /usr/local/bin/matlab -nojvm -nodisplay -nosplash 

PLOT_EMU = $(RSCRIPT) $(HOME)/share/scripts/plotEmuSignatures.R


NMF_DIR = $(HOME)/usr/wtsi_mut_sig
PREPARE_NMF_INPUT = $(NMF_DIR)/createWTSIinput.m


.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: all

ALL := nmf_mutsig/mutations.txt.mut.matrix

all : $(ALL)

include ~/share/modules/variant_callers/somatic/mutect.inc

nmf_mutsig/mutations.txt : alltables/allTN.mutect.$(MUTECT_FILTER_SUFFIX).tab.txt
	$(INIT) awk 'NR > 1 { sub("X", "23", $$3); sub("Y", "24", $$3); sub("MT", "25", $$3); print $$1 "_" $$2, $$3, $$4, $$6 ">" $$7 }' $< > $@

nmf_mutsig/mutations.txt.mut.matrix : nmf_mutsig/mutations.txt
	$(INIT) $(EMU_PREPARE) --chr $(EMU_REF_DIR) --mut $< --pre $(@D) --regions $(EMU_TARGETS_FILE)

nmf_mutsig/input.mat : nmf_mutsig/mutations.txt.mut.matrix
	$(INIT) 
