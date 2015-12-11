# Run wtsi NMF mutation sig on tumour/normal data
# Detect mutation signatures using mutect calls
##### DEFAULTS ######

include modules/Makefile.inc

LOGDIR = log/nmf_mutsig.$(NOW)

EMU_PREPARE = $(HOME)/usr/bin/EMu-prepare
MATLABPATH := modules/mut_sigs
ifeq ($(HOSTNAME),ika.cbio.mskcc.org)
export MATLAB_BIN := /usr/local/MATLAB/R2013a/bin/matlab
else
export MATLAB_BIN := /usr/local/bin/matlab
endif
MATLAB = export MATLABPATH=$(MATLABPATH); $(MATLAB_BIN) -nodisplay -nosplash 

NMF_DIR = $(HOME)/usr/nmf_mut_sig
NMF_TYPES_FILE = $(NMF_DIR)/types.mat

NMF_MIN_SIG = 1
NMF_MAX_SIG = 4

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: all

ALL := nmf_mutsig/mutations.txt.mut.matrix nmf_mutsig/results.mat nmf_mutsig/plot.timestamp

all : $(ALL)

include modules/variant_callers/somatic/somaticVariantCaller.inc

nmf_mutsig/mutations.txt : alltables/allTN.$(call SOMATIC_FILTER_SUFFIX,mutect).tab.txt
	$(INIT) awk 'NR > 1 { sub("X", "23", $$3); sub("Y", "24", $$3); sub("MT", "25", $$3); print $$1 "_" $$2, $$3, $$4, $$6 ">" $$7 }' $< > $@

nmf_mutsig/mutations.txt.mut.matrix : nmf_mutsig/mutations.txt
	$(INIT) $(EMU_PREPARE) --chr $(EMU_REF_DIR) --mut $< --pre $(@D) --regions $(EMU_TARGETS_FILE)

nmf_mutsig/input.mat : nmf_mutsig/mutations.txt.mut.matrix
	$(INIT) $(MATLAB) -r "createNMFinput $< $(<:.mut.matrix=.samples) $(NMF_TYPES_FILE) $(PROJECT_NAME) $@"

nmf_mutsig/results.mat : nmf_mutsig/input.mat
	$(INIT) $(MATLAB) -r "runNMF $< $(@:.mat=) $(NMF_DIR) $(NMF_MIN_SIG) $(NMF_MAX_SIG)"

nmf_mutsig/plot.timestamp : nmf_mutsig/results.mat
	$(INIT) $(MATLAB) -r "plotNMF $(<:.mat=) $(NMF_DIR) $(NMF_MIN_SIG) $(NMF_MAX_SIG)" && touch $@
