include modules/Makefile.inc
include modules/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR = log/recurrent_mutations.$(NOW)

ALLTABLES_NONSYNONYMOUS_MUTECT = alltables/allTN.$(call SOMATIC_VCF_SUFFIXES,mutect).tab.nonsynonymous.txt
ALLTABLES_NONSYNONYMOUS_STRELKA_VARSCAN = alltables/allTN.strelka_varscan_indels.tab.nonsynonymous.txt

recurrent_mutations: recurrent_mutations/recurrent_mutations.tsv

recurrent_mutations/recurrent_mutations.tsv: $(ALLTABLES_NONSYNONYMOUS_MUTECT) $(ALLTABLES_NONSYNONYMOUS_STRELKA_VARSCAN)
	$(INIT) unset PYTHONPATH; \
	source $(ANACONDA_27_ENV)/bin/activate $(ANACONDA_27_ENV); \
	python modules/scripts/recurrent_mutations_plot.py $^ $(@D)
