
# vim: set ft=make :
# sub module containing vcf related tools

include modules/Makefile.inc
include modules/variant_callers/variantCaller.inc
include modules/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR = log/gemini.$(NOW)
GEMINI = unset PYTHONPATH; $(HOME)/share/usr/bin/gemini
GEMINI_LOAD_OPTS = -t snpEff

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: gemini

GEMINI_DB = gemini/gemini.db
gemini : $(if $(SAMPLE_PAIRS),gemini/mutect_gemini.timestamp gemini/strelka_gemini.timestamp gemini/varscan_gemini.timestamp)

gemini/samples.ped : $(SAMPLE_SET_FILE)
	$(INIT) perl -t$$' ' -lane 'BEGIN { print "#Family_ID\tIndividual_ID\tPaternal_ID\tMaternal_ID\tSex\tPhenotype\tEthnicity"; } $$i = 0; while ($$f = pop @F) { print "$$.\t$$f\t-9\t-9\t0\t" . (($$i++ == 0)? "1" : "2") . "\t-9"; }' $< > $@


GATK_VCFS = $(foreach sample,$(SAMPLES),\
			$(foreach type,gatk_snps gatk_indels,\
			vcf/$(sample).$(call FILTER_SUFFIX,$(type)).norm.vcf.gz))
gemini/gatk_gemini.timestamp : $(if $(SAMPLE_PAIRS),gemini/samples.ped) $(GATK_VCFS) $(addsuffix .tbi,$(GATK_VCFS))
	$(call LSCRIPT_PARALLEL_MEM,8,1G,2G,"for vcf in $(filter %.vcf.gz,$^); do $(GEMINI) load --cores 8 $(GEMINI_LOAD_OPTS) -v $$vcf $(if $(SAMPLE_PAIRS),-p $(filter %.ped,$^)) $(GEMINI_DB) ; done && touch $@")

ifdef SAMPLE_PAIRS
MUTECT_VCFS = $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).$(call SOMATIC_FILTER_SUFFIX,mutect).norm.vcf.gz)
VARSCAN_VCFS = $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).$(call SOMATIC_FILTER_SUFFIX,varscan_indels).norm.vcf.gz)
STRELKA_VCFS = $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).$(call SOMATIC_FILTER_SUFFIX,strelka_indels).rn.norm.vcf.gz)
gemini/mutect_gemini.timestamp : gemini/samples.ped $(MUTECT_VCFS) $(addsuffix .tbi,$(MUTECT_VCFS))
	$(call LSCRIPT_PARALLEL_MEM,8,1G,2G,"for vcf in $(filter %.vcf.gz,$^); do $(GEMINI) load --cores 8 $(GEMINI_LOAD_OPTS) -v \$$vcf -p $< $(GEMINI_DB); done && touch $@")
gemini/varscan_gemini.timestamp : gemini/samples.ped $(VARSCAN_VCFS) $(addsuffix .tbi,$(VARSCAN_VCFS))
	$(call LSCRIPT_PARALLEL_MEM,8,1G,2G,"for vcf in $(filter %.vcf.gz,$^); do $(GEMINI) load --cores 8 $(GEMINI_LOAD_OPTS) -v \$$vcf -p $< $(GEMINI_DB); done && touch $@")
gemini/strelka_gemini.timestamp : gemini/samples.ped $(STRELKA_VCFS) $(addsuffix .tbi,$(STRELKA_VCFS))
	$(call LSCRIPT_PARALLEL_MEM,8,1G,2G,"for vcf in $(filter %.vcf.gz,$^); do $(GEMINI) load --cores 8 $(GEMINI_LOAD_OPTS) -v \$$vcf -p $< $(GEMINI_DB); done && touch $@")
endif


include modules/vcf_tools/vcftools.mk
