# run pyrohmmvar: realignment-based variant calling method for 454 and ion torrent

include ~/share/modules/Makefile.inc
include ~/share/modules/variant_callers/gatk.inc

LOGDIR = log/pyrohmm.$(NOW)

PYROHMMVAR = $(HOME)/share/usr/bin/pyrohmmvar
PYROHMMVAR_MODEL = $(HOME)/share/reference/pyrohmm_parameter_config
PYROHMM2VCF = $(PERL) $(HOME)/share/scripts/pyroHMMVcf.pl

.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY: all

FILTER_SUFFIX := 

tables : $(foreach sample,$(SAMPLES),pyrohmm/tables/$(sample).pyrohmm.txt)
vcfs : $(foreach sample,$(SAMPLES),vcf/$(sample).pyrohmm.ann.vcf)

define pyrohmm-chr
pyrohmm/chr_tables/%.$1.pyrohmm.txt : bam/%.bam bam/%.bam.bai
	$$(call LSCRIPT_MEM,6G,10G,"$$(PYROHMMVAR) -r $1 -b $$< -f $$(REF_FASTA) -m $$(PYROHMMVAR_MODEL) > $$@")
endef
$(foreach chr,$(CHROMOSOMES), $(eval $(call pyrohmm-chr,$(chr))))

pyrohmm/tables/%.pyrohmm.txt : $(foreach chr,$(CHROMOSOMES),pyrohmm/chr_tables/%.$(chr).pyrohmm.txt)
	$(INIT) cat $^ > $@

vcf/%.pyrohmm.vcf : pyrohmm/tables/%.pyrohmm.txt
	$(INIT) $(PYROHMM2VCF) -f $(REF_FASTA) -n $* < $< | $(VCF_SORT) $(REF_DICT) - > $@ 2> $(LOG)

include ~/share/modules/vcf_tools/vcftools.mk
