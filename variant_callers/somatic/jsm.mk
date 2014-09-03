# Run JointSNVMix and then Museq on Tumour-Normal paired samples
include ~/share/modules/Makefile.inc
include ~/share/modules/variant_callers/gatk.inc

#export PYTHONPATH = ~/share/usr/lib/python
#export PYTHONHOME = ~/share/usr/lib/python2.7
#export PYTHONHOME = ~/share/usr/lib/python
#JSM = PYTHONPATH=~/share/usr/lib/python2.7:~/share/usr/lib/python ~/share/usr/bin/python ~/share/usr/bin/jsm.py
JSM = $(HOME)/share/usr/bin/jsm.py

VPATH = bam jsm/vcf

LOGDIR = log/jsm.$(NOW)

SAMPLE_PAIR_FILE = sample_pairs.txt

TUMOR_SAMPLES ?= $(shell cut -f 1 $(SAMPLE_PAIR_FILE))
NORMAL_SAMPLES ?= $(shell cut -f 2 $(SAMPLE_PAIR_FILE))
SAMPLES ?= $(TUMOR_SAMPLES) $(NORMAL_SAMPLES)

$(foreach i,$(shell seq 1 $(words $(TUMOR_SAMPLES))),$(eval normal_lookup.$(word $i,$(TUMOR_SAMPLES)) := $(word $i,$(NORMAL_SAMPLES))))
$(foreach i,$(shell seq 1 $(words $(TUMOR_SAMPLES))),$(eval tumor_lookup.$(word $i,$(NORMAL_SAMPLES)) := $(word $i,$(TUMOR_SAMPLES))))

JSM_PRIORS_FILE = ~/share/usr/jsm/joint_priors.cfg
JSM_PARAMS_FILE = ~/share/usr/jsm/joint_params.cfg
JSM_THRESHOLD = 0.8

VCF_SORT = $(HOME)/share/usr/bin/vcfsorter.pl

VARIANT_ANNOTATOR = $(PERL) $(HOME)/share/scripts/variantAnnotator.pl

NSAMPLES = $(words $(TUMOR_SAMPLES))

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : all jsm_tables museq_vcf museq_tables

#all : museq_vcf museq_tables museq_table
all : jsm_tables


jsm_tables : $(foreach tumor,$(TUMOR_SAMPLES),jsm/tables/$(tumor)_$(normal_lookup.$(tumor)).jsm.filtered.txt)

museq_vcf : $(foreach tumor,$(TUMOR_SAMPLES),museq/vcf/$(tumor)_$(normal_lookup.$(tumor)).museq.filtered.annotated.vcf)

museq_tables : $(foreach tumor,$(TUMOR_SAMPLES),museq/tables/$(tumor)_$(normal_lookup.$(tumor)).museq.filtered.annotated.txt)
	
museq_table : museq/tables/all.museq.filtered.annotated.txt

#$(call jsm-train,normal.bam,tumor.bam)
define jsm-train
jsm/params/$1_$2.params : $1.bam $2.bam $1.bam.bai $2.bam.bai
	$$(call INIT_MEM,8G,10G) $$(JSM) train --priors_file $(JSM_PRIORS_FILE) --initial_parameters_file $(JSM_PARAMS_FILE) --skip_size 300 --model snvmix2 $$(REF_FASTA) $$(word 2,$$^) $$(word 1,$$^) $$@ &> $$(LOGDIR)/$$(@F).log
endef
$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call jsm-train,$(tumor),$(normal_lookup.$(tumor)))))

#$(call jsm-classify-chr,normal.bam,tumor.bam,chr)
define jsm-classify-chr
jsm/chr_tables/$1_$2.jsm.$3.txt : jsm/params/$1_$2.params $1.bam $2.bam
	$$(call INIT_MEM,8G,10G) $$(JSM) classify --chromosome $3 --model snvmix2 --parameters_file $$< --out_file $$@ --post_process $$(REF_FASTA) $$(word 3,$$^) $$(word 2,$$^) &> $$(LOGDIR)/$$(@F).log
endef
$(foreach chr,$(CHROMOSOMES),$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call jsm-classify-chr,$(tumor),$(normal_lookup.$(tumor)),$(chr)))))

jsm/tables/%.jsm.txt : $(foreach chr,$(CHROMOSOMES),jsm/chr_tables/%.jsm.$(chr).txt)
	$(INIT) head -1 $< > $@; \
	for i in $^; do \
		sed '1d' $$i >> $@; \
	done

#filter jsm results
#jsm/tables/%.jsm.filtered.txt : jsm/tables/%.jsm.txt
#	$(INIT) head -1 $< > $@ && \
#	awk '$$10 + $$11 > $(JSM_THRESHOLD) { print }' $< >> $@
jsm/tables/%.jsm.filtered.txt : jsm/tables/%.jsm.txt
	$(INIT) head -1 $< > $@ && \
	awk '$$18 > $(JSM_THRESHOLD) { print }' $< >> $@

include ~/share/modules/museqTN.mk
