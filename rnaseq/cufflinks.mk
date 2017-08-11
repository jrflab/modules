# This module is used for running cufflinks
# input: $(SAMPLES) 
# Options: BAM_PHRED64 = true/false
# Authors: Fong Chun Chan <fongchunchan@gmail.com>
#
include modules/Makefile.inc

LOGDIR = log/cufflinks.$(NOW)


NUM_CORES ?= 8
CUFFLINKS = $(HOME)/share/usr/bin/cufflinks
CUFFLINKS_OPTS = -b $(REF_FASTA) -u -g $(GENES_GTF) -p $(NUM_CORES) -u --no-update-check -v
CUFFCOMPARE = $(HOME)/share/usr/bin/cuffcompare
CUFFCOMPARE_OPTS = --no-update-check
CUFFMERGE = $(HOME)/share/usr/bin/cuffmerge
CUFFMERGE_OPTS = --no-update-check
CUFFDIFF = $(HOME)/share/usr/bin/cuffdiff
CUFFDIFF_OPTS = --no-update-check -v
CUFFQUANT = $(HOME)/share/usr/bin/cuffquant
CUFFQUANT_OPTS = --no-update-check -v
CUFFNORM = $(HOME)/share/usr/bin/cuffnorm
CUFFNORM_OPTS = --no-update-check -v
CUFFCOMPARE_OPTS = --no-update-check -s $(REF_FASTA) -r $(GENES_GTF) -V -v

PHENO_FILE ?= pheno.txt
ifneq ($(wildcard $(PHENO_FILE)),)
  A = $(shell sed '1d' $(PHENO_FILE) | cut -f1)
  B = $(shell sed '1d' $(PHENO_FILE) | cut -f2)
  $(foreach i,$(shell seq 1 $(words $(A))),$(eval pheno.$(word $i,$(B)) += $(word $i,$(A))))
  PHENOTYPES = $(shell sed '1d' $(PHENO_FILE) | cut -f2 | sort | uniq)
endif

..DUMMY := $(shell mkdir -p version; $(CUFFLINKS) &> version/tophat.txt; echo "options: $(CUFFLINKS_OPTS)" >> version/cufflinks.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : all_cufflinks cufflinks cuffcmp cuffmerge cuffdiff cuffnorm

all_cufflinks : cufflinks cuffcmp cuffmerge cuffdiff cuffnorm
cufflinks : $(foreach sample,$(SAMPLES),cufflinks/gtf/$(sample).transcripts.gtf)
cuffcmp : cufflinks/cuffcmp/cc.stats
cuffmerge : cufflinks/gtf/merged.gtf
cuffdiff : cufflinks/cuffdiff/gene_exp.diff
cuffnorm : cufflinks/cuffnorm/gene_exp.txt

cufflinks/gtf/%.transcripts.gtf cufflinks/fpkm_tracking/%.isoforms.fpkm_tracking cufflinks/fpkm_tracking/%.genes.fpkm_tracking : bam/%.bam
	$(call RUN,-n $(NUM_CORES) -s 2G -m 4G,"${CUFFLINKS} ${CUFFLINKS_OPTS} -o cufflinks/$* $<  && \
		mkdir -p cufflinks/gtf cufflinks/fpkm_tracking && \
		ln cufflinks/$*/transcripts.gtf cufflinks/gtf/$*.transcripts.gtf && \
		ln cufflinks/$*/isoforms.fpkm_tracking cufflinks/fpkm_tracking/$*.isoforms.fpkm_tracking && \
		ln cufflinks/$*/genes.fpkm_tracking cufflinks/fpkm_tracking/$*.genes.fpkm_tracking")

cufflinks/cuffcmp/cc.stats : $(foreach sample,$(SAMPLES),cufflinks/gtf/$(sample).transcripts.gtf)
	$(call RUN,-s 10G -m 20G,"$(CUFFCOMPARE) $(CUFFCOMPARE_OPTS) -o $(@:.stats=) $^")

cufflinks/assembly_list.txt : $(foreach sample,$(SAMPLES),cufflinks/gtf/$(sample).transcripts.gtf)
	$(INIT) echo "$^" | tr ' ' '\n' > $@

cufflinks/gtf/merged.gtf : cufflinks/assembly_list.txt
	$(call RUN,-n 8 -s 1G -m 2.5G,"$(CUFFMERGE) $(CUFFMERGE_OPTS) -o $(@D) -g $(GENES_GTF) -p 8 $<")

cufflinks/cxb/%.cxb : cufflinks/gtf/merged.gtf bam/%.bam
	$(call RUN,-n 4 -s 1G -m 2.5G,"mkdir -p cufflinks/$* && \
	   	$(CUFFQUANT) $(CUFFQUANT_OPTS) -o cufflinks/$* -b $(REF_FASTA) -p 4 $^ && \
		ln cufflinks/$*/abundances.cxb $@")

cufflinks/cuffdiff/gene_exp.diff : cufflinks/gtf/merged.gtf $(foreach sample,$(SAMPLES),cufflinks/cxb/$(sample).cxb)
	$(call RUN,-n 8 -s 1G -m 4G,"$(CUFFDIFF) $(CUFFDIFF_OPTS) -o $(@D) -p 8 $< \
		$(foreach pheno,$(PHENOTYPES),$(subst $( ),$(,),$(foreach s,$(pheno.$(pheno)),cufflinks/cxb/$s.cxb))) \
		-L $(subst $( ),$(,),$(PHENOTYPES))")

cufflinks/cuffnorm/gene_exp.txt : cufflinks/gtf/merged.gtf $(foreach sample,$(SAMPLES),cufflinks/cxb/$(sample).cxb)
	$(call RUN,-n 8 -s 1G -m 2G,"$(CUFFNORM) $(CUFFNORM_OPTS) -o $(@D) -p 8 $< \
		$(foreach pheno,$(PHENOTYPES),$(subst $( ),$(,),$(foreach s,$(pheno.$(pheno)),cufflinks/cxb/$s.cxb))) \
		-L $(subst $( ),$(,),$(PHENOTYPES))")
