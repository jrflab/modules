# This module is used for running cufflinks
# input: $(SAMPLES) 
# Options: PHRED64 = true/false
# Authors: Fong Chun Chan <fongchunchan@gmail.com>
#
include modules/Makefile.inc

LOGDIR = log/cufflinks.$(NOW)


NUM_CORES ?= 8
CUFFLINKS = $(HOME)/share/usr/bin/cufflinks
CUFFCOMPARE = $(HOME)/share/usr/bin/cuffcompare
CUFFMERGE = $(HOME)/share/usr/bin/cuffmerge
CUFFDIFF = $(HOME)/share/usr/bin/cuffdiff
CUFFQUANT = $(HOME)/share/usr/bin/cuffquant
CUFFNORM = $(HOME)/share/usr/bin/cuffnorm
CUFFLINKS_OPTS = -b $(REF_FASTA) -u -g $(GENES_GTF) -p $(NUM_CORES) -u --no-update-check -v
CUFFCOMPARE_OPTS = -s $(REF_FASTA) -r $(GENES_GTF) -V

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
	$(call LSCRIPT_PARALLEL_MEM,$(NUM_CORES),2G,4G,"${CUFFLINKS} ${CUFFLINKS_OPTS} -o cufflinks/$* $<  && \
		mkdir -p cufflinks/gtf cufflinks/fpkm_tracking && \
		ln cufflinks/$*/transcripts.gtf cufflinks/gtf/$*.transcripts.gtf && \
		ln cufflinks/$*/isoforms.fpkm_tracking cufflinks/fpkm_tracking/$*.isoforms.fpkm_tracking && \
		ln cufflinks/$*/genes.fpkm_tracking cufflinks/fpkm_tracking/$*.genes.fpkm_tracking")

cufflinks/cuffcmp/cc.stats : $(foreach sample,$(SAMPLES),cufflinks/gtf/$(sample).transcripts.gtf)
	$(call LSCRIPT_MEM,10G,20G,"$(CUFFCOMPARE) $(CUFFCOMPARE_OPTS) -o $(@:.stats=) $^")

cufflinks/assembly_list.txt : $(foreach sample,$(SAMPLES),cufflinks/gtf/$(sample).transcripts.gtf)
	$(INIT) echo "$^" | tr ' ' '\n' > $@

cufflinks/gtf/merged.gtf : cufflinks/assembly_list.txt
	$(call LSCRIPT_PARALLEL_MEM,8,1G,2.5G,"$(CUFFMERGE) -o $(@D) -g $(GENES_GTF) -p 8 $<")

cufflinks/cxb/%.cxb : cufflinks/gtf/merged.gtf bam/%.bam
	$(call LSCRIPT_PARALLEL_MEM,4,1G,2.5G,"mkdir -p cufflinks/$* && \
	   	$(CUFFQUANT) -o cufflinks/$* -b $(REF_FASTA) -p 4 $^ && \
		ln cufflinks/$*/abundances.cxb $@")

cufflinks/cuffdiff/gene_exp.diff : cufflinks/gtf/merged.gtf $(foreach sample,$(SAMPLES),cufflinks/cxb/$(sample).cxb)
	$(call LSCRIPT_PARALLEL_MEM,8,1G,2.5G,"$(CUFFDIFF) -o $(@D) -p 8 $< \
		$(foreach pheno,$(PHENOTYPES),$(subst $( ),$(,),\
			$(foreach s,$(pheno.$(pheno)),cufflinks/cxb/$s.cxb))) \
		-L $(subst $( ),$(,),$(PHENOTYPES))")

cufflinks/cuffnorm/gene_exp.txt : cufflinks/gtf/merged.gtf $(foreach sample,$(SAMPLES),cufflinks/cxb/$(sample).cxb)
	$(call LSCRIPT_PARALLEL_MEM,8,1G,2.5G,"$(CUFFNORM) -o $(@D) -p 8 $< \
		$(foreach pheno,$(PHENOTYPES),$(subst $( ),$(,),\
			$(foreach s,$(pheno.$(pheno)),cufflinks/cxb/$s.cxb))) \
		-L $(subst $( ),$(,),$(PHENOTYPES))")
