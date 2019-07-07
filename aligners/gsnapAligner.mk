include modules/Makefile.inc
include modules/hg19.inc

SNP = snp135
NPARTS = 5 
SEQ = $(shell seq 0 $(shell expr $(NPARTS) \- 1))
OPTS = -d $(REF) -D ${GSNAP_REF_DIR} -B 4 -t 3 -A sam --novelsplicing=1 --pairexpect=150 -v $(SNP) -s $(REF).splicesites.iit -n 1 --quiet-if-excessive --nofails --gunzip
ifeq ($(BAM_PHRED64),true)
	OPTS += -J 64 -j -31
endif
# note: GSNAP actually uses 2 + number of specified threads

SAMPLE_FILE = samples.txt
SAMPLES = $(shell cat $(SAMPLE_FILE))

LOGDIR = gsnap/log

VPATH ?= unprocessed_bam

# Indicates how PCR-duplicates should be handled
# Options: rmdup, markdup, none
BAM_DUP_TYPE ?= rmdup
BAM_NO_RECAL ?= false
BAM_NO_REALN ?= false

GSNAP_BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)

gsnap_bams : $(GSNAP_BAMS) $(addsuffix .bai,$(GSNAP_BAMS))
gsnap_unmerged : $(foreach sample,$(SAMPLES),$(foreach i,$(SEQ),gsnap/sam/$(sample).gsnap.$i.sam))

.DELETE_ON_ERROR:

.SECONDARY:

#fastq/%.1.fastq fastq/%.2.fastq : gsc_bam/%.bam
#SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,10G,40G)" $(BAM2FASTQ) -o fastq/$*#.fastq $< &> $(LOGDIR)/$(@F).log && mv fastq/$*_1.fastq fastq/$*.1.fastq && mv fastq/$*_2.fastq fastq/$*.2.fastq

#$(call gsnap-part,part)
define gsnap-part
gsnap/bam/$(1).gsnap.$(2).bam : fastq/$1.1.fastq.gz fastq/$1.2.fastq.gz
	$$(call INIT_PARALLEL_MEM,4,3G,3.5G) $$(GSNAP) $$(OPTS) --read-group-id=$1 --part=$2/$$(NPARTS) $$^ 2> $(LOGDIR)/$1/$1.gsnap.$2.log | $$(SAMTOOLS) view -b -S - 2> $(LOGDIR)/$1/$1.gsnap.$2.samview.log | $$(SAMTOOLS) sort - $$(basename $$@) 2> $$(LOG)
endef
$(foreach sample,$(SAMPLES),$(eval $(foreach i,$(SEQ),$(eval $(call gsnap-part,$(sample),$i)))))
#$(foreach sample,$(SAMPLES),$(eval $(foreach i,$(SEQ),$(eval $(call gsnap-part,$(sample).readtrim,$i)))))

gsnap/bam/%.gsnap.bam : $(foreach i,$(SEQ),gsnap/bam/%.gsnap.$i.bam)
	SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,2G,3G)" $(SAMTOOLS) merge $@ $^ && rm -f $^

ifeq ($(BAM_DUP_TYPE),rmdup)
bam/%.bam : gsnap/bam/%.gsnap.sorted.filtered.rmdup.bam
	$(INIT) ln -f $< $@
else ifeq ($(BAM_DUP_TYPE),markdup) 
bam/%.bam : gsnap/bam/%.gsnap.sorted.filtered.markdup.bam
	$(INIT) ln -f $< $@
else ifeq ($(BAM_DUP_TYPE),none)
bam/%.bam : gsnap/bam/%.gsnap.sorted.filtered.bam
	$(INIT) ln -f $< $@
endif

bam/%.readtrim.bam : gsnap/bam/%.readtrim.gsnap.sorted.filtered.markdup.bam
	$(INIT) $(MKDIR) $(@D); ln -f $< $@

BAM_SUFFIX := gsnap.sorted.filtered

ifeq ($(BAM_DUP_TYPE),rmdup)
BAM_SUFFIX := $(BAM_SUFFIX).rmdup
else ifeq ($(BAM_DUP_TYPE),markdup) 
BAM_SUFFIX := $(BAM_SUFFIX).markdup
endif

ifeq ($(BAM_NO_RECAL),false)
BAM_SUFFIX := $(BAM_SUFFIX).recal
endif

ifeq ($(BAM_NO_REALN),false)
BAM_SUFFIX := $(BAM_SUFFIX).realn
endif

BAM_SUFFIX := $(BAM_SUFFIX).bam

bam/%.bam : gsnap/bam/%.$(BAM_SUFFIX)
	$(INIT) ln -f $< $@

#maps : $(GSNAP_MAPS)/$(REF).splicesites.iit $(GSNAP_MAPS)/$(SNP).iit

#$(GSNAP_MAPS)/%.splicesites.iit : refGene.txt.gz
#	zcat $< | psl_splicesites -s 1 | iit_store -o $@
#
#refGene.txt.gz :
#	wget ftp://hgdownload.cse.ucsc.edu/goldenPath/$(REF)/database/refGene.txt.gz
#
#$(SNP).txt.gz : 
#	wget ftp://hgdownload.cse.ucsc.edu/goldenPath/$(REF)/database/$(SNP).txt.gz
#
#$(GSNAP_MAPS)/$(SNP).iit : $(SNP).txt.gz
#	zcat $< | dbsnp_iit -w 1 | iit_store -o $@;\
#		snpindex -d $(REF) -v $(SNP) $@

#clean :
#	rm -f gsnap/sam/*.sam
#	rm -f gsnap/bam/*.bam
#	rm -fr gsnap/log/*

include modules/bam_tools/processBam.mk
include modules/fastq_tools/fastq.mk
