include modules/Makefile.inc

SAM_TO_FASTQ = $(JAVA) -Xmx2G -jar $(JARDIR)/SamToFastq.jar VALIDATION_STRINGENCY=LENIENT
SAMPLE_FILE = samples.txt
SAMPLES = $(shell cat $(SAMPLE_FILE))

LOGDIR = log/pBwa.$(NOW)

BWA_BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)
SAMTOOLS_SORT_MEM = 2000000000
SEQ_PLATFORM = illumina

.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY : all bwa_bams

all : bwa_bams
bwa_bams : $(BWA_BAMS) $(addsuffix .bai,$(BWA_BAMS))

bwa/sai/%.1.sai bwa/sai/%.2.sai : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,8G,9G) -pe openmpi 10-100" $(MKDIR) $(@D) $(LOGDIR); $(PBWA) aln -f $(@D)/$* $(REF_FASTA) $^ 2> $(LOGDIR)/$(@F).log

%.bam.bai : %.bam
	SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,1G,2G)" $(MKDIR) $(@D) $(LOGDIR); $(SAMTOOLS) index $< &> $(LOGDIR)/$(@F).log

fastq/%.1.fastq.gz fastq/%.2.fastq.gz : gsc_bam/%.bam
	SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,10G,12G)" $(MKDIR) $(@D) $(LOGDIR); $(BAM2FASTQ) -o fastq/$*#.fastq $< &> $(LOGDIR)/$(@F).log && mv fastq/$*_1.fastq fastq/$*.1.fastq && mv fastq/$*_2.fastq fastq/$*.2.fastq && gzip fastq/$*.1.fastq fastq/$*.2.fastq

bwa/bam/%.bwa.sam : bwa/sai/%.1.sai bwa/sai/%.2.sai fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,7G,9G) -pe openmpi 10-100" $(MKDIR) $(@D) $(LOGDIR); \
			 LBID=`echo "$*" | sed 's/_[0-9]\+//'`; \
			 $(PBWA) sampe -f $@ -P -r "@RG\tID:$*\tLB:$${LBID}\tPL:${SEQ_PLATFORM}\tSM:$*" $(REF_FASTA) $(basename $(word 1,$^)) $(basename $(word 2,$^)) $(word 3,$^) $(word 4,$^) 2> $(LOGDIR)/$(@F).log


bam/%.bam : bwa/bam/%.bwa.sorted.filtered.fixmate.markdup.bam
	$(MKDIR) $(@D); ln -f $< $@

include modules/bam_tools/processBam.mk
