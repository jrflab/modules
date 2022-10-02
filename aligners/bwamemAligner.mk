include modules/Makefile.inc
include modules/variant_callers/gatk.inc
include modules/aligners/align.inc

ALIGNER := bwamem

LOGDIR ?= log/bwamem.$(NOW)

SAMTOOLS_SORT_MEM = 2000000000
SEQ_PLATFORM = illumina

VPATH ?= unprocessed_bam

# use fastq; otherwise use bams
FASTQ_CHUNKS := 10
FASTQ_CHUNK_SEQ := $(shell seq 1 $(FASTQ_CHUNKS))
FASTQUTILS = $(HOME)/share/usr/ngsutils/bin/fastqutils

BWA_ALN_OPTS ?= -M
BWAMEM_REF_FASTA ?= $(REF_FASTA)
BWAMEM_THREADS = 8
BWAMEM_MEM_PER_THREAD = $(if $(findstring true,$(PDX)),4G,2G)

BWA_BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)

bwamem : $(BWA_BAMS) $(addsuffix .bai,$(BWA_BAMS)) \
	 $(foreach sample,$(SAMPLES),metrics/$(sample).dedup_metrics.txt) \
	 metrics/dedup_metrics.txt

bam/%.bam : bwamem/bam/%.bwamem.$(BAM_SUFFIX)
	$(call RUN,,"ln -f $(<) $(@)")

define align-split-fastq
bwamem/bam/$2.bwamem.bam : $3
	$$(call RUN,-n $$(BWAMEM_THREADS) -s 1G -m $$(BWAMEM_MEM_PER_THREAD),"$$(BWA) mem -t $$(BWAMEM_THREADS) $$(BWA_ALN_OPTS) -R \"@RG\tID:$2\tLB:$1\tPL:$${SEQ_PLATFORM}\tSM:$1\" $$(BWAMEM_REF_FASTA) $$^ | $$(SAMTOOLS) view -bhS - > $$@")
endef
$(foreach ss,$(SPLIT_SAMPLES),$(if $(fq.$(ss)),$(eval $(call align-split-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))

bwamem/bam/%.bwamem.bam : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	LBID=`echo "$*" | sed 's/_[A-Za-z0-9]\+//'`; \
	$(call RUN,-n $(BWAMEM_THREADS) -s 1G -m $(BWAMEM_MEM_PER_THREAD),"$(BWA) mem -t $(BWAMEM_THREADS) $(BWA_ALN_OPTS) -R \"@RG\tID:$*\tLB:$${LBID}\tPL:${SEQ_PLATFORM}\tSM:$${LBID}\" $(BWAMEM_REF_FASTA) $^ | $(SAMTOOLS) view -bhS - > $@")

bwamem/bam/%.bwamem.bam : fastq/%.fastq.gz
	LBID=`echo "$*" | sed 's/_[A-Za-z0-9]\+//'`; \
	$(call RUN,-n $(BWAMEM_THREADS) -s 1G -m $(BWAMEM_MEM_PER_THREAD),"$(BWA) mem -t $(BWAMEM_THREADS) $(BWA_ALN_OPTS) -R \"@RG\tID:$*\tLB:$${LBID}\tPL:${SEQ_PLATFORM}\tSM:$${LBID}\" $(BWAMEM_REF_FASTA) $^ | $(SAMTOOLS) view -bhS - > $@")

fastq/%.fastq.gz : fastq/%.fastq
	$(call RUN,,"gzip -c $< > $(@) && $(RM) $<")
	
define dedup-metrics
metrics/$1.dedup_metrics.txt : bam/$1.bam
	$$(call RUN, -c -n 1 -s 16G -m 24G -v $(INNOVATION_ENV) -w 24:00:00, "set -o pipefail && \
									      picard \
									      -Xmx16G \
									      MarkDuplicates \
									      VALIDATION_STRINGENCY=LENIENT \
									      MAX_RECORDS_IN_RAM=4000000 \
									      TMP_DIR=$(TMPDIR) \
									      INPUT=$$(<) \
									      OUTPUT=/dev/null \
									      METRICS_FILE=$$(@)")

endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call dedup-metrics,$(sample))))
	
metrics/dedup_metrics.txt : $(foreach sample,$(SAMPLES),metrics/$(sample).dedup_metrics.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G -v $(INNOVATION_ENV),"set -o pipefail && \
							       $(RSCRIPT) $(SCRIPTS_DIR)/dedup_summary.R --option 1 --sample_names '$(SAMPLES)'")



..DUMMY := $(shell mkdir -p version; $(BWA) &> version/bwamem.txt; echo "options: $(BWA_ALN_OPTS)" >> version/bwamem.txt )
.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY: bwamem

include modules/bam_tools/processBam.mk
include modules/fastq_tools/fastq.mk
include modules/aligners/align.mk
