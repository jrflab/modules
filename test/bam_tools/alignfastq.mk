include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/align_fastq.$(NOW)
PHONY += marianas

BWA_ALN_OPTS ?= -M
BWAMEM_REF_FASTA ?= $(REF_FASTA)
BWAMEM_THREADS = 12
BWAMEM_MEM_PER_THREAD = 2G

align_fastq : $(foreach sample,$(SAMPLES),marianas/$(sample)/$(sample).bam)

define fastq-to-bam
marianas/$1/$1.bam : marianas/$1/$1_R1_umi-clipped.fastq.gz marianas/$1/$1_R2_umi-clipped.fastq.gz
	$$(call RUN,-c -n $(BWAMEM_THREADS) -s 1G -m $(BWAMEM_MEM_PER_THREAD),"set -o pipefail && \
																		   $(BWA) mem -t $(BWAMEM_THREADS) $(BWA_ALN_OPTS) -R \"@RG\tID:$1\tLB:$1\tPL:${SEQ_PLATFORM}\tSM:$1\" $(BWAMEM_REF_FASTA) $$(^) | $(SAMTOOLS) view -bhS - > $$(@)")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call fastq-to-bam,$(sample))))


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
