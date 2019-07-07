include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/merge_alignments.$(NOW)
PHONY += fgbio

mearge_alignments : $(foreach sample,$(SAMPLES),fgbio/$(sample).sorted.bam)

TMPDIR ?= /home/${USER}/share/data/${USER}/tmp
JAVA = /home/${USER}/share/usr/jdk1.8.0_74/bin/java
PICARD = /home/${USER}/share/usr/picard/bin/picard.jar
REF_FASTA ?= /home/${USER}/share/reference/GATK_bundle/2.3/human_g1k_v37.fasta
POOL_A_INTERVAL ?= /home/${USER}/share/reference/target_panels/MSK-ACCESS-v1_0-probe-A.sorted.list
POOL_B_INTERVAL ?= /home/${USER}/share/reference/target_panels/MSK-ACCESS-v1_0-probe-B.sorted.list

define merge-alignments
fgbio/%.sorted.bam : fgbio/%.qn.sorted.bam
	$$(call RUN,-c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(JAVA) -Xmx8G -jar $(PICARD) SamToFastq \
									  I=$$^ \
									  FASTQ=/dev/stdout \
									  CLIPPING_ATTRIBUTE=XT \
									  CLIPPING_ACTION=2 \
									  INTERLEAVE=true \
									  NON_PF=true \
									  TMP_DIR=$(TMPDIR) | \
									  bwa mem -M -t 12 -R \"@RG\tID:$1\tLB:$1\tPL:illumina\tSM:$1\" \
									  -p $(REF_FASTA) /dev/stdin | \
									  $(JAVA) -Xmx8G -jar $(PICARD) SortSam \
									  I=/dev/stdin \
									  O=$1.sorted.bam \
									  SORT_ORDER=queryname \
									  TMP_DIR=$TMPDIR")
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call merge-alignments,$(sample))))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
