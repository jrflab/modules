include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/align_consensus.$(NOW)
PHONY += fgbio bam

align_consensus : $(foreach sample,$(SAMPLES),bam/$(sample).bam)

TMPDIR ?= /home/${USER}/share/data/${USER}/tmp
JAVA = /home/${USER}/share/usr/jdk1.8.0_74/bin/java
PICARD = /home/${USER}/share/usr/picard/bin/picard.jar
REF_FASTA ?= /home/${USER}/share/reference/GATK_bundle/2.3/human_g1k_v37.fasta
POOL_A_INTERVAL ?= /home/${USER}/share/reference/target_panels/MSK-ACCESS-v1_0-probe-A.sorted.list
POOL_B_INTERVAL ?= /home/${USER}/share/reference/target_panels/MSK-ACCESS-v1_0-probe-B.sorted.list

define bam-to-bam
fgbio/%.bam : fgbio/%.filtered.bam
	$$(call RUN,-c -n 12 -s 2G -m 4G,"set -o pipefail && \
									  $(JAVA) -Xmx8G -jar $(PICARD) SortSam \
									  I=fgbio/$$(*).filtered.bam \
									  O=/dev/stdout \
									  SORT_ORDER=queryname \
									  TMP_DIR=$(TMPDIR) | \
									  $(JAVA) -Xmx8G -jar $(PICARD) SamToFastq \
									  I=/dev/stdin \
									  FASTQ=/dev/stdout \
									  CLIPPING_ATTRIBUTE=XT \
									  CLIPPING_ACTION=2 \
									  INTERLEAVE=true \
									  NON_PF=true \
									  TMP_DIR=$(TMPDIR) | \
									  bwa mem -M -t 12 -R \"@RG\tID:$$(*)\tLB:$$(*)\tPL:illumina\tSM:$$(*)\" \
									  -p $REF /dev/stdin | \
									  $(JAVA) -Xmx8G -jar $(PICARD) SortSam \
									  I=/dev/stdin \
									  O=fgbio/$$(*).bam \
									  SORT_ORDER=coordinate \
									  TMP_DIR=$(TMPDIR)")

bam/%.bam : fgbio/%.bam
	$$(call RUN,-c -n 1 -s 4G -m 8G,"set -o pipefail && \
									 mkdir -p bam && \
									 cp fgbio/$$(*).bam bam/$$(*).bam && \
									 samtools index bam/$$(*).bam && \
									 cp bam/$$(*).bam.bai bam/$$(*).bai")
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call bam-to-bam,$(sample))))


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
