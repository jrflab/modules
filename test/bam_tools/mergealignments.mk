include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/merge_alignments.$(NOW)
PHONY += fgbio

mearge_alignments : $(foreach sample,$(SAMPLES),fgbio/$(sample).regrouped.bam)

TMPDIR ?= /home/${USER}/share/data/${USER}/tmp
JAVA = /home/${USER}/share/usr/jdk1.8.0_74/bin/java
PICARD = /home/${USER}/share/usr/picard/bin/picard.jar
REF_FASTA ?= /home/${USER}/share/reference/GATK_bundle/2.3/human_g1k_v37.fasta
POOL_A_INTERVAL ?= /home/${USER}/share/reference/target_panels/MSK-ACCESS-v1_0-probe-A.sorted.list
POOL_B_INTERVAL ?= /home/${USER}/share/reference/target_panels/MSK-ACCESS-v1_0-probe-B.sorted.list

define merge-alignments
fgbio/%.qn.sorted.bam : fgbio/%.qn.sorted.ubam
	$$(call RUN,-c -n 12 -s 2G -m 4G,"set -o pipefail && \
									  $(JAVA) -Xmx8g -jar $(PICARD) SamToFastq \
									  I=$$^ \
									  FASTQ=/dev/stdout \
									  CLIPPING_ATTRIBUTE=XT \
									  CLIPPING_ACTION=N \
									  INTERLEAVE=true \
									  NON_PF=true \
									  TMP_DIR=$(TMPDIR) | \
									  bwa mem -M -t 12 -R \"@RG\tID:$1\tLB:$1\tPL:illumina\tSM:$1\" \
									  -p $(REF_FASTA) /dev/stdin | \
									  $(JAVA) -Xmx8g -jar $(PICARD) SortSam \
									  I=/dev/stdin \
									  O=$1.qn.sorted.bam \
									  SORT_ORDER=queryname \
									  TMP_DIR=$(TMPDIR)")
									  
fgbio/%.merged.bam : fgbio/%.qn.sorted.bam
	$$(call RUN,-c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(JAVA) -Xmx12g -jar $PICARD MergeBamAlignment \
									  VALIDATION_STRINGENCY=SILENT \
									  R=$(REF_FASTA) \
									  UNMAPPED_BAM=$1.qn.sorted.ubam \
									  ALIGNED_BAM=$1.qn.sorted.bam \
									  OUTPUT=$1.merged.bam \
									  CREATE_INDEX=true \
									  ADD_MATE_CIGAR=true \
									  CLIP_ADAPTERS=false \
									  CLIP_OVERLAPPING_READS=false \
									  INCLUDE_SECONDARY_ALIGNMENTS=true \
									  MAX_INSERTIONS_OR_DELETIONS=-1 \
									  PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
									  ATTRIBUTES_TO_RETAIN=XS \
									  SO=coordinate \
									  TMP_DIR=$(TMPDIR)")
									  
fgbio/%.regrouped.bam : fgbio/%.merged.bam
	$$(call RUN,-c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  samtools view -H $1.merged.bam > $1.sam && \
									  grep \"^@RG\" $1.sam | sed \"s/ID:A/ID:$1/g\" >> $1.sam && \
									  samtools reheader -P $1.sam $1.bam > $1.regrouped.bam && \
									  samtools index $1.regrouped.bam && \
									  mv $1.bam.bai $1.regrouped.bai")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call merge-alignments,$(sample))))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
