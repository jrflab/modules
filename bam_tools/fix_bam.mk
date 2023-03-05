include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR = log/fix_bam.$(NOW)

PICARD_JAR = ~/share/usr/picard/bin/picard.jar

fix_bam : $(foreach sample,$(SAMPLES),fixed_bam/$(sample).bam)

define fix-bam
unprocessed_bam/$1.ubam : unprocessed_bam/$1.bam
	$$(call RUN,-c -n 1 -s 12G -m 18G -w 72:00:00,"java -Djava.io.tmpdir=$(TMPDIR) -Xmx16G -jar $$(PICARD_JAR) RevertSam \
						       I=$$(<) \
						       O=$$(@) \
						       SANITIZE=true \
						       MAX_DISCARD_FRACTION=0.005 \
						       ATTRIBUTE_TO_CLEAR=XT \
						       ATTRIBUTE_TO_CLEAR=XN \
						       ATTRIBUTE_TO_CLEAR=AS \
						       ATTRIBUTE_TO_CLEAR=OC \
						       ATTRIBUTE_TO_CLEAR=OP \
						       SORT_ORDER=queryname \
						       RESTORE_ORIGINAL_QUALITIES=true \
						       REMOVE_DUPLICATE_INFORMATION=true \
						       REMOVE_ALIGNMENT_INFORMATION=true \
						       TMP_DIR=$(TMPDIR)")

unprocessed_bam/$1.fixed.bam : unprocessed_bam/$1.bam unprocessed_bam/$1.ubam
	$$(call RUN, -c -n 1 -s 12G -m 18G -w 72:00:00,"java -Djava.io.tmpdir=$(TMPDIR) -Xmx16G -jar $$(PICARD_JAR) MergeBamAlignment \
							R=$$(DMP_FASTA) \
							ALIGNED_BAM=$$(<).bam \
							UNMAPPED_BAM=$$(<<) \
							O=$$(@).fixed.bam \
							CREATE_INDEX=true \
							ADD_MATE_CIGAR=true \
							CLIP_ADAPTERS=true \
							CLIP_OVERLAPPING_READS=true \
							INCLUDE_SECONDARY_ALIGNMENTS=false \
							MAX_INSERTIONS_OR_DELETIONS=-1 \
							TMP_DIR=$(TMPDIR)")

unprocessed_bam/$1.dedup.bam : unprocessed_bam/$1.fixed.bam
	$$(call RUN, -c -n 1 -s 12G -m 18G -w 72:00:00,"java -Djava.io.tmpdir=$$(TMPDIR) -Xmx16G -jar $$(PICARD_JAR) MarkDuplicates \
							I=$$(<) \
							O=$$(@) \
							M=unprocessed_bam/$1.txt \
							TMP_DIR=$$(TMPDIR)")

fixed_bam/$1.bam : unprocessed_bam/$1.dedup.bam
	$$(call RUN, -c -n 1 -s 12G -m 18G -w 72:00:00,"java -Djava.io.tmpdir=$(TMPDIR) -Xmx16G -jar $$(PICARD_JAR) AddOrReplaceReadGroups \
							I=$$(<) \
							O=$$(@) \
							RGID=$1 \
							RGLB=$1 \
							RGPL=illumina \
							RGPU=NA \
							RGSM=$1 \
							TMP_DIR=$(TMPDIR) && \
							samtools index $$(@) && \
							cp fixed_bam/$1.bam.bai fixed_bam/$1.bai && \
							rm -rf unprocessed_bam/$1.ubam && \
							rm -rf unprocessed_bam/$1.fixed.bam && \
							rm -rf unprocessed_bam/$1.dedup.bam && \
							rm -rf unprocessed_bam/$1.fixed.bai && \
							rm -rf unprocessed_bam/$1.dedup.bai && \
							rm -rf unprocessed_bam/$1.txt")
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call fix-bam,$(sample))))

..DUMMY := $(shell mkdir -p version; \
             echo "picard" > version/fix_bam.txt)
.SECONDARY: 
.DELETE_ON_ERROR:
.PHONY: fix_bam
