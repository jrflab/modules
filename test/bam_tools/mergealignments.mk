include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/merge_alignments.$(NOW)
PHONY += fgbio

mearge_alignments : $(foreach sample,$(SAMPLES),fgbio/$(sample).regrouped.bam)

JAVA = /home/${USER}/share/usr/jdk1.8.0_74/bin/java
PICARD = /home/${USER}/share/usr/picard/bin/picard.jar
POOL_A_INTERVAL ?= /home/${USER}/share/reference/target_panels/MSK-ACCESS-v1_0-probe-A.sorted.list
POOL_B_INTERVAL ?= /home/${USER}/share/reference/target_panels/MSK-ACCESS-v1_0-probe-B.sorted.list

define merge-alignments
fgbio/%.qn.sorted.bam : fgbio/%.qn.sorted.ubam
	$$(call RUN,-c -n 12 -s 2G -m 4G,"set -o pipefail && \
									  $(JAVA) -Xmx16G -jar $(PICARD) SamToFastq \
									  I=fgbio/$$(*).qn.sorted.ubam \
									  FASTQ=/dev/stdout \
									  CLIPPING_ATTRIBUTE=XT \
									  CLIPPING_ACTION=N \
									  INTERLEAVE=true \
									  NON_PF=true \
									  TMP_DIR=$(TMPDIR) | \
									  bwa mem -M -t 12 -R \"@RG\tID:$$(*)\tLB:$$(*)\tPL:illumina\tSM:$$(*)\" \
									  -p $(REF_FASTA) /dev/stdin | \
									  $(JAVA) -Xmx16G -jar $(PICARD) FixMateInformation \
									  I=/dev/stdin \
									  O=fgbio/$$(*).mate.fixed.bam \
									  SORT_ORDER=coordinate \
									  TMP_DIR=$(TMPDIR) \
									  VALIDATION_STRINGENCY=LENIENT && \
									  $(JAVA) -Xmx8G -jar $(PICARD) SortSam \
									  I=fgbio/$$(*).mate.fixed.bam \
									  O=fgbio/$$(*).qn.sorted.bam \
									  SORT_ORDER=queryname \
									  TMP_DIR=$(TMPDIR)")
									  
fgbio/%.merged.bam : fgbio/%.qn.sorted.bam
	$$(call RUN,-c -n 1 -s 12G -m 24G -w 2880,"set -o pipefail && \
									  		   $(JAVA) -Xmx12G -jar $(PICARD) MergeBamAlignment \
									  		   VALIDATION_STRINGENCY=SILENT \
									  		   R=$(REF_FASTA) \
									  		   UNMAPPED_BAM=fgbio/$$(*).qn.sorted.ubam \
									  		   ALIGNED_BAM=fgbio/$$(*).qn.sorted.bam \
									  		   OUTPUT=fgbio/$$(*).merged.bam \
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
									  samtools view -H fgbio/$$(*).merged.bam > fgbio/$$(*).sam && \
									  grep \"^@RG\" fgbio/$$(*).sam | sed \"s/ID:A/ID:$$(*)/g\" >> fgbio/$$(*).sam && \
									  samtools reheader -P fgbio/$$(*).sam fgbio/$$(*).merged.bam > fgbio/$$(*).regrouped.bam && \
									  samtools index fgbio/$$(*).regrouped.bam && \
									  mv fgbio/$$(*).regrouped.bam.bai fgbio/$$(*).regrouped.bai")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call merge-alignments,$(sample))))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
