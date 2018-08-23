# Author: David Brown <ndbronw6@gmail.com>

include modules/Makefile.inc

LOGDIR ?= log/fix_bam.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: fix_bam

VPATH = fixed_bam unprocessed_bam

fixed_bam : $(foreach sample,$(SAMPLES),unprocessed_bam/$(sample).bam)

define fix-bam
fixed_bam/%.bam : %.bam
	$$(call RUN,-c -n 1 -s 1G -m 2G,"echo $$<")
  #java -Xmx8G -jar picard.jar RevertSam I=../data/${ID}.bam O=../res/${ID}.ubam SANITIZE=true MAX_DISCARD_FRACTION=0.005 ATTRIBUTE_TO_CLEAR=XT ATTRIBUTE_TO_CLEAR=XN ATTRIBUTE_TO_CLEAR=AS ATTRIBUTE_TO_CLEAR=OC ATTRIBUTE_TO_CLEAR=OP SORT_ORDER=queryname RESTORE_ORIGINAL_QUALITIES=true REMOVE_DUPLICATE_INFORMATION=true REMOVE_ALIGNMENT_INFORMATION=true
  #java -jar -Xmx8G -jar picard.jar MergeBamAlignment R=/ifs/depot/assemblies/H.sapiens/b37_dmp/b37.fasta UNMAPPED_BAM=../res/${ID}.ubam ALIGNED_BAM=../data/${ID}.bam O=../res/${ID}.fixed.bam CREATE_INDEX=true ADD_MATE_CIGAR=true CLIP_ADAPTERS=true CLIP_OVERLAPPING_READS=true INCLUDE_SECONDARY_ALIGNMENTS=false MAX_INSERTIONS_OR_DELETIONS=-1
  #java -jar -Xmx8G picard.jar MarkDuplicates I=../res/${ID}.fixed.bam O=../res/${ID}.dedup.bam M=../res/${ID}.txt
  #java -jar picard.jar AddOrReplaceReadGroups I=../res/${ID}.dedup.bam O=../res/${ID}.bam RGID=${ID} RGLB=${LB} RGPL=${PL} RGPU=${PU} RGSM=${SM}
  #samtools index ../res/${ID}.bam
  #cp ../res/${ID}.bam.bai ../res/${ID}.bai
 endef
 $(foreach pair,$(SAMPLES),\
		$(eval $(call fix-bam,$(sample)))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
