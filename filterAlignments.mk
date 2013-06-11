# This module is for filtering sequencing alignments
# Author: Fong Chun Chan <fongchunchan@gmail.com>

#FILTERED_BAM_FILES=$(foreach sample,$(SAMPLES),tophat/filtered_bam/$(sample).tophat.filtered.rmdup.bam)

%.filtered.rmdup.bam : %.bam
	SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,1G,2G)" \
	$(MKDIR) $(@D)/logs;\
	${SAMTOOLS} view -F 1796 -b $< 2> $(@D)/logs/$(@F).filter.log | ${SAMTOOLS} rmdup - $@ 2> $(@D)/logs/$(@F).rmdup.log

%.filtered.bam : %.bam
	SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,1G,2G)" \
	$(MKDIR) $(@D)/logs;\
	${SAMTOOLS} view -F 1796 -b $< 2> $(@D)/logs/$(@F).filter.log | $(@D)/logs/$(@F).rmdup.log

%.rmdup.bam : %.bam
	SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,1G,2G)" \
	$(MKDIR) $(@D)/logs;\
	${SAMTOOLS} rmdup $< $@ 2> $(@D)/logs/$(@F).rmdup.log

%.bam.bai : %.bam
	SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,1G,2G)" \
	$(MKDIR) $(@D)/logs;\
	${SAMTOOLS} index $< $@ 2> $(@D)/logs/$(@F).index.log
