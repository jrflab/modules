# various bam processing steps
##### MAKE INCLUDES #####
include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc

SAMPLE_FILE ?= samples.txt
SAMPLES ?= $(shell cat $(SAMPLE_FILE))


#LOGDIR = bam/log

# not primary alignment
# read fails platform/vendor quality checks
BAM_FILTER_FLAGS ?= 768

.DELETE_ON_ERROR:

.SECONDARY: 
#.PHONY: all sort


#BAM_FILES = $(foreach sample,$(SAMPLES),processed_bam/$(sample).bam)
#BAM_FILES = $(foreach sample,$(SAMPLES),bam/$(sample).sorted.reord.rg.rmdup.recal.bam)

#all : $(BAM_FILES) $(addsuffix .bai,$(BAM_FILES))

#bam/%.bam : bam/%.sorted.reord.filtered.fixmate.markdup.bam
#	$(MKDIR) $(@D); ln -v $< $@

# coordinate sorting
sort : $(foreach sample,$(SAMPLES),bam/$(sample).sorted.bam)

# indices
# if bam file is a symlink, need to create a symlink to index
%.bam.bai : %.bam
	$(call INIT_MEM,4G,8G)  $(SAMTOOLS) index $<

#&& if readlink $< > /dev/null; then ln -fs $$(readlink -f $<).bai $(<D)/$(@F); fi &> $(LOGDIR)/$(@F).log

# sam to bam
#%.bam : %.sam
#	$(call INIT_MEM,2G,3G) $(SAMTOOLS) view -bS $< > $@

# filter
%.filtered.bam : %.bam
	$(call INIT_MEM,6G,7G) $(SAMTOOLS) view -bF $(BAM_FILTER_FLAGS) $< > $@ 2> $(LOG) && $(RM) $<

%.fixmate.bam : %.bam
	$(call INIT_MEM,9G,14G) $(call FIX_MATE_MEM,8G) I=$< O=$@ &> $(LOG) && $(RM) $<

# recalibrate base quality
%.recal_report.grp : %.bam %.bam.bai
	$(call INIT_MEM,11G,15G) $(call GATK_MEM,10G) -T BaseRecalibrator -R $(REF_FASTA) -knownSites $(DBSNP) -I $< -o $@ &> $(LOGDIR)/$@.log

# recalibration
%.recal.bam : %.bam %.recal_report.grp
	$(call INIT_MEM,11G,15G) $(call GATK_MEM,10G) -T PrintReads -R $(REF_FASTA) -I $< -BQSR $(word 2,$^) -o $@ &> $(LOGDIR)/$@.log && $(RM) $<

# sort only if necessary
#%.sorted.bam : %.bam
#	$(call INIT_MEM,20G,25G) if ! $(SAMTOOLS) view -H $< | grep -q 'SO:coordinate' -; then $(call SORT_SAM_MEM,19G) I=$< O=$@ SO=coordinate &> $(LOG); else ln -v $< $@; fi && $(RM) $<


%.sorted.bam : %.bam
	$(call INIT_MEM,9G,14G) $(call SORT_SAM_MEM,8G) I=$< O=$@ SO=coordinate &> $(LOG) && $(RM) $<

# reorder only if necessary
%.reord.bam : %.bam
	$(call INIT_MEM,6G,7G) if ! $(SAMTOOLS) view -H $< | grep '@SQ' | cut -f 2 | sed 's/^SN://' | diff -q - $(REF_CHR_ORD) > /dev/null; then $(call REORDER_SAM_MEM,6G) I=$< O=$@ R=$(REF_FASTA) &> $(LOG); else ln -v $< $@; fi && $(RM) $<

# mark duplicates
%.markdup.bam : %.bam
	$(call INIT_MEM,14G,18G) $(call MARK_DUP_MEM,14G) I=$< O=$@ METRICS_FILE=metrics/$(call strip-suffix,$(@F)).dup_metrics.txt &> $(LOG) && $(RM) $<

%.rmdup.bam : %.bam
	$(call INIT_MEM,4G,7G) $(SAMTOOLS) rmdup $< $@ &> $(LOG) && $(RM) $<

# clean sam files
%.clean.bam : %.bam
	$(call INIT_MEM,6G,12G) $(call CLEANBAM_MEM,6G) I=$< O=$@ &> $(LOG)
	
# add rg
%.rg.bam : %.bam
	$(call INIT_MEM,10G,12G) $(call ADD_RG_MEM,10G) I=$< O=$@ RGLB=$(call strip-suffix,$(@F)) RGPL=illumina RGPU=00000000 RGSM=$(call strip-suffix,$(@F)) RGID=$(call strip-suffix,$(@F)) &> $(LOG) && $(RM) $<
