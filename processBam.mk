# various bam processing steps
##### MAKE INCLUDES #####
SAMPLE_FILE ?= samples.txt
SAMPLES ?= $(shell cat $(SAMPLE_FILE))

SPLIT_CHR ?= true

# not primary alignment
# read fails platform/vendor quality checks
BAM_FILTER_FLAGS ?= 768

.DELETE_ON_ERROR:

.SECONDARY: 
#.PHONY: all sort


#BAM_FILES = $(foreach sample,$(SAMPLES),processed_$(sample).bam)
#BAM_FILES = $(foreach sample,$(SAMPLES),$(sample).sorted.reord.rg.rmdup.recal.bam)

#all : $(BAM_FILES) $(addsuffix .bai,$(BAM_FILES))

#%.bam : %.sorted.reord.filtered.fixmate.markdup.bam
#	$(MKDIR) $(@D); ln -v $< $@

# coordinate sorting
sort : $(foreach sample,$(SAMPLES),$(sample).sorted.bam)

# indices
# if bam file is a symlink, need to create a symlink to index
%.bam.bai : %.bam
	$(call LAUNCH_MEM,4G,8G)  $(SAMTOOLS) index $<

%.bai : %.bam
	$(call LAUNCH_MEM,4G,8G) $(SAMTOOLS) index $< $@

# sam to bam
#%.bam : %.sam
#	$(call INIT_MEM,2G,3G) $(SAMTOOLS) view -bS $< > $@

# filter
%.filtered.bam : %.bam
	$(call LAUNCH_MEM,6G,7G) $(SAMTOOLS) view -bF $(BAM_FILTER_FLAGS) $< > $@ 2> $(LOG) && $(RM) $<

#%.fixmate.bam : %.bam
#	$(call LAUNCH_MEM,9G,14G) $(call FIX_MATE_MEM,8G) I=$< O=$@ &> $(LOG) && $(RM) $<

# recalibrate base quality
%.recal_report.grp : %.bam %.bai
	 $(call LAUNCH_MEM,11G,15G) $(call GATK_MEM,10G) -T BaseRecalibrator -R $(REF_FASTA) -knownSites $(DBSNP) -I $< -o $@ &> $(LOGDIR)/$@.log

# recalibration
%.recal.bam : %.bam %.recal_report.grp
	$(call LAUNCH_MEM,11G,15G) $(call GATK_MEM,10G) -T PrintReads -R $(REF_FASTA) -I $< -BQSR $(word 2,$^) -o $@ &> $(LOGDIR)/$@.log && $(RM) $<

# sort only if necessary
#%.sorted.bam : %.bam
#	$(call INIT_MEM,20G,25G) if ! $(SAMTOOLS) view -H $< | grep -q 'SO:coordinate' -; then $(call SORT_SAM_MEM,19G) I=$< O=$@ SO=coordinate &> $(LOG); else ln -v $< $@; fi && $(RM) $<


%.sorted.bam : %.bam
	$(call LAUNCH_MEM,9G,14G) $(call SORT_SAM_MEM,8G) I=$< O=$@ SO=coordinate &> $(LOG) && $(RM) $<

# reorder only if necessary
#%.reord.bam : %.bam
#	if ! $(SAMTOOLS) view -H $< | grep '@SQ' | cut -f 2 | sed 's/^SN://' | diff -q - $(REF_CHR_ORD) > /dev/null; then $(call LAUNCH_MEM,6G,7G) $(call REORDER_SAM_MEM,6G) I=$< O=$@ R=$(REF_FASTA) &> $(LOG); else ln -v $< $@; fi && $(RM) $<


# mark duplicates
%.markdup.bam : %.bam
	$(call LAUNCH_MEM,14G,18G) $(call MARK_DUP_MEM,14G) I=$< O=$@ METRICS_FILE=metrics/$(call strip-suffix,$(@F)).dup_metrics.txt &> $(LOG) && $(RM) $<

%.rmdup.bam : %.bam
	$(call LAUNCH_MEM,4G,7G) $(SAMTOOLS) rmdup $< $@ &> $(LOG) && $(RM) $<

# clean sam files
#%.clean.bam : %.bam
#	$(call LAUNCH_MEM,6G,12G) $(call CLEANBAM_MEM,6G) I=$< O=$@ &> $(LOG)

# add rg
#%.rg.bam : %.bam
#	$(call LAUNCH_MEM,10G,12G) $(call ADD_RG_MEM,10G) I=$< O=$@ RGLB=$(call strip-suffix,$(@F)) RGPL=illumina RGPU=00000000 RGSM=$(call strip-suffix,$(@F)) RGID=$(call strip-suffix,$(@F)) &> $(LOG) && $(RM) $<


# if SPLIT_CHR is set to true, we will split realn processing by chromosome
ifeq ($(SPLIT_CHR),true)
# indel realignment intervals (i.e. where to do MSA)
# split by samples and chromosomes
# %=sample
# $(eval $(call chr-target-aln,chromosome))
define chr-target-realn
%.$1.chr_split.intervals : %.bam %.bam.bai
	$$(call LAUNCH_PARALLEL_MEM,6,2G,2.5G) $$(call GATK_MEM,8G) -T RealignerTargetCreator \
	-I $$< \
	-L $1 \
	-nt 6 -R $$(REF_FASTA)  -o $$@  --known $$(KNOWN_INDELS) \
	&> $$(LOG)
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call chr-target-realn,$(chr))))

# indel realignment per chromosome
# only realign if intervals is non-empty
# %=sample
# $(eval $(call chr-aln,chromosome))
define chr-realn
%.$(1).chr_realn.bam : %.bam %.$(1).chr_split.intervals %.bam.bai
	if [[ -s $$(word 2,$$^) ]]; then $$(call LAUNCH_MEM,9G,12G) $$(call GATK_MEM,8G) -T IndelRealigner \
	-I $$< -R $$(REF_FASTA) -L $1 -targetIntervals $$(word 2,$$^) \
	-o $$@ --knownAlleles $$(KNOWN_INDELS) &> $$(LOG); \
	else $$(call LAUNCH_MEM,9G,12G) $$(call GATK_MEM,8G) -T PrintReads -R $$(REF_FASTA) -I $$< -L $1 -o $$*.$1.chr_realn.bam &> $$(LOG); fi
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call chr-realn,$(chr))))

# merge sample chromosome bams
%.realn.bam : $(foreach chr,$(CHROMOSOMES),%.$(chr).chr_realn.bam) $(foreach chr,$(CHROMOSOMES),%.$(chr).chr_realn.bai)
	$(call LAUNCH_PARALLEL_MEM,2,10G,11G) $(MERGE_SAMS) $(foreach i,$(filter %.bam,$^), I=$i) SORT_ORDER=coordinate O=$@ USE_THREADING=true &> $(LOGDIR)/$@.log && $(RM) $^ 

else # no splitting by chr
%.realn.bam : %.bam %.intervals %.bam.bai
	if [[ -s $(word 2,$^) ]]; then $(call LAUNCH_MEM,9G,12G) $(call GATK_MEM,8G) -T IndelRealigner \
	-I $< -R $(REF_FASTA) -targetIntervals $(word 2,$^) \
	-o $@ --knownAlleles $(KNOWN_INDELS) &> $(LOG) && $(RM) $< ; \
	else mv $< $@ && mv $(word 3,$^) $*.realn.bai &> $(LOG) ; fi

%.intervals : %.bam %.bam.bai
	$(call LAUNCH_PARALLEL_MEM,6,2G,3G) $(call GATK_MEM,8G) -T RealignerTargetCreator \
	-I $< \
	-nt 6 -R $(REF_FASTA)  -o $@  --known $(KNOWN_INDELS) \
	&> $(LOG)
endif

