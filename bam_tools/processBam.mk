# various bam processing steps
# can be used to reprocess bam files, merge them, or merge and reprocess bam files
# possible post-processing steps are defined in modules/aligners/align.inc
##### MAKE INCLUDES #####

ifndef PROCESS_BAM_MK

include modules/Makefile.inc
include modules/variant_callers/gatk.inc
include modules/aligners/align.inc


LOGDIR ?= log/process_bam.$(NOW)

MERGE_SPLIT_BAMS ?= false  # merge processed split bams
NUM_SORT_SPLITS ?= 50
SORT_SPLIT_SEQ = $(shell seq 0 $$(($(NUM_SORT_SPLITS) - 1)))

CHR1_BASE_RECAL ?= false
BASE_RECAL_OPTS = -knownSites $(DBSNP) $(if $(findstring true,$(CHR1_BASE_RECAL)),-L $(word 1,$(CHROMOSOMES)))

ifneq ($(KNOWN_INDELS),)
REALN_OPTS = --knownAlleles $(KNOWN_INDELS)
REALN_TARGET_OPTS = --known $(KNOWN_INDELS)
endif

# not primary alignment
# read fails platform/vendor quality checks
BAM_FILTER_FLAGS ?= 768

.DELETE_ON_ERROR:
.SECONDARY: 

REPROCESS ?= false

BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)
ifeq ($(REPROCESS),true)
processed_bams : $(BAMS) $(addsuffix .bai,$(BAMS))
bam/%.bam : unprocessed_bam/%.$(BAM_SUFFIX)
	$(INIT) ln -f $< $@
else
ifeq ($(MERGE_SPLIT_BAMS),true)
merged_bams : $(BAMS) $(addsuffix .bai,$(BAMS))
bam/%.bam : unprocessed_bam/%$(if $(findstring true,$(FIX_RG)),.rg).bam
	$(INIT) ln -f $< $@
endif
endif


ifeq ($(MERGE_SPLIT_BAMS),true)
define bam-header
unprocessed_bam/$1.header.sam : $$(foreach split,$2,unprocessed_bam/$$(split).bam)
	$$(INIT) $$(SAMTOOLS) view -H $$< | grep -v '^@RG' > $$@.tmp; \
	for bam in $$(^M); do $$(SAMTOOLS) view -H $$$$bam | grep '^@RG' >> $$@.tmp; done; \
	uniq $$@.tmp > $$@ && $(RM) $$@.tmp
endef
$(foreach sample,$(SPLIT_SAMPLES),$(eval $(call bam-header,$(sample),$(split_lookup.$(sample)))))

define merged-bam
unprocessed_bam/$1.bam : unprocessed_bam/$1.header.sam $$(foreach split,$2,unprocessed_bam/$$(split).bam)
	$$(call LSCRIPT_MEM,12G,15G,"$$(SAMTOOLS) merge -f -h $$< $$@ $$(filter %.bam,$$^)")
endef
$(foreach sample,$(SPLIT_SAMPLES),$(eval $(call merged-bam,$(sample),$(split_lookup.$(sample)))))
endif


# indices
# if bam file is a symlink, need to create a symlink to index
%.bam.bai : %.bam
	$(call LSCRIPT_CHECK_MEM,4G,8G,"$(SAMTOOLS) index $<")

%.bai : %.bam.bai
	$(INIT) cp $< $@

# limit coverage
%.dcov.bam : %.bam
	$(call LSCRIPT_MEM,18G,24G,"$(call GATK_MEM,18G) -T PrintReads -R $(REF_FASTA) -I $< -dcov 50 -o $@")

# filter
%.filtered.bam : %.bam
	$(call LSCRIPT_MEM,6G,7G,"$(SAMTOOLS) view -bF $(BAM_FILTER_FLAGS) $< > $@ && $(RM) $<")

%.fixmate.bam : %.bam
	$(call LSCRIPT_MEM,9G,14G,"$(call FIX_MATE_MEM,8G) I=$< O=$@ && $(RM) $<")

# recalibrate base quality
%.recal_report.grp : %.bam %.bai
	$(call LSCRIPT_MEM,11G,15G,"$(call GATK_MEM,10G) -T BaseRecalibrator -R $(REF_FASTA) $(BASE_RECAL_OPTS) -I $< -o $@")

#sort only if necessary
#%.sorted.bam : %.bam
#	 if ! $(SAMTOOLS) view -H $< | grep -q 'SO:coordinate' -; then $(call LSCRIPT_MEM,20G,25G,"$(call SORT_SAM_MEM,19G) I=$< O=$@ SO=coordinate"); else cp $< $@ && ln -v $< $@; fi && $(RM) $<


# reorder only if necessary
#%.reord.bam : %.bam
#	if ! $(SAMTOOLS) view -H $< | grep '@SQ' | cut -f 2 | sed 's/^SN://' | diff -q - $(REF_CHR_ORD) > /dev/null; then $(call LAUNCH_MEM,6G,7G) $(call REORDER_SAM_MEM,6G) I=$< O=$@ R=$(REF_FASTA) &> $(LOG); else ln -v $< $@; fi && $(RM) $<


# mark duplicates
%.markdup.bam : %.bam
	$(call LSCRIPT_MEM,14G,18G,"$(MKDIR) metrics; $(call MARK_DUP_MEM,14G) I=$< O=$@ METRICS_FILE=metrics/$(call strip-suffix,$(@F)).dup_metrics.txt && $(RM) $<")

%.rmdup.bam : %.bam
	$(call LSCRIPT_MEM,4G,7G,"$(SAMTOOLS) rmdup $< $@ && $(RM) $<")

# clean sam files
%.clean.bam : %.bam
	$(call LSCRIPT_MEM,6G,12G,"$(call CLEANBAM_MEM,6G) I=$< O=$@")

# add rg
%.rg.bam : %.bam
	$(call LSCRIPT_MEM,12G,16G,"$(call ADD_RG_MEM,10G) I=$< O=$@ RGLB=$(call strip-suffix,$(@F)) RGPL=$(SEQ_PLATFORM) RGPU=00000000 RGSM=$(call strip-suffix,$(@F)) RGID=$(call strip-suffix,$(@F)) VERBOSITY=ERROR && $(RM) $<")

# if SPLIT_CHR is set to true, we will split realn processing by chromosome
ifeq ($(SPLIT_CHR),true)
# indel realignment intervals (i.e. where to do MSA)
# split by samples and chromosomes
# %=sample
# $(eval $(call chr-target-aln,chromosome))
define chr-target-realn
%.$1.chr_split.intervals : %.bam %.bam.bai
	$$(call LSCRIPT_PARALLEL_MEM,4,3G,4G,"$$(call GATK_MEM,11G) -T RealignerTargetCreator \
		-I $$(<) \
		-L $1 \
		-nt 4 -R $$(REF_FASTA)  -o $$@ $$(REALN_TARGET_OPTS)")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call chr-target-realn,$(chr))))

# indel realignment per chromosome
# only realign if intervals is non-empty
# %=sample
# $(eval $(call chr-aln,chromosome))
define chr-realn
%.$(1).chr_realn.bam : %.bam %.$(1).chr_split.intervals %.bam.bai
	$$(call LSCRIPT_MEM,9G,12G,"if [[ -s $$(word 2,$$^) ]]; then $$(call GATK_MEM,8G) -T IndelRealigner \
	-I $$(<) -R $$(REF_FASTA) -L $1 -targetIntervals $$(word 2,$$^) \
	-o $$(@) $$(REALN_OPTS); \
	else $$(call GATK_MEM,8G) -T PrintReads -R $$(REF_FASTA) -I $$< -L $1 -o $$@ ; fi")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call chr-realn,$(chr))))

# merge sample realn chromosome bams
%.realn.bam : $(foreach chr,$(CHROMOSOMES),%.$(chr).chr_realn.bam) $(foreach chr,$(CHROMOSOMES),%.$(chr).chr_realn.bai)
	$(call LSCRIPT_PARALLEL_MEM,2,10G,11G,"$(MERGE_SAMS) $(foreach i,$(filter %.bam,$^), I=$(i)) SORT_ORDER=coordinate O=$@ USE_THREADING=true && $(RM) $^ $(@:.realn.bam=.bam)")

# merge sample recal chromosome bams
%.recal.bam : $(foreach chr,$(CHROMOSOMES),%.$(chr).chr_recal.bam) $(foreach chr,$(CHROMOSOMES),%.$(chr).chr_recal.bai)
	$(call LSCRIPT_PARALLEL_MEM,2,10G,11G,"$(MERGE_SAMS) $(foreach i,$(filter %.bam,$^), I=$(i)) SORT_ORDER=coordinate O=$@ USE_THREADING=true && $(RM) $^ $(@:.recal.bam=.bam)")

define chr-recal
%.$1.chr_recal.bam : %.bam %.recal_report.grp
	$$(call LSCRIPT_MEM,11G,15G,"$$(call GATK_MEM,10G) -T PrintReads -L $1 -R $$(REF_FASTA) -I $$< -BQSR $$(<<) -o $$@")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call chr-recal,$(chr))))

else # no splitting by chr

# recalibration
%.recal.bam : %.bam %.recal_report.grp
	$(call LSCRIPT_MEM,11G,15G,"$(call GATK_MEM,10G) -T PrintReads -R $(REF_FASTA) -I $< -BQSR $(word 2,$^) -o $@ && $(RM) $<")

%.realn.bam : %.bam %.intervals %.bam.bai
	if [[ -s $(word 2,$^) ]]; then $(call LSCRIPT_MEM,9G,12G,"$(call GATK_MEM,8G) -T IndelRealigner \
	-I $< -R $(REF_FASTA) -targetIntervals $(<<) \
	-o $@ $(REALN_OPTS) && $(RM) $<") ; \
	else mv $< $@ ; fi

%.intervals : %.bam %.bam.bai
	$(call LSCRIPT_PARALLEL_MEM,4,2.5G,3G,"$(call GATK_MEM,8G) -T RealignerTargetCreator \
	-I $< \
	-nt 4 -R $(REF_FASTA) -o $@ $(REALN_TARGET_OPTS)")
endif

endif
PROCESS_BAM_MK = true


